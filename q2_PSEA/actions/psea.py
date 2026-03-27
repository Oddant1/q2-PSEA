import os
import tempfile
import time

import numpy as np
import pandas as pd
import q2_PSEA.actions.splines as splines
import q2_PSEA.utils as utils
import rpy2.robjects as ro

from math import log, pow
from q2_PSEA.actions.r_functions import INTERNAL
from q2_PSEA.formats import (
    PSEAPairsFormat,
    PeptideSetsFormat,
    SpeciesColorsFormat,
    SpeciesTaxonomyFormat,
)
from q2_pepsirf.format_types import PepsirfContingencyTSVFormat
from rpy2.robjects import pandas2ri


pandas2ri.activate()


def _pair_key(pair):
    return f"{pair[0]}~{pair[1]}"


def _load_pairs(pairs_file):
    pairs = []
    pair_2_title = {}
    with open(pairs_file, "r") as fh:
        fh.readline()
        for line in fh.readlines():
            line_tup = tuple(line.replace("\n", "").split("\t"))
            pair = line_tup[0:2]
            pairs.append(pair)

            if len(line_tup) > 2:
                pair_2_title[pair] = line_tup[2]
            else:
                pair_2_title[pair] = ""

    return pairs, pair_2_title


def _read_gmt_as_dict(gmt_file):
    gmt_dict = {}
    with open(gmt_file, "r") as file:
        for line in file.readlines():
            if not line.strip():
                continue

            row = line.strip().split("\t")
            if len(row) < 3:
                raise AssertionError(
                    "Iterative analysis requires GMT-formatted peptide sets."
                )

            gmt_dict[row[0]] = set([pep for pep in row[2:] if pep])

    return gmt_dict


def _serialize_gmt_dict(gmt_dict):
    return {species: sorted(peptides) for species, peptides in gmt_dict.items()}


def _deserialize_gmt_dict(gmt_dict):
    return {species: set(peptides) for species, peptides in gmt_dict.items()}


def _matrix_df_from_event_dict(event_dict, pairs):
    pair_index = [_pair_key(pair) for pair in pairs]

    if not event_dict:
        return pd.DataFrame(index=pair_index)

    matrix_df = pd.DataFrame.from_dict(event_dict)
    matrix_df.index = pair_index
    matrix_df.sort_index(inplace=True)

    if len(matrix_df.columns) > 0:
        matrix_df = matrix_df[
            sorted(
                matrix_df.columns.tolist(),
                key=lambda col: matrix_df[col].sum(),
                reverse=True
            )
        ]

    return matrix_df


def make_psea_table(
        ctx,
        scores,
        pairs,
        peptide_sets,
        species_taxa,
        species_colors,
        threshold,
        p_val_thresh=0.05,
        nes_thresh=1,
        min_size=15,
        max_size=2000,
        permutation_num=10000,
        spline_type="r-smooth",
        degree=3,
        dof=0,
        table_dir="./psea_table_outdir",
        pepsirf_binary="pepsirf",
        iterative_analysis=True,
        iter_tables_dir="",
        max_workers=None,
        summary_tables_dir="./psea_ae_summary_tables",
        vis_outputs_dir=None,
        seed=149
):
    del pepsirf_binary
    del max_workers

    start_time = time.perf_counter()

    volcano = ctx.get_action("ps-plot", "volcano")
    zscatter = ctx.get_action("ps-plot", "zscatter")
    aeplots = ctx.get_action("ps-plot", "aeplots")

    create_fgsea = ctx.get_action("psea", "create_fgsea_table_for_pair")
    run_iterative = ctx.get_action("psea", "run_iterative_peptide_analysis")

    assert spline_type in splines.SPLINE_TYPES, \
        f"'{spline_type}' is not a valid spline method!"
    assert not os.path.exists(table_dir), \
        f"'{table_dir}' already exists! Please move or remove this directory."
    assert not os.path.exists(summary_tables_dir), \
        f"'{summary_tables_dir}' already exists! Please move or remove this directory."
    if vis_outputs_dir is not None:
        assert not os.path.exists(vis_outputs_dir), \
            f"'{vis_outputs_dir}' already exists! Please move or remove this directory."
        os.mkdir(vis_outputs_dir)

    if iterative_analysis:
        try:
            _read_gmt_as_dict(str(peptide_sets.view(PeptideSetsFormat)))
        except AssertionError:
            raise
        except Exception as e:
            raise AssertionError(
                "You are running iterative analysis without a GMT peptide sets file."
            ) from e

    os.mkdir(table_dir)
    os.mkdir(summary_tables_dir)

    pairs_file = str(pairs.view(PSEAPairsFormat))
    species_color_file = str(species_colors.view(SpeciesColorsFormat))
    pairs_list, _ = _load_pairs(pairs_file)

    scores_file = str(scores.view(PepsirfContingencyTSVFormat))
    scores_df = pd.read_csv(scores_file, sep="\t", index_col=0)
    processed_scores = process_scores(scores_df, pairs_list)

    with tempfile.TemporaryDirectory() as tempdir:
        processed_scores_file = os.path.join(tempdir, "transformed_scores.tsv")
        processed_scores.to_csv(processed_scores_file, sep="\t")

        processed_scores_art = ctx.make_artifact(
            type="FeatureTable[Zscore]",
            view=processed_scores_file,
            view_type=PepsirfContingencyTSVFormat
        )

        if iterative_analysis:
            pair_peptide_sets_art, = run_iterative(
                pairs=pairs,
                processed_scores=processed_scores_art,
                og_peptide_sets=peptide_sets,
                species_taxa=species_taxa,
                threshold=threshold,
                permutation_num=permutation_num,
                min_size=min_size,
                max_size=max_size,
                spline_type=spline_type,
                degree=degree,
                dof=dof,
                p_val_thresh=p_val_thresh,
                nes_thresh=nes_thresh,
                iter_tables_dir=iter_tables_dir,
                seed=seed
            )
            pair_peptide_sets_dict = pair_peptide_sets_art.view(dict)
        else:
            pair_peptide_sets_dict = {}

        pos_nes_event_matrix = {}
        neg_nes_event_matrix = {}
        zero_nes_event_matrix = {}
        empty_pair_row = [0] * len(pairs_list)
        pos_nes_count_dict = {}
        neg_nes_count_dict = {}
        zero_nes_count_dict = {}

        pair_spline_frames = []
        taxa_access = "species_name"

        for pair_idx, pair in enumerate(pairs_list):
            table_prefix = _pair_key(pair)

            if iterative_analysis:
                pair_gmt_dict = _deserialize_gmt_dict(pair_peptide_sets_dict[table_prefix])
                pair_pep_sets_file = os.path.join(
                    tempdir,
                    f"{pair[0]}_{pair[1]}".replace(".", "-") + ".gmt"
                )
                write_gmt_from_dict(pair_pep_sets_file, pair_gmt_dict)
                pair_peptide_sets_art = ctx.make_artifact(
                    type="PeptideSets",
                    view=pair_pep_sets_file,
                    view_type=PeptideSetsFormat
                )
            else:
                pair_peptide_sets_art = peptide_sets

            psea_table_art, spline_data_art = create_fgsea(
                processed_scores=processed_scores_art,
                peptide_sets=pair_peptide_sets_art,
                species_taxa=species_taxa,
                pair_a=pair[0],
                pair_b=pair[1],
                threshold=threshold,
                permutation_num=permutation_num,
                min_size=min_size,
                max_size=max_size,
                spline_type=spline_type,
                degree=degree,
                dof=dof,
                seed=seed
            )

            table_df = psea_table_art.view(pd.DataFrame)
            if "species_name" not in table_df.columns:
                taxa_access = "ID"

            table_df.to_csv(
                os.path.join(table_dir, f"{table_prefix}_psea_table.tsv"),
                sep="\t",
                index=False
            )

            spline_df = spline_data_art.view(pd.DataFrame)
            pair_spline_frames.append(spline_df)

            for _, row in table_df.iterrows():
                taxa = row[taxa_access]

                if row["p.adjust"] < p_val_thresh and np.absolute(row["NES"]) > nes_thresh:
                    if row["NES"] > 0:
                        if taxa not in pos_nes_event_matrix:
                            pos_nes_event_matrix[taxa] = empty_pair_row.copy()
                        pos_nes_event_matrix[taxa][pair_idx] = 1

                        if taxa not in pos_nes_count_dict:
                            pos_nes_count_dict[taxa] = 0
                        pos_nes_count_dict[taxa] += 1

                    elif row["NES"] < 0:
                        if taxa not in neg_nes_event_matrix:
                            neg_nes_event_matrix[taxa] = empty_pair_row.copy()
                        neg_nes_event_matrix[taxa][pair_idx] = 1

                        if taxa not in neg_nes_count_dict:
                            neg_nes_count_dict[taxa] = 0
                        neg_nes_count_dict[taxa] += 1

                    else:
                        if taxa not in zero_nes_event_matrix:
                            zero_nes_event_matrix[taxa] = empty_pair_row.copy()
                        zero_nes_event_matrix[taxa][pair_idx] = 1

                        if taxa not in zero_nes_count_dict:
                            zero_nes_count_dict[taxa] = 0
                        zero_nes_count_dict[taxa] += 1

        pos_nes_event_matrix_df = _matrix_df_from_event_dict(pos_nes_event_matrix, pairs_list)
        pos_nes_event_matrix_df.to_csv(
            os.path.join(summary_tables_dir, "Positive_NES_taxa_matrix.tsv"),
            sep="\t"
        )

        neg_nes_event_matrix_df = _matrix_df_from_event_dict(neg_nes_event_matrix, pairs_list)
        neg_nes_event_matrix_df.to_csv(
            os.path.join(summary_tables_dir, "Negative_NES_taxa_matrix.tsv"),
            sep="\t"
        )

        pos_nes_count_dict = {
            k: v for k, v in sorted(
                pos_nes_count_dict.items(),
                key=lambda item: item[1],
                reverse=True
            )
        }
        neg_nes_count_dict = {
            k: v for k, v in sorted(
                neg_nes_count_dict.items(),
                key=lambda item: item[1],
                reverse=True
            )
        }

        with open(os.path.join(summary_tables_dir, "Positive_NES_AE.tsv"), "w") as pos_file:
            pos_file.write("Species\tEvents\n")
            for taxa in pos_nes_count_dict.keys():
                pos_file.write(f"{taxa}\t{pos_nes_count_dict[taxa]}\n")

        with open(os.path.join(summary_tables_dir, "Negative_NES_AE.tsv"), "w") as neg_file:
            neg_file.write("Species\tEvents\n")
            for taxa in neg_nes_count_dict.keys():
                neg_file.write(f"{taxa}\t{neg_nes_count_dict[taxa]}\n")

        if len(zero_nes_count_dict) > 0:
            print("\n")
            for taxa in zero_nes_count_dict.keys():
                if zero_nes_count_dict[taxa] > 1:
                    print(f"{zero_nes_count_dict[taxa]} events for {taxa}, which has an NES of 0")
                else:
                    print(f"{zero_nes_count_dict[taxa]} event for {taxa}, which has an NES of 0")
                print("The pairs which this occurred are: ")
                for pair_index in range(len(zero_nes_event_matrix[taxa])):
                    if zero_nes_event_matrix[taxa][pair_index] == 1:
                        print(_pair_key(pairs_list[pair_index]))

        spline_data_file = os.path.join(tempdir, "spline_data.tsv")
        if pair_spline_frames:
            pair_spline_data = pd.concat(pair_spline_frames, axis=0, ignore_index=True)
        else:
            pair_spline_data = pd.DataFrame({"x": [], "y": [], "pair": []})
        pair_spline_data.to_csv(spline_data_file, sep="\t", index=False)

        scatter_plot, = zscatter(
            zscores=processed_scores_art,
            pairs_file=pairs_file,
            spline_file=spline_data_file,
            p_val_access="p.adjust",
            le_peps_access="core_enrichment",
            taxa_access=taxa_access,
            highlight_data=table_dir,
            highlight_threshold=p_val_thresh,
            colors_file=species_color_file,
            vis_outputs_dir=vis_outputs_dir
        )

        volcano_plot, = volcano(
            xy_dir=table_dir,
            xy_access=["NES", "p.adjust"],
            taxa_access=taxa_access,
            x_threshold=nes_thresh,
            y_threshold=p_val_thresh,
            xy_labels=["Enrichment score", "Adjusted p-values"],
            pairs_file=pairs_file,
            colors_file=species_color_file,
            vis_outputs_dir=vis_outputs_dir
        )

        ae_plot, = aeplots(
            pos_nes_ae_file=os.path.join(summary_tables_dir, "Positive_NES_AE.tsv"),
            neg_nes_ae_file=os.path.join(summary_tables_dir, "Negative_NES_AE.tsv"),
            xy_access=["Events", "Species"],
            xy_labels=["Number of AEs in cohort", "Species"],
            colors_file=species_color_file,
            vis_outputs_dir=vis_outputs_dir
        )

    end_time = time.perf_counter()
    print(f"\nFinished in {round(end_time - start_time, 2)} seconds")

    return scatter_plot, volcano_plot, ae_plot


def create_fgsea_table_for_pair(
    processed_scores: PepsirfContingencyTSVFormat,
    peptide_sets: PeptideSetsFormat,
    species_taxa: SpeciesTaxonomyFormat,
    pair_a: str,
    pair_b: str,
    threshold: float,
    permutation_num: int,
    min_size: int,
    max_size: int,
    spline_type: str,
    degree: int,
    dof: int,
    seed: int
) -> (pd.DataFrame, pd.DataFrame):
    print(f"Working on pair ({pair_a}, {pair_b})...")

    processed_scores_df = pd.read_csv(str(processed_scores), sep="\t", index_col=0)
    species_taxa_file = str(species_taxa)

    processed_scores_df, peptide_sets_df = utils.remove_peptides(
        processed_scores_df,
        str(peptide_sets)
    )

    pair = [pair_a, pair_b]
    data_sorted = processed_scores_df.loc[:, pair].sort_values(by=pair_a)
    x = data_sorted.loc[:, pair_a].to_numpy()
    y = data_sorted.loc[:, pair_b].to_numpy()

    if not dof:
        dof = ro.NULL

    if spline_type == "py-smooth":
        yfit = splines.smooth_spline(x, y)
    elif spline_type == "cubic":
        yfit = splines.R_SPLINES.cubic_spline(x, y, degree, dof)
    else:
        yfit = splines.R_SPLINES.smooth_spline(x, y)

    max_z = np.apply_over_axes(np.max, data_sorted.loc[:, pair], 1)
    max_z = pd.Series(
        data=[num for elem in max_z for num in elem],
        index=data_sorted.index
    )
    delta_z = pd.Series(data=y - yfit, index=data_sorted.index)

    table = INTERNAL.psea(
        max_z,
        delta_z,
        peptide_sets_df,
        species_taxa_file,
        threshold,
        permutation_num,
        min_size,
        max_size,
        seed
    )

    with (ro.default_converter + pandas2ri.converter).context():
        table = ro.conversion.get_conversion().rpy2py(table)

    spline_df = pd.DataFrame({
        "x": x.tolist(),
        "y": yfit.tolist(),
        "pair": [_pair_key((pair_a, pair_b))] * len(x)
    })

    return table, spline_df


def process_scores(scores, pairs) -> pd.DataFrame:
    base = 2
    offset = 3
    power = pow(base, offset)

    reps_list = []
    for pair in pairs:
        for rep in pair:
            reps_list.append(rep)

    reps_list = list(np.unique(reps_list))
    processed_scores = scores.loc[:, reps_list]

    processed_scores = processed_scores.apply(lambda row: power + row, axis=0)
    processed_scores = processed_scores.apply(
        lambda row: row.apply(lambda val: 1 if val < 1 else val),
        axis=0
    )

    return processed_scores.apply(
        lambda row: row.apply(lambda val: log(val, base) - offset)
    )


def run_iterative_peptide_analysis(
    ctx,
    pairs,
    processed_scores,
    og_peptide_sets,
    species_taxa,
    threshold,
    permutation_num,
    min_size,
    max_size,
    spline_type,
    degree,
    dof,
    p_val_thresh,
    nes_thresh,
    iter_tables_dir,
    seed
):
    iterative_single_pair = ctx.get_action("psea", "run_iterative_process_single_pair")

    pairs_file = str(pairs.view(PSEAPairsFormat))
    pairs_list, _ = _load_pairs(pairs_file)

    gmt_dict = _read_gmt_as_dict(str(og_peptide_sets.view(PeptideSetsFormat)))

    pair_state_artifacts = {}
    sig_species_found_dict = {}
    for pair in pairs_list:
        initial_state = {
            "tested_species": [],
            "gmt_dict": _serialize_gmt_dict(gmt_dict),
            "sig_species_found": True,
        }
        pair_state_artifacts[pair] = ctx.make_artifact(
            type="IterativePairState",
            view=initial_state,
            view_type=dict
        )
        sig_species_found_dict[pair] = True

    iteration_num = 1
    while any(sig_species_found_dict.values()):
        print(f"\nIteration: {iteration_num}")

        if iter_tables_dir:
            iter_out_dir = os.path.join(iter_tables_dir, f"Iteration_{iteration_num}")
            os.makedirs(iter_out_dir, exist_ok=True)
        else:
            iter_out_dir = ""

        for pair in pairs_list:
            if not sig_species_found_dict[pair]:
                continue

            updated_state_art, = iterative_single_pair(
                processed_scores=processed_scores,
                species_taxa=species_taxa,
                pair_state=pair_state_artifacts[pair],
                pair_a=pair[0],
                pair_b=pair[1],
                threshold=threshold,
                permutation_num=permutation_num,
                min_size=min_size,
                max_size=max_size,
                spline_type=spline_type,
                degree=degree,
                dof=dof,
                p_val_thresh=p_val_thresh,
                nes_thresh=nes_thresh,
                iter_out_dir=iter_out_dir,
                seed=seed
            )

            pair_state_artifacts[pair] = updated_state_art
            updated_state = updated_state_art.view(dict)
            sig_species_found_dict[pair] = updated_state["sig_species_found"]

        iteration_num += 1

    print("\nEnd of Iterative Peptide Analysis\n")

    pair_sets = {}
    for pair in pairs_list:
        state = pair_state_artifacts[pair].view(dict)
        pair_sets[_pair_key(pair)] = state["gmt_dict"]

    pair_peptide_sets = ctx.make_artifact(
        type="IterativePeptideSets",
        view=pair_sets,
        view_type=dict
    )

    return pair_peptide_sets,


def run_iterative_process_single_pair(
    ctx,
    processed_scores,
    species_taxa,
    pair_state,
    pair_a,
    pair_b,
    threshold,
    permutation_num,
    min_size,
    max_size,
    spline_type,
    degree,
    dof,
    p_val_thresh,
    nes_thresh,
    iter_out_dir,
    seed
):
    create_fgsea = ctx.get_action("psea", "create_fgsea_table_for_pair")

    pair_state_dict = pair_state.view(dict)
    tested_species = set([str(s) for s in pair_state_dict["tested_species"]])
    gmt_dict = _deserialize_gmt_dict(pair_state_dict["gmt_dict"])

    with tempfile.TemporaryDirectory() as tempdir:
        pair_sets_filename = os.path.join(
            tempdir,
            f"{pair_a}_{pair_b}".replace(".", "-") + ".gmt"
        )
        write_gmt_from_dict(pair_sets_filename, gmt_dict)

        pair_sets_art = ctx.make_artifact(
            type="PeptideSets",
            view=pair_sets_filename,
            view_type=PeptideSetsFormat
        )

        table_art, _ = create_fgsea(
            processed_scores=processed_scores,
            peptide_sets=pair_sets_art,
            species_taxa=species_taxa,
            pair_a=pair_a,
            pair_b=pair_b,
            threshold=threshold,
            permutation_num=permutation_num,
            min_size=min_size,
            max_size=max_size,
            spline_type=spline_type,
            degree=degree,
            dof=dof,
            seed=seed
        )

        table = table_art.view(pd.DataFrame)

    table.sort_values(by=["p.adjust"], ascending=True, inplace=True)

    if iter_out_dir:
        os.makedirs(iter_out_dir, exist_ok=True)
        table.to_csv(os.path.join(iter_out_dir, f"{_pair_key((pair_a, pair_b))}.tsv"), sep="\t", index=False)

    sig_species_found = False

    for _, row in table.iterrows():
        row_id = str(row["ID"])
        if row["p.adjust"] < p_val_thresh and np.absolute(row["NES"]) > nes_thresh \
                and row_id not in tested_species:

            species_name = row["species_name"] if "species_name" in row else row_id
            print(f"Found {species_name} in ({pair_a}, {pair_b}) to be significant")

            sig_species_found = True
            tested_species.add(row_id)

            all_tested_peps = set(str(row["all_tested_peptides"]).split("/"))
            for gmt_species in gmt_dict.keys():
                if gmt_species != row_id:
                    gmt_dict[gmt_species] = gmt_dict[gmt_species] - all_tested_peps

            break

    updated_state = {
        "tested_species": sorted(tested_species),
        "gmt_dict": _serialize_gmt_dict(gmt_dict),
        "sig_species_found": sig_species_found,
    }

    updated_pair_state = ctx.make_artifact(
        type="IterativePairState",
        view=updated_state,
        view_type=dict
    )

    return updated_pair_state,


def write_gmt_from_dict(outfile_name, gmt_dict) -> None:
    with open(outfile_name, "w") as gmt_file:
        for species in gmt_dict.keys():
            gmt_file.write(f"{species}\t\t")

            for peptide in gmt_dict[species]:
                gmt_file.write(f"{peptide}\t")

            gmt_file.write("\n")
