#! /usr/bin/env python

import json

import pandas as pd
import q2_PSEA

from q2_PSEA.actions.psea import (
    create_fgsea_table_for_pair,
    make_psea_table,
    run_iterative_peptide_analysis,
    run_iterative_process_single_pair,
)
from q2_PSEA.formats import (
    IterativePairStateDirFmt,
    IterativePairStateFormat,
    IterativePeptideSetsDirFmt,
    IterativePeptideSetsFormat,
    PSEAPairsDirFmt,
    PSEAPairsFormat,
    PSEATableDirFmt,
    PSEATableFormat,
    PeptideSetsDirFmt,
    PeptideSetsFormat,
    SpeciesColorsDirFmt,
    SpeciesColorsFormat,
    SpeciesTaxonomyDirFmt,
    SpeciesTaxonomyFormat,
    SplineDataDirFmt,
    SplineDataFormat,
)
from q2_PSEA.types import (
    IterativePairState,
    IterativePeptideSets,
    PSEAPairs,
    PSEATable,
    PeptideSets,
    SpeciesColors,
    SpeciesTaxonomy,
    SplineData,
)
from q2_pepsirf.format_types import FeatureTable, Zscore
from qiime2.plugin import Bool, Float, Int, Plugin, Str, Visualization


plugin = Plugin(
    "psea",
    version=q2_PSEA.__version__,
    website="https://github.com/LadnerLab/q2-PSEA.git",
    description="Qiime2 Plugin for PSEA."
)


plugin.register_semantic_types(
    PSEAPairs,
    PeptideSets,
    SpeciesTaxonomy,
    SpeciesColors,
    PSEATable,
    SplineData,
    IterativePairState,
    IterativePeptideSets,
)

plugin.register_formats(
    PSEAPairsFormat,
    PeptideSetsFormat,
    SpeciesTaxonomyFormat,
    SpeciesColorsFormat,
    PSEATableFormat,
    SplineDataFormat,
    IterativePairStateFormat,
    IterativePeptideSetsFormat,
    PSEAPairsDirFmt,
    PeptideSetsDirFmt,
    SpeciesTaxonomyDirFmt,
    SpeciesColorsDirFmt,
    PSEATableDirFmt,
    SplineDataDirFmt,
    IterativePairStateDirFmt,
    IterativePeptideSetsDirFmt,
)

plugin.register_artifact_class(PSEAPairs, PSEAPairsDirFmt)
plugin.register_artifact_class(PeptideSets, PeptideSetsDirFmt)
plugin.register_artifact_class(SpeciesTaxonomy, SpeciesTaxonomyDirFmt)
plugin.register_artifact_class(SpeciesColors, SpeciesColorsDirFmt)
plugin.register_artifact_class(PSEATable, PSEATableDirFmt)
plugin.register_artifact_class(SplineData, SplineDataDirFmt)
plugin.register_artifact_class(IterativePairState, IterativePairStateDirFmt)
plugin.register_artifact_class(IterativePeptideSets, IterativePeptideSetsDirFmt)


@plugin.register_transformer
def _psea_table_df_to_format(data: pd.DataFrame) -> PSEATableFormat:
    ff = PSEATableFormat()
    data.to_csv(str(ff), sep="\t", index=False)
    return ff


@plugin.register_transformer
def _psea_table_format_to_df(ff: PSEATableFormat) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep="\t")


@plugin.register_transformer
def _spline_data_df_to_format(data: pd.DataFrame) -> SplineDataFormat:
    ff = SplineDataFormat()
    data.to_csv(str(ff), sep="\t", index=False)
    return ff


@plugin.register_transformer
def _spline_data_format_to_df(ff: SplineDataFormat) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep="\t")


@plugin.register_transformer
def _iterative_pair_state_dict_to_format(data: dict) -> IterativePairStateFormat:
    ff = IterativePairStateFormat()
    with ff.open() as fh:
        json.dump(data, fh)
    return ff


@plugin.register_transformer
def _iterative_pair_state_format_to_dict(ff: IterativePairStateFormat) -> dict:
    with ff.open() as fh:
        return json.load(fh)


@plugin.register_transformer
def _iterative_peptide_sets_dict_to_format(data: dict) -> IterativePeptideSetsFormat:
    ff = IterativePeptideSetsFormat()
    with ff.open() as fh:
        json.dump(data, fh)
    return ff


@plugin.register_transformer
def _iterative_peptide_sets_format_to_dict(ff: IterativePeptideSetsFormat) -> dict:
    with ff.open() as fh:
        return json.load(fh)


plugin.methods.register_function(
    function=create_fgsea_table_for_pair,
    inputs={
        "processed_scores": FeatureTable[Zscore],
        "peptide_sets": PeptideSets,
        "species_taxa": SpeciesTaxonomy,
    },
    parameters={
        "pair_a": Str,
        "pair_b": Str,
        "threshold": Float,
        "permutation_num": Int,
        "min_size": Int,
        "max_size": Int,
        "spline_type": Str,
        "degree": Int,
        "dof": Int,
        "seed": Int,
    },
    outputs=[("psea_table", PSEATable), ("spline_data", SplineData)],
    name="Create FGSEA Table For Pair",
    description=(
        "Create FGSEA and spline data tables for a single sample pair."
    ),
)


plugin.pipelines.register_function(
    function=run_iterative_process_single_pair,
    inputs={
        "processed_scores": FeatureTable[Zscore],
        "species_taxa": SpeciesTaxonomy,
        "pair_state": IterativePairState,
    },
    parameters={
        "pair_a": Str,
        "pair_b": Str,
        "threshold": Float,
        "permutation_num": Int,
        "min_size": Int,
        "max_size": Int,
        "spline_type": Str,
        "degree": Int,
        "dof": Int,
        "p_val_thresh": Float,
        "nes_thresh": Float,
        "iter_out_dir": Str,
        "seed": Int,
    },
    outputs=[("updated_pair_state", IterativePairState)],
    name="Run Iterative Process Single Pair",
    description=(
        "Run one iterative filtering step for a single sample pair."
    ),
)


plugin.pipelines.register_function(
    function=run_iterative_peptide_analysis,
    inputs={
        "pairs": PSEAPairs,
        "processed_scores": FeatureTable[Zscore],
        "og_peptide_sets": PeptideSets,
        "species_taxa": SpeciesTaxonomy,
    },
    parameters={
        "threshold": Float,
        "permutation_num": Int,
        "min_size": Int,
        "max_size": Int,
        "spline_type": Str,
        "degree": Int,
        "dof": Int,
        "p_val_thresh": Float,
        "nes_thresh": Float,
        "iter_tables_dir": Str,
        "seed": Int,
    },
    outputs=[("pair_peptide_sets", IterativePeptideSets)],
    name="Run Iterative Peptide Analysis",
    description=(
        "Run iterative peptide filtering across all sample pairs."
    ),
)


plugin.pipelines.register_function(
    function=make_psea_table,
    inputs={
        "scores": FeatureTable[Zscore],
        "pairs": PSEAPairs,
        "peptide_sets": PeptideSets,
        "species_taxa": SpeciesTaxonomy,
        "species_colors": SpeciesColors,
    },
    parameters={
        "threshold": Float,
        "p_val_thresh": Float,
        "nes_thresh": Float,
        "min_size": Int,
        "max_size": Int,
        "permutation_num": Int,
        "spline_type": Str,
        "degree": Int,
        "dof": Int,
        "table_dir": Str,
        "pepsirf_binary": Str,
        "iterative_analysis": Bool,
        "iter_tables_dir": Str,
        "max_workers": Int,
        "summary_tables_dir": Str,
        "vis_outputs_dir": Str,
        "seed": Int,
    },
    outputs=[
        ("scatter_plot", Visualization),
        ("volcano_plot", Visualization),
        ("ae_plots", Visualization),
    ],
    name="Make PSEA Table",
    description=(
        "Python wrapper around R functions to perform PSEA and create output plots."
    ),
)
