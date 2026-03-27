#! /usr/bin/env python
import importlib

from q2_types.feature_data import FeatureData
from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import Epitope, GMT, MappedEpitope, PSEAScores, Zscore

import q2_PSEA

from q2_PSEA.actions.psea import (
    create_fgsea_table_for_pair,
    make_psea_table,
    run_iterative_peptide_analysis,
    run_iterative_process_single_pair,
)
from q2_PSEA.actions.splines import SPLINE_TYPES
from q2_PSEA.format_types import (
    PSEAPairs,
    PSEAPairsDirFmt,
    PSEAPairsFormat,
    PSEASpeciesColor,
    PSEASpeciesColorDirFmt,
    PSEASpeciesColorFormat,
    PSEASpeciesTaxa,
    PSEASpeciesTaxaDirFmt,
    PSEASpeciesTaxaFormat,
)
from qiime2.plugin import (
    Bool,
    Choices,
    Collection,
    Float,
    Int,
    List,
    Plugin,
    Str,
    Visualization,
)


plugin = Plugin(
    "psea",
    version=q2_PSEA.__version__,
    website="https://github.com/LadnerLab/q2-PSEA.git",
    description="Qiime2 plugin for peptide set enrichment analysis.",
)

plugin.register_formats(
    PSEAPairsFormat,
    PSEAPairsDirFmt,
    PSEASpeciesTaxaFormat,
    PSEASpeciesTaxaDirFmt,
    PSEASpeciesColorFormat,
    PSEASpeciesColorDirFmt,
)

plugin.register_semantic_types(
    PSEAPairs,
    PSEASpeciesTaxa,
    PSEASpeciesColor,
)

plugin.register_semantic_type_to_format(
    PSEAPairs,
    PSEAPairsDirFmt,
)
plugin.register_semantic_type_to_format(
    PSEASpeciesTaxa,
    PSEASpeciesTaxaDirFmt,
)
plugin.register_semantic_type_to_format(
    PSEASpeciesColor,
    PSEASpeciesColorDirFmt,
)

plugin.methods.register_function(
    function=create_fgsea_table_for_pair,
    inputs={
        "scores": FeatureTable[Zscore],
        "peptide_sets": GMT,
        "species_taxa": PSEASpeciesTaxa,
        "epitope_map": FeatureData[MappedEpitope],
        "mapped_processed_scores": FeatureTable[Zscore],
    },
    parameters={
        "pair": List[Str],
        "threshold": Float,
        "permutation_num": Int,
        "min_size": Int,
        "max_size": Int,
        "spline_type": Str % Choices(SPLINE_TYPES),
        "degree": Int,
        "dof": Int,
        "seed": Int,
    },
    outputs=[("psea_table", FeatureData[PSEAScores])],
    input_descriptions={
        "scores": "Processed Z-score table.",
        "peptide_sets": "GMT mapping taxa IDs to peptide IDs.",
        "species_taxa": "Optional species-name to taxa-ID mapping.",
        "epitope_map": "Optional epitope mapping used for residual collapse.",
        "mapped_processed_scores": (
            "Optional processed epitope-level scores, required when epitope_map is provided."
        ),
    },
    parameter_descriptions={
        "pair": "Two sample IDs to compare.",
        "threshold": "Minimum Z score required to include peptides in GSEA.",
        "permutation_num": "Number of permutations in fgsea.",
        "min_size": "Minimum peptide set size.",
        "max_size": "Maximum peptide set size.",
        "spline_type": "Spline implementation used for residual fitting.",
        "degree": "Degree used by cubic spline mode.",
        "dof": "Degrees of freedom for cubic spline mode.",
        "seed": "Random seed used by fgsea.",
    },
    output_descriptions={
        "psea_table": "PSEA table for the specified pair.",
    },
    name="Create FGSEA Table For Pair",
    description=(
        "Create a single pair-level PSEA table as a FeatureData[PSEAScores] artifact."
    ),
)

plugin.pipelines.register_function(
    function=run_iterative_process_single_pair,
    inputs={
        "scores": FeatureTable[Zscore],
        "peptide_sets": GMT,
        "species_taxa": PSEASpeciesTaxa,
        "epitope_map": FeatureData[MappedEpitope],
        "mapped_processed_scores": FeatureTable[Zscore],
    },
    parameters={
        "pair": List[Str],
        "threshold": Float,
        "permutation_num": Int,
        "min_size": Int,
        "max_size": Int,
        "spline_type": Str % Choices(SPLINE_TYPES),
        "degree": Int,
        "dof": Int,
        "p_val_thresh": Float,
        "nes_thresh": Float,
        "iter_tables_dir": Str,
        "max_workers": Int,
        "seed": Int,
    },
    outputs=[
        ("filtered_peptide_sets", GMT),
        ("last_psea_table", FeatureData[PSEAScores]),
    ],
    input_descriptions={
        "scores": "Processed Z-score table.",
        "peptide_sets": "Initial GMT for this pair.",
        "species_taxa": "Optional species-name to taxa-ID mapping.",
        "epitope_map": "Optional epitope mapping used for residual collapse.",
        "mapped_processed_scores": (
            "Optional processed epitope-level scores, required when epitope_map is provided."
        ),
    },
    parameter_descriptions={
        "pair": "Two sample IDs to compare.",
        "threshold": "Minimum Z score required to include peptides in GSEA.",
        "permutation_num": "Number of permutations in fgsea.",
        "min_size": "Minimum peptide set size.",
        "max_size": "Maximum peptide set size.",
        "spline_type": "Spline implementation used for residual fitting.",
        "degree": "Degree used by cubic spline mode.",
        "dof": "Degrees of freedom for cubic spline mode.",
        "p_val_thresh": "Adjusted p-value significance threshold.",
        "nes_thresh": "Absolute NES significance threshold.",
        "iter_tables_dir": "Optional directory for per-iteration tables.",
        "max_workers": "Reserved worker-count parameter for compatibility.",
        "seed": "Random seed used by fgsea.",
    },
    output_descriptions={
        "filtered_peptide_sets": "Pair-specific GMT after iterative filtering.",
        "last_psea_table": "Final PSEA table generated for the pair.",
    },
    name="Run Iterative Process Single Pair",
    description=(
        "Run iterative peptide-set filtering for a single pair by calling the "
        "registered pair FGSEA method."
    ),
)

plugin.pipelines.register_function(
    function=run_iterative_peptide_analysis,
    inputs={
        "scores": FeatureTable[Zscore],
        "pairs": PSEAPairs,
        "peptide_sets": GMT,
        "species_taxa": PSEASpeciesTaxa,
        "epitope_map": FeatureData[MappedEpitope],
        "mapped_processed_scores": FeatureTable[Zscore],
    },
    parameters={
        "threshold": Float,
        "permutation_num": Int,
        "min_size": Int,
        "max_size": Int,
        "spline_type": Str % Choices(SPLINE_TYPES),
        "degree": Int,
        "dof": Int,
        "p_val_thresh": Float,
        "nes_thresh": Float,
        "iter_tables_dir": Str,
        "max_workers": Int,
        "seed": Int,
    },
    outputs=[
        ("pair_peptide_sets", Collection[GMT]),
    ],
    input_descriptions={
        "scores": "Processed Z-score table.",
        "pairs": "Pairs table describing sample comparisons.",
        "peptide_sets": "Initial GMT for all pairs.",
        "species_taxa": "Optional species-name to taxa-ID mapping.",
        "epitope_map": "Optional epitope mapping used for residual collapse.",
        "mapped_processed_scores": (
            "Optional processed epitope-level scores, required when epitope_map is provided."
        ),
    },
    parameter_descriptions={
        "threshold": "Minimum Z score required to include peptides in GSEA.",
        "permutation_num": "Number of permutations in fgsea.",
        "min_size": "Minimum peptide set size.",
        "max_size": "Maximum peptide set size.",
        "spline_type": "Spline implementation used for residual fitting.",
        "degree": "Degree used by cubic spline mode.",
        "dof": "Degrees of freedom for cubic spline mode.",
        "p_val_thresh": "Adjusted p-value significance threshold.",
        "nes_thresh": "Absolute NES significance threshold.",
        "iter_tables_dir": "Optional directory for per-iteration tables.",
        "max_workers": "Reserved worker-count parameter for compatibility.",
        "seed": "Random seed used by fgsea.",
    },
    output_descriptions={
        "pair_peptide_sets": "Pair-specific GMTs after iterative filtering.",
    },
    name="Run Iterative Peptide Analysis",
    description=(
        "Run iterative peptide filtering across all pairs by calling the registered "
        "single-pair iterative pipeline."
    ),
)

plugin.pipelines.register_function(
    function=make_psea_table,
    inputs={
        "scores": FeatureTable[Zscore],
        "pairs": PSEAPairs,
        "peptide_sets": GMT,
        "species_taxa": PSEASpeciesTaxa,
        "species_color": PSEASpeciesColor,
        "epitope": FeatureData[Epitope],
    },
    parameters={
        "threshold": Float,
        "collapse": Str % Choices(["Bacterial", "Viral", "Both"]),
        "p_val_thresh": Float,
        "nes_thresh": Float,
        "min_size": Int,
        "max_size": Int,
        "permutation_num": Int,
        "spline_type": Str % Choices(SPLINE_TYPES),
        "degree": Int,
        "dof": Int,
        "table_dir": Str,
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
        ("psea_tables", Collection[FeatureData[PSEAScores]]),
    ],
    input_descriptions={
        "scores": "FeatureTable[Zscore] with peptide-level Z scores.",
        "pairs": "Tab-delimited table containing sample pairs.",
        "peptide_sets": "GMT mapping taxa IDs to peptide IDs.",
        "species_taxa": "Optional species-name to taxa-ID mapping.",
        "species_color": "Optional species color mapping used by plotting actions.",
        "epitope": "Optional epitope mapping input for epitope-collapsed analysis.",
    },
    parameter_descriptions={
        "threshold": "Minimum Z score required to include peptides in GSEA.",
        "collapse": "Epitope collapse strategy when epitope input is provided.",
        "p_val_thresh": (
            "Adjusted p-value threshold for highlighting taxa in output plots."
        ),
        "nes_thresh": "Absolute NES threshold for highlighting taxa in output plots.",
        "min_size": "Minimum peptide set size.",
        "max_size": "Maximum peptide set size.",
        "permutation_num": "Number of fgsea permutations.",
        "spline_type": "Spline implementation used for residual fitting.",
        "degree": "Degree used by cubic spline mode.",
        "dof": "Degrees of freedom for cubic spline mode.",
        "table_dir": "Directory where pair-level PSEA tables are written.",
        "iterative_analysis": "Whether to run iterative peptide filtering.",
        "iter_tables_dir": "Optional directory where iterative tables are written.",
        "max_workers": "Reserved worker-count parameter for compatibility.",
        "summary_tables_dir": "Directory where AE summary tables are written.",
        "vis_outputs_dir": "Optional directory for non-QIIME sidecar image outputs.",
        "seed": "Random seed used by fgsea.",
    },
    output_descriptions={
        "scatter_plot": "Scatter plot visualization across all pairs.",
        "volcano_plot": "Volcano plot visualization across all pairs.",
        "ae_plots": "Antibody event summary visualization.",
        "psea_tables": "Collection of pair-level PSEA tables.",
    },
    name="Make PSEA Table",
    description=(
        "Run PSEA across all sample pairs and produce pair-level tables and summary "
        "visualizations."
    ),
)

importlib.import_module("q2_pepsirf.transformers")
importlib.import_module("q2_PSEA.transformers")
