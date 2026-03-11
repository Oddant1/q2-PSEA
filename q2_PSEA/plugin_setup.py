#! /usr/bin/env python
from q2_types.feature_data import FeatureData
from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import PSEAScores
from q2_pepsirf.types import Zscore

from q2_PSEA.formats import (
    PairsTSVFormat,
    PeptideSetsFormat,
    SpeciesTaxaTSVFormat,
    SpeciesColorsTSVFormat,
)
from q2_PSEA.types import Pairs, PeptideSets, SpeciesTaxa, SpeciesColors

import q2_PSEA

from q2_PSEA.actions.psea import make_psea_table
from qiime2.plugin import (
    Bool, Float, Int, Plugin, Str, Visualization, Choices, Collection
)


# q2-PSEA plugin object
plugin = Plugin(
    "psea", version=q2_PSEA.__version__,
    website="https://github.com/LadnerLab/q2-PSEA.git",
    description="Qiime2 Plugin for PSEA."  # TODO: get a description
)

plugin.register_semantic_types(Pairs, PeptideSets, SpeciesTaxa, SpeciesColors)
plugin.register_formats(
    PairsTSVFormat,
    PeptideSetsFormat,
    SpeciesTaxaTSVFormat,
    SpeciesColorsTSVFormat,
)
plugin.register_semantic_type_to_format(Pairs, PairsTSVFormat)
plugin.register_semantic_type_to_format(PeptideSets, PeptideSetsFormat)
plugin.register_semantic_type_to_format(SpeciesTaxa, SpeciesTaxaTSVFormat)
plugin.register_semantic_type_to_format(SpeciesColors, SpeciesColorsTSVFormat)


# register make_psea_table function
plugin.pipelines.register_function(
    function=make_psea_table,
    inputs={
        "scores": FeatureTable[Zscore],
        "pairs": Pairs,
        "peptide_sets": PeptideSets,
    },
    parameters={
        "threshold": Float,
        "epitope": Str,
        "collapse": Str % Choices(['Bacterial', 'Viral', 'Both']),
        "p_val_thresh": Float,
        "nes_thresh": Float,
        "species_taxa": Str,
        "species_colors": Str,
        "min_size": Int,
        "max_size": Int,
        "permutation_num": Int,
        "spline_type": Str % Choices(q2_PSEA.actions.splines.SPLINE_TYPES),
        "degree": Int,
        "dof": Int,
        "iterative_analysis": Bool,
        "iter_tables_dir": Str,
        "max_workers": Int,
        "seed": Int
    },
    parameter_descriptions={
        "threshold": "Minimum Z score a peptide must maintain to be"
            " considered in Gene Set Enrichment Analysis.",
        "epitope": "File containing information relating peptides to"
            " epitopes and species/subtype. If this argument is passed in the"
            " peptide level residuals will be calculated then collapsed to"
            " epitope level prior to further analysis.",
        "collapse": "Which Category we should collapse to epitope level. Note that"
                    " this parameter means nothing if the epitope input isn't used.",
        "p_val_thresh": "Specifies the value adjusted p-values must meet to be"
            " considered for highlighting in volcano and scatter plots.",
        "nes_thresh": "Specifies the value ",
        "species_taxa": "Tab-delimited file containing species name and taxonomy ID associations.",
        "species_colors": "Tab-delimited file containing species name and HEX color code for charts.",
        "min_size": "Minimum allowed number of peptides from peptide set also"
            " the data set.",
        "max_size": "Maximum allowed number of peptides from peptide set also"
            " the data set.",
        "permutation_num": "Number of permutations. Minimal possible nominal"
            " p-value is about 1/perm.",
        "spline_type": "Specifies which spline operation to use.",
        "degree": "Specifies the degree of the piecewise polynomial. Note: at"
            " the moment, this will only affect the `cubic` spline approach.",
        "dof": "Degree of freedom to use when fitting the spline. Note: at the"
            " moment, this will only affect the `cubic` spline approach.",
        "iterative_analysis": "Boolean value, whether or not to use iterative approach"
            " to filter cross-reactive peptides from less significant species."
            " GMT peptide_sets recommended.",
        "iter_tables_dir": "Directory name to output iteration tables to. Only generates if name is provided.",
        "max_workers": "Maximum number of processes to run at a time. If none set,"
            " defaults to the number of processors on the machine.",
        "seed": "Seed for permutation. Seed used to generate a random number for phenotype and gene_set permutations when running GSEA.",
    },
    outputs=[
        ("scatter_plot", Visualization), ("volcano_plot", Visualization),
        ("ae_plots", Visualization),
        ("psea_tables", Collection[FeatureData[PSEAScores]])
    ],
    output_descriptions={
        "scatter_plot": "Name of plot file visualization comparison between"
            " two samples. This plot includes the smooth spline fit to the"
            " given data and highlights the leading edge peptides for all"
            " significant taxa.",
        "volcano_plot": "Name of plot file visualization comparison between"
            " enrichment scores (ES) and p-values.",
        "ae_plots": "Name of plot of visualization for event summary, if created."
    },
    name="Make PSEA Table",
    description="Qiime2 plug-in which provides a Python wrapper around R"
        " functions to perform Peptide Set Enrichment Analysis. A Z score"
        " scatter plot for each sample replicate pair complete with a spline"
        " and highlighting of significant taxa. As well as a volcano plot"
        " showing relationship between adjusted p-values and normalized"
        " enrichment scores (NES) with highlighting of points which pass"
        " provided thresholds."
)
