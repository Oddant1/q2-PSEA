import pandas as pd

from q2_PSEA.format_types import (
    PSEAPairsFormat,
    PSEASpeciesColorFormat,
    PSEASpeciesTaxaFormat,
)
from q2_PSEA.plugin_setup import plugin


@plugin.register_transformer
def _pairs_to_df(ff: PSEAPairsFormat) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep="\t")


@plugin.register_transformer
def _df_to_pairs(df: pd.DataFrame) -> PSEAPairsFormat:
    result = PSEAPairsFormat()
    df.to_csv(str(result), sep="\t", index=False)
    return result


@plugin.register_transformer
def _species_taxa_to_df(ff: PSEASpeciesTaxaFormat) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep="\t", header=None)


@plugin.register_transformer
def _df_to_species_taxa(df: pd.DataFrame) -> PSEASpeciesTaxaFormat:
    result = PSEASpeciesTaxaFormat()
    df.to_csv(str(result), sep="\t", index=False, header=False)
    return result


@plugin.register_transformer
def _species_color_to_df(ff: PSEASpeciesColorFormat) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep="\t", header=None)


@plugin.register_transformer
def _df_to_species_color(df: pd.DataFrame) -> PSEASpeciesColorFormat:
    result = PSEASpeciesColorFormat()
    df.to_csv(str(result), sep="\t", index=False, header=False)
    return result
