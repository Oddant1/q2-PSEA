from qiime2.plugin import SemanticType

import pandas as pd
import qiime2.plugin.model as model


PSEAPairs = SemanticType("PSEAPairs")
PSEASpeciesTaxa = SemanticType("PSEASpeciesTaxa")
PSEASpeciesColor = SemanticType("PSEASpeciesColor")


def _validate_tsv(path, min_columns, format_name, header):
    try:
        df = pd.read_csv(path, sep="\t", header=header)
    except Exception as exc:
        raise model.ValidationError(
            f"Could not parse {format_name} as a tab-delimited table: {exc}"
        ) from exc

    if df.shape[0] == 0:
        raise model.ValidationError(f"{format_name} is empty.")

    if df.shape[1] < min_columns:
        raise model.ValidationError(
            f"{format_name} must have at least {min_columns} columns."
        )


class PSEAPairsFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        _validate_tsv(self.path, min_columns=2, format_name="PSEAPairsFormat", header=0)


PSEAPairsDirFmt = model.SingleFileDirectoryFormat(
    "PSEAPairsDirFmt",
    "pairs.tsv",
    PSEAPairsFormat,
)


class PSEASpeciesTaxaFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        _validate_tsv(
            self.path,
            min_columns=2,
            format_name="PSEASpeciesTaxaFormat",
            header=None,
        )


PSEASpeciesTaxaDirFmt = model.SingleFileDirectoryFormat(
    "PSEASpeciesTaxaDirFmt",
    "species-taxa.tsv",
    PSEASpeciesTaxaFormat,
)


class PSEASpeciesColorFormat(model.TextFileFormat):
    def _validate_(self, level="min"):
        _validate_tsv(
            self.path,
            min_columns=2,
            format_name="PSEASpeciesColorFormat",
            header=None,
        )


PSEASpeciesColorDirFmt = model.SingleFileDirectoryFormat(
    "PSEASpeciesColorDirFmt",
    "species-colors.tsv",
    PSEASpeciesColorFormat,
)
