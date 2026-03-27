import json

from qiime2.plugin import ValidationError, model


class PSEAPairsFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            header = fh.readline()
            if not header:
                raise ValidationError("Pairs file is empty.")

            observed_data_row = False
            for line_number, line in enumerate(fh, start=2):
                stripped = line.strip()
                if not stripped:
                    continue

                if len(stripped.split("\t")) < 2:
                    raise ValidationError(
                        f"Pairs row {line_number} has fewer than two tab-delimited columns."
                    )

                observed_data_row = True
                if level == "min" and line_number >= 5:
                    break

            if not observed_data_row:
                raise ValidationError("Pairs file must contain at least one pair row.")


class PeptideSetsFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            for line in fh:
                if line.strip():
                    return
        raise ValidationError("Peptide sets file is empty.")


class SpeciesTaxonomyFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            observed_data_row = False
            for line_number, line in enumerate(fh, start=1):
                stripped = line.strip()
                if not stripped:
                    continue

                if len(stripped.split("\t")) < 2:
                    raise ValidationError(
                        f"Species taxonomy row {line_number} has fewer than two tab-delimited columns."
                    )

                observed_data_row = True
                if level == "min" and line_number >= 5:
                    break

            if not observed_data_row:
                raise ValidationError("Species taxonomy file is empty.")


class SpeciesColorsFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            observed_data_row = False
            for line_number, line in enumerate(fh, start=1):
                stripped = line.strip()
                if not stripped:
                    continue

                if len(stripped.split("\t")) < 2:
                    raise ValidationError(
                        f"Species colors row {line_number} has fewer than two tab-delimited columns."
                    )

                observed_data_row = True
                if level == "min" and line_number >= 5:
                    break

            if not observed_data_row:
                raise ValidationError("Species colors file is empty.")


class PSEATableFormat(model.TextFileFormat):
    def _validate_(self, level):
        required_columns = {"ID", "NES", "p.adjust", "all_tested_peptides"}
        with self.open() as fh:
            header = fh.readline().strip()
            if not header:
                raise ValidationError("PSEA table is empty.")

            observed = set(header.split("\t"))
            missing = required_columns - observed
            if missing:
                raise ValidationError(
                    f"PSEA table header is missing required columns: {sorted(missing)}"
                )


class SplineDataFormat(model.TextFileFormat):
    def _validate_(self, level):
        required_columns = {"x", "y", "pair"}
        with self.open() as fh:
            header = fh.readline().strip()
            if not header:
                raise ValidationError("Spline data table is empty.")

            observed = set(header.split("\t"))
            missing = required_columns - observed
            if missing:
                raise ValidationError(
                    f"Spline data header is missing required columns: {sorted(missing)}"
                )


class IterativePairStateFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            try:
                state = json.load(fh)
            except json.JSONDecodeError as e:
                raise ValidationError(f"Iterative pair state is not valid JSON: {e}")

        if not isinstance(state, dict):
            raise ValidationError("Iterative pair state must be a JSON object.")
        if "tested_species" not in state or "gmt_dict" not in state:
            raise ValidationError(
                "Iterative pair state JSON must contain 'tested_species' and 'gmt_dict'."
            )


class IterativePeptideSetsFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            try:
                state = json.load(fh)
            except json.JSONDecodeError as e:
                raise ValidationError(f"Iterative peptide sets are not valid JSON: {e}")

        if not isinstance(state, dict):
            raise ValidationError("Iterative peptide sets must be a JSON object.")


PSEAPairsDirFmt = model.SingleFileDirectoryFormat(
    "PSEAPairsDirFmt",
    "pairs.tsv",
    PSEAPairsFormat
)

PeptideSetsDirFmt = model.SingleFileDirectoryFormat(
    "PeptideSetsDirFmt",
    "peptide_sets.txt",
    PeptideSetsFormat
)

SpeciesTaxonomyDirFmt = model.SingleFileDirectoryFormat(
    "SpeciesTaxonomyDirFmt",
    "species_taxonomy.tsv",
    SpeciesTaxonomyFormat
)

SpeciesColorsDirFmt = model.SingleFileDirectoryFormat(
    "SpeciesColorsDirFmt",
    "species_colors.tsv",
    SpeciesColorsFormat
)

PSEATableDirFmt = model.SingleFileDirectoryFormat(
    "PSEATableDirFmt",
    "psea_table.tsv",
    PSEATableFormat
)

SplineDataDirFmt = model.SingleFileDirectoryFormat(
    "SplineDataDirFmt",
    "spline_data.tsv",
    SplineDataFormat
)

IterativePairStateDirFmt = model.SingleFileDirectoryFormat(
    "IterativePairStateDirFmt",
    "iterative_pair_state.json",
    IterativePairStateFormat
)

IterativePeptideSetsDirFmt = model.SingleFileDirectoryFormat(
    "IterativePeptideSetsDirFmt",
    "iterative_peptide_sets.json",
    IterativePeptideSetsFormat
)
