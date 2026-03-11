from qiime2.plugin import TextFileFormat, ValidationError


class PairsTSVFormat(TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            header = fh.readline()
            if not header:
                raise ValidationError("Pairs file is empty.")


class PeptideSetsFormat(TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            first = fh.readline()
            if not first:
                raise ValidationError("Peptide sets file is empty.")


class SpeciesTaxaTSVFormat(TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            first = fh.readline()
            if not first:
                raise ValidationError("Species taxa file is empty.")


class SpeciesColorsTSVFormat(TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            first = fh.readline()
            if not first:
                raise ValidationError("Species colors file is empty.")
