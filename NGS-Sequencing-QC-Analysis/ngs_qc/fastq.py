"""Small, strict four-line FASTQ reader for a reproducible teaching example."""
from dataclasses import dataclass
import gzip
from pathlib import Path

DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class Read:
    name: str
    sequence: str
    encoded_quality: str

    @property
    def scores(self) -> list[int]:
        return [ord(char) - 33 for char in self.encoded_quality]

    def to_fastq(self) -> str:
        return f"@{self.name}\n{self.sequence}\n+\n{self.encoded_quality}\n"


def read_fastq(path: str | Path):
    """Yield FASTQ reads; reject truncated/invalid four-line records.

    Uses Phred+33. No auto-detection of Phred+64 or wrapped FASTQ records.
    Supports .fastq and .fastq.gz input.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"FASTQ input not found: {source}")
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    with opener(source, mode="rt", encoding="ascii", newline=None) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            record_number += 1
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            where = f"record {record_number} in {source.name}"
            if not sequence or not plus or not quality:
                raise ValueError(f"Truncated FASTQ {where}")
            if not header.startswith("@") or not header[1:].strip():
                raise ValueError(f"Invalid FASTQ header in {where}")
            if not plus.startswith("+"):
                raise ValueError(f"Invalid FASTQ separator in {where}")
            sequence = sequence.rstrip("\r\n").upper()
            quality = quality.rstrip("\r\n")
            if not sequence or len(sequence) != len(quality):
                raise ValueError(f"Sequence / quality length mismatch in {where}")
            if set(sequence) - DNA_IUPAC:
                raise ValueError(f"Non-IUPAC nucleotide in {where}")
            if any(not (33 <= ord(char) <= 126) for char in quality):
                raise ValueError(f"Invalid Phred+33 quality character in {where}")
            yield Read(header[1:].strip(), sequence, quality)
