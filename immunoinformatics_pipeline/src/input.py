from __future__ import annotations

from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(fasta_path: str | Path) -> str:
    """Parse a FASTA file and return the concatenated protein sequence."""
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    lines = path.read_text(encoding="utf-8").splitlines()
    sequence_parts: list[str] = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        sequence_parts.append(line)

    if not sequence_parts:
        raise ValueError("No protein sequence found in FASTA file.")

    sequence = "".join(sequence_parts).upper()
    validate_sequence(sequence)
    return sequence


def validate_sequence(sequence: str) -> None:
    """Validate that the sequence contains only standard amino acid residues."""
    invalid_residues = sorted(set(sequence) - VALID_AMINO_ACIDS)
    if invalid_residues:
        raise ValueError(
            "Invalid amino acid residues detected: " + ", ".join(invalid_residues)
        )


def fetch_uniprot_fasta(accession: str) -> str:
    """
    Fetch a UniProt FASTA record using the public REST endpoint and return it as text.
    """
    accession = accession.strip()
    if not accession:
        raise ValueError("UniProt accession cannot be empty.")

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        with urlopen(url, timeout=30) as response:
            payload = response.read().decode("utf-8")
    except HTTPError as exc:
        raise ValueError(f"UniProt accession lookup failed for {accession}: HTTP {exc.code}") from exc
    except URLError as exc:
        raise ConnectionError(f"Unable to reach UniProt for accession {accession}.") from exc

    if not payload.startswith(">"):
        raise ValueError(f"UniProt did not return a FASTA record for {accession}.")
    return payload


def save_fasta_text(fasta_text: str, output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(fasta_text, encoding="utf-8")
    return path
