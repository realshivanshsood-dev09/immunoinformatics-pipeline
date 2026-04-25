from __future__ import annotations

from dataclasses import dataclass

VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class PreprocessResult:
    sequence: str
    removed_residues: int
    segments: list[tuple[int, str]]
    antigenicity_prefilter: float


def heuristic_antigenicity(sequence: str) -> float:
    """
    Estimate antigenicity using composition and residue diversity.
    Hydrophilic, charged, and aromatic exposure-prone residues increase score.
    """
    if not sequence:
        return 0.0

    exposure_favoring = set("DEHKRSTNQYP")
    hydrophobic_penalty = set("AILMFWV")

    exposure_fraction = sum(aa in exposure_favoring for aa in sequence) / len(sequence)
    hydrophobic_fraction = sum(aa in hydrophobic_penalty for aa in sequence) / len(sequence)
    diversity_fraction = len(set(sequence)) / 20.0

    score = 0.55 * exposure_fraction + 0.30 * diversity_fraction + 0.15 * (1 - hydrophobic_fraction)
    return round(max(0.0, min(1.0, score)), 4)


def clean_sequence(sequence: str) -> tuple[str, int]:
    cleaned = "".join(residue for residue in sequence.upper() if residue in VALID_AMINO_ACIDS)
    removed = len(sequence) - len(cleaned)
    return cleaned, removed


def segment_sequence(sequence: str, segment_length: int, overlap: int) -> list[tuple[int, str]]:
    if len(sequence) <= segment_length:
        return [(1, sequence)]

    segments: list[tuple[int, str]] = []
    step = max(1, segment_length - overlap)
    for start in range(0, len(sequence), step):
        end = min(start + segment_length, len(sequence))
        segment = sequence[start:end]
        if segment:
            segments.append((start + 1, segment))
        if end == len(sequence):
            break
    return segments


def preprocess_sequence(
    sequence: str,
    segment_length: int = 1200,
    overlap: int = 50,
) -> PreprocessResult:
    cleaned_sequence, removed = clean_sequence(sequence)
    prefilter_score = heuristic_antigenicity(cleaned_sequence)
    segments = segment_sequence(cleaned_sequence, segment_length=segment_length, overlap=overlap)
    return PreprocessResult(
        sequence=cleaned_sequence,
        removed_residues=removed,
        segments=segments,
        antigenicity_prefilter=prefilter_score,
    )
