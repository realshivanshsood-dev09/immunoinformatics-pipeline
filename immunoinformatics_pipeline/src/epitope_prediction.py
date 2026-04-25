from __future__ import annotations

import csv
from pathlib import Path

from src.preprocess import heuristic_antigenicity


EpitopeRecord = dict[str, object]


def load_netmhc_predictions(file_path: str | Path, epitope_class: str) -> list[EpitopeRecord]:
    """Load NetMHCpan/NetMHCIIpan outputs into a unified record structure."""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Prediction file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required_columns = {"peptide", "position", "allele", "ic50", "rank"}
        missing = required_columns - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {path.name}: {sorted(missing)}")

        source_tool = "NetMHCpan" if epitope_class == "MHC_I" else "NetMHCIIpan"
        predictions: list[EpitopeRecord] = []
        for row in reader:
            peptide = row["peptide"].strip().upper()
            predictions.append(
                {
                    "peptide": peptide,
                    "position": int(float(row["position"])),
                    "allele": row["allele"].strip(),
                    "binding_affinity": float(row["ic50"]),
                    "binding_rank": float(row["rank"]),
                    "epitope_class": epitope_class,
                    "source_tool": source_tool,
                    "peptide_length": len(peptide),
                    "bepipred_score": 0.0,
                }
            )
    return predictions


def _bcell_window_score(window: str) -> float:
    hydrophilic = set("DEHKNQRSTYP")
    flexible = set("GSPDNT")
    cysteine_penalty = window.count("C") / max(1, len(window))
    hydrophilic_fraction = sum(aa in hydrophilic for aa in window) / len(window)
    flexible_fraction = sum(aa in flexible for aa in window) / len(window)
    antigenicity = heuristic_antigenicity(window)
    score = 0.45 * hydrophilic_fraction + 0.25 * flexible_fraction + 0.30 * antigenicity - 0.10 * cysteine_penalty
    return round(max(0.0, min(1.0, score)), 4)


def predict_bcell_epitopes(
    sequence: str,
    window_size: int = 14,
    min_score: float = 0.58,
) -> list[EpitopeRecord]:
    """
    Deterministic BepiPred-style simulation using hydrophilicity, flexibility,
    and exposure-associated composition.
    """
    if len(sequence) < window_size:
        return []

    predictions: list[EpitopeRecord] = []
    for start in range(0, len(sequence) - window_size + 1):
        peptide = sequence[start : start + window_size]
        score = _bcell_window_score(peptide)
        if score >= min_score:
            predictions.append(
                {
                    "peptide": peptide,
                    "position": start + 1,
                    "allele": "B-cell-linear",
                    "binding_affinity": round(1000.0 - (score * 700.0), 3),
                    "binding_rank": round(max(0.1, (1 - score) * 12), 3),
                    "epitope_class": "B_CELL",
                    "source_tool": "BepiPred_simulated",
                    "peptide_length": len(peptide),
                    "bepipred_score": score,
                }
            )

    predictions.sort(key=lambda record: (-float(record["bepipred_score"]), int(record["position"])))
    return predictions


def merge_epitope_predictions(*record_sets: list[EpitopeRecord]) -> list[EpitopeRecord]:
    merged: list[EpitopeRecord] = []
    for record_set in record_sets:
        merged.extend(record_set)
    merged.sort(key=lambda record: (int(record["position"]), str(record["epitope_class"])))
    return merged
