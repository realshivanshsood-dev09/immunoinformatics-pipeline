from __future__ import annotations

from src.epitope_prediction import EpitopeRecord
from src.preprocess import heuristic_antigenicity


def simulate_vaxijen(peptide: str) -> float:
    return heuristic_antigenicity(peptide)


def simulate_allertop(peptide: str) -> float:
    if not peptide:
        return 1.0
    aromatic = sum(aa in "FWY" for aa in peptide) / len(peptide)
    cysteine = peptide.count("C") / len(peptide)
    motif_repetition = 1.0 - (len(set(peptide)) / len(peptide))
    charge_balance = abs(sum(aa in "KRH" for aa in peptide) - sum(aa in "DE" for aa in peptide)) / len(peptide)
    score = 0.35 * aromatic + 0.30 * cysteine + 0.20 * motif_repetition + 0.15 * charge_balance
    return round(max(0.0, min(1.0, score)), 4)


def simulate_toxinpred(peptide: str) -> float:
    if not peptide:
        return 1.0
    hydrophobic = sum(aa in "AILMFWV" for aa in peptide) / len(peptide)
    basic = sum(aa in "KRH" for aa in peptide) / len(peptide)
    proline_penalty = peptide.count("P") / len(peptide)
    glycine_serine_relief = sum(aa in "GS" for aa in peptide) / len(peptide)
    score = 0.45 * hydrophobic + 0.30 * basic + 0.15 * proline_penalty - 0.20 * glycine_serine_relief
    return round(max(0.0, min(1.0, score)), 4)


def annotate_biological_filters(epitopes: list[EpitopeRecord]) -> list[EpitopeRecord]:
    annotated: list[EpitopeRecord] = []
    for record in epitopes:
        item = dict(record)
        peptide = str(item["peptide"])
        item["antigenicity"] = simulate_vaxijen(peptide)
        item["allergenicity"] = simulate_allertop(peptide)
        item["toxicity"] = simulate_toxinpred(peptide)
        annotated.append(item)
    return annotated


def filter_epitopes(epitopes: list[EpitopeRecord], thresholds: dict) -> list[EpitopeRecord]:
    annotated = annotate_biological_filters(epitopes)
    class_rank_limits = thresholds.get("max_rank", {})
    max_mhc_i_rank = float(class_rank_limits.get("mhc_i", 2.0))
    max_mhc_ii_rank = float(class_rank_limits.get("mhc_ii", 10.0))
    min_antigenicity = float(thresholds.get("min_antigenicity", 0.45))
    max_toxicity = float(thresholds.get("max_toxicity", 0.45))
    max_allergenicity = float(thresholds.get("max_allergenicity", 0.55))

    filtered: list[EpitopeRecord] = []
    for item in annotated:
        epitope_class = str(item["epitope_class"])
        binding_rank = float(item["binding_rank"])
        passes_binding = True
        if epitope_class == "MHC_I":
            passes_binding = binding_rank <= max_mhc_i_rank
        elif epitope_class == "MHC_II":
            passes_binding = binding_rank <= max_mhc_ii_rank

        item["passes_antigenicity"] = float(item["antigenicity"]) >= min_antigenicity
        item["passes_toxicity"] = float(item["toxicity"]) <= max_toxicity
        item["passes_allergenicity"] = float(item["allergenicity"]) <= max_allergenicity
        item["passes_binding"] = passes_binding
        item["passes_filters"] = (
            bool(item["passes_antigenicity"])
            and bool(item["passes_toxicity"])
            and bool(item["passes_allergenicity"])
            and bool(item["passes_binding"])
        )

        if item["passes_filters"]:
            filtered.append(item)

    filtered.sort(
        key=lambda record: (
            str(record["epitope_class"]),
            float(record["binding_rank"]),
            -float(record["antigenicity"]),
        )
    )
    return filtered
