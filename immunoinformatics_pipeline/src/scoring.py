from __future__ import annotations

from src.epitope_prediction import EpitopeRecord


def _min_max_normalize(values: list[float], invert: bool = False) -> list[float]:
    if not values:
        return []
    minimum = min(values)
    maximum = max(values)
    if maximum == minimum:
        normalized = [1.0 for _ in values]
    else:
        normalized = [(value - minimum) / (maximum - minimum) for value in values]
    if invert:
        normalized = [1.0 - value for value in normalized]
    return [max(0.0, min(1.0, value)) for value in normalized]


def score_epitopes(epitopes: list[EpitopeRecord], weights: dict) -> list[EpitopeRecord]:
    if not epitopes:
        return []

    scored = [dict(record) for record in epitopes]
    antigenicity_values = [float(record["antigenicity"]) for record in scored]
    affinity_values = [float(record["binding_affinity"]) for record in scored]
    toxicity_values = [float(record["toxicity"]) for record in scored]
    allergenicity_values = [float(record["allergenicity"]) for record in scored]

    norm_antigenicity = _min_max_normalize(antigenicity_values)
    norm_affinity = _min_max_normalize(affinity_values, invert=True)
    norm_toxicity = _min_max_normalize(toxicity_values)
    norm_allergenicity = _min_max_normalize(allergenicity_values)

    w_antigenicity = float(weights.get("antigenicity", 0.35))
    w_binding = float(weights.get("binding_affinity", 0.40))
    w_toxicity = float(weights.get("toxicity", 0.15))
    w_allergenicity = float(weights.get("allergenicity", 0.10))

    for index, record in enumerate(scored):
        record["norm_antigenicity"] = round(norm_antigenicity[index], 4)
        record["norm_binding_affinity"] = round(norm_affinity[index], 4)
        record["norm_toxicity"] = round(norm_toxicity[index], 4)
        record["norm_allergenicity"] = round(norm_allergenicity[index], 4)
        score = (
            w_antigenicity * norm_antigenicity[index]
            + w_binding * norm_affinity[index]
            - w_toxicity * norm_toxicity[index]
            - w_allergenicity * norm_allergenicity[index]
        )
        record["composite_score"] = round(score, 4)

    scored.sort(
        key=lambda record: (
            -float(record["composite_score"]),
            float(record["binding_rank"]),
            -float(record["antigenicity"]),
        )
    )
    return scored
