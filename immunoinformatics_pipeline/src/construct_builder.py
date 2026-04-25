from __future__ import annotations

from src.epitope_prediction import EpitopeRecord


def select_top_epitopes(
    ranked_epitopes: list[EpitopeRecord],
    top_n: int,
    max_bcell_epitopes: int = 3,
) -> list[EpitopeRecord]:
    selected: list[EpitopeRecord] = []
    seen_peptides: set[str] = set()
    bcell_count = 0

    for record in ranked_epitopes:
        peptide = str(record["peptide"])
        if peptide in seen_peptides:
            continue
        if record["epitope_class"] == "B_CELL":
            if bcell_count >= max_bcell_epitopes:
                continue
            bcell_count += 1
        selected.append(dict(record))
        seen_peptides.add(peptide)
        if len(selected) >= top_n:
            break
    return selected


def build_vaccine_construct(selected_epitopes: list[EpitopeRecord], construct_config: dict) -> str:
    if not selected_epitopes:
        return ""

    ctl_linker = str(construct_config.get("ctl_linker", "AAY"))
    htl_linker = str(construct_config.get("htl_linker", "GPGPG"))
    adjuvant_sequence = str(construct_config.get("adjuvant_sequence", ""))
    add_adjuvant = bool(construct_config.get("add_adjuvant", True))
    adjuvant_linker = str(construct_config.get("adjuvant_linker", "EAAAK"))

    ctl_epitopes = [str(item["peptide"]) for item in selected_epitopes if item["epitope_class"] == "MHC_I"]
    htl_epitopes = [str(item["peptide"]) for item in selected_epitopes if item["epitope_class"] == "MHC_II"]
    bcell_epitopes = [str(item["peptide"]) for item in selected_epitopes if item["epitope_class"] == "B_CELL"]

    construct_sections: list[str] = []
    if add_adjuvant and adjuvant_sequence:
        construct_sections.append(adjuvant_sequence)
        construct_sections.append(adjuvant_linker)
    if ctl_epitopes:
        construct_sections.append(ctl_linker.join(ctl_epitopes))
    if htl_epitopes:
        construct_sections.append(htl_linker.join(htl_epitopes))
    if bcell_epitopes:
        construct_sections.append("KK".join(bcell_epitopes))

    return "GPGPG".join(section for section in construct_sections if section)
