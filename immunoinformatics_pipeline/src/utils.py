from __future__ import annotations

import csv
from pathlib import Path

from src.epitope_prediction import EpitopeRecord


def ensure_output_dir(output_dir: str | Path) -> Path:
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_records(records: list[EpitopeRecord], output_path: str | Path) -> None:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        path.write_text("", encoding="utf-8")
        return

    fieldnames: list[str] = []
    for record in records:
        for key in record.keys():
            if key not in fieldnames:
                fieldnames.append(key)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def _svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white" />',
    ]


def plot_score_distribution(ranked_epitopes: list[EpitopeRecord], output_path: str | Path) -> None:
    path = Path(output_path)
    if not ranked_epitopes:
        path.write_text("", encoding="utf-8")
        return

    scores = [float(record["composite_score"]) for record in ranked_epitopes]
    width, height = 820, 420
    margin = 50
    bins = 10
    minimum = min(scores)
    maximum = max(scores)
    if maximum == minimum:
        maximum = minimum + 1e-6

    counts = [0] * bins
    for score in scores:
        index = min(int(((score - minimum) / (maximum - minimum)) * bins), bins - 1)
        counts[index] += 1

    chart_width = width - 2 * margin
    chart_height = height - 2 * margin
    bar_width = chart_width / bins
    max_count = max(counts) or 1

    svg = _svg_header(width, height)
    svg.extend(
        [
            f'<line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" stroke="black" />',
            f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height - margin}" stroke="black" />',
            f'<text x="{width / 2}" y="28" text-anchor="middle" font-size="18">Ranked Epitope Score Distribution</text>',
        ]
    )
    for idx, count in enumerate(counts):
        bar_height = (count / max_count) * (chart_height - 10)
        x = margin + idx * bar_width + 4
        y = height - margin - bar_height
        svg.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_width - 8:.2f}" height="{bar_height:.2f}" fill="#2f6f8f" />'
        )
    svg.append("</svg>")
    path.write_text("\n".join(svg), encoding="utf-8")


def plot_epitope_positions(raw_epitopes: list[EpitopeRecord], sequence_length: int, output_path: str | Path) -> None:
    path = Path(output_path)
    if not raw_epitopes:
        path.write_text("", encoding="utf-8")
        return

    width, height = 1100, 300
    margin = 70
    class_y = {"MHC_I": 90, "MHC_II": 150, "B_CELL": 210}
    colors = {"MHC_I": "#d95f02", "MHC_II": "#1b9e77", "B_CELL": "#7570b3"}

    svg = _svg_header(width, height)
    svg.extend(
        [
            f'<text x="{width / 2}" y="28" text-anchor="middle" font-size="18">Epitope Position Mapping Across the Protein</text>',
            f'<line x1="{margin}" y1="250" x2="{width - margin}" y2="250" stroke="black" />',
        ]
    )
    for label, y in class_y.items():
        svg.append(f'<text x="20" y="{y + 5}" font-size="14">{label}</text>')
        svg.append(f'<line x1="{margin}" y1="{y}" x2="{width - margin}" y2="{y}" stroke="#dddddd" stroke-dasharray="4,4" />')

    for record in raw_epitopes:
        position = int(record["position"])
        epitope_class = str(record["epitope_class"])
        x = margin + ((position - 1) / max(1, sequence_length - 1)) * (width - 2 * margin)
        y = class_y.get(epitope_class, 120)
        color = colors.get(epitope_class, "#666666")
        svg.append(f'<circle cx="{x:.2f}" cy="{y}" r="5" fill="{color}" />')

    svg.append("</svg>")
    path.write_text("\n".join(svg), encoding="utf-8")


def write_markdown_report(
    output_path: str | Path,
    preprocess_summary: dict,
    raw_epitopes: list[EpitopeRecord],
    filtered_epitopes: list[EpitopeRecord],
    ranked_epitopes: list[EpitopeRecord],
    vaccine_construct: str,
) -> None:
    top_records = ranked_epitopes[:10]
    lines = [
        "# Immunoinformatics Epitope-to-Vaccine Report",
        "",
        "## Sequence Summary",
        f"- Clean sequence length: {preprocess_summary['sequence_length']}",
        f"- Removed residues: {preprocess_summary['removed_residues']}",
        f"- Antigenicity pre-filter score: {preprocess_summary['antigenicity_prefilter']:.4f}",
        f"- Segments generated: {preprocess_summary['segment_count']}",
        "",
        "## Pipeline Yield",
        f"- Raw epitopes: {len(raw_epitopes)}",
        f"- Filtered epitopes: {len(filtered_epitopes)}",
        f"- Ranked epitopes: {len(ranked_epitopes)}",
        "",
        "## Top Ranked Epitopes",
        "| Peptide | Class | Position | Rank | Antigenicity | Toxicity | Allergenicity | Composite |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for record in top_records:
        lines.append(
            "| {peptide} | {epitope_class} | {position} | {binding_rank} | {antigenicity} | {toxicity} | {allergenicity} | {composite_score} |".format(
                **record
            )
        )
    if not top_records:
        lines.append("| No ranked epitopes available | - | - | - | - | - | - | - |")

    lines.extend(
        [
            "",
            "## Vaccine Construct",
            f"- Construct length: {len(vaccine_construct)} aa",
            f"- Sequence: `{vaccine_construct}`",
        ]
    )
    Path(output_path).write_text("\n".join(lines), encoding="utf-8")
