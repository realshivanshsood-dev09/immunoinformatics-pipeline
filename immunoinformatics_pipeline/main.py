from __future__ import annotations

import ast
import argparse
from pathlib import Path

from src.construct_builder import build_vaccine_construct, select_top_epitopes
from src.epitope_prediction import load_netmhc_predictions, merge_epitope_predictions, predict_bcell_epitopes
from src.filtering import filter_epitopes
from src.input import fetch_uniprot_fasta, parse_fasta, save_fasta_text
from src.preprocess import preprocess_sequence
from src.scoring import score_epitopes
from src.utils import ensure_output_dir, plot_epitope_positions, plot_score_distribution, save_records, write_markdown_report


def _parse_scalar(value: str):
    lowered = value.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    try:
        return ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return value.strip("'\"")


def load_config(config_path: str | Path) -> dict:
    """
    Lightweight YAML loader for the bundled config schema.
    Supports nested mappings and inline lists, which keeps the project runnable
    without external dependencies.
    """
    path = Path(config_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    root: dict = {}
    stack: list[tuple[int, dict]] = [(-1, root)]

    for raw_line in lines:
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        key, _, remainder = raw_line.strip().partition(":")

        while len(stack) > 1 and indent <= stack[-1][0]:
            stack.pop()
        current = stack[-1][1]

        if remainder.strip() == "":
            child: dict = {}
            current[key] = child
            stack.append((indent, child))
        else:
            current[key] = _parse_scalar(remainder.strip())

    return root


def _resolve_config_path(base_dir: Path, value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        return str(path)
    return str((base_dir / path).resolve())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Immunoinformatics epitope-to-vaccine pipeline"
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to the pipeline configuration file.",
    )
    parser.add_argument(
        "--fasta",
        help="Path to a local protein FASTA file. Overrides config.yaml input_fasta.",
    )
    parser.add_argument(
        "--uniprot-accession",
        help="UniProt accession to fetch as FASTA before running the pipeline.",
    )
    parser.add_argument(
        "--netmhc-i",
        help="Path to NetMHCpan class I CSV file. Overrides config.yaml netmhc_i_file.",
    )
    parser.add_argument(
        "--netmhc-ii",
        help="Path to NetMHCIIpan class II CSV file. Overrides config.yaml netmhc_ii_file.",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for results. Overrides config.yaml output_dir.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.fasta and args.uniprot_accession:
        raise ValueError("Use either --fasta or --uniprot-accession, not both.")

    config_path = args.config
    config_path = Path(config_path).resolve()
    config = load_config(config_path)
    base_dir = config_path.parent

    config["input_fasta"] = _resolve_config_path(base_dir, str(config["input_fasta"]))
    config["netmhc_i_file"] = _resolve_config_path(base_dir, str(config["netmhc_i_file"]))
    config["netmhc_ii_file"] = _resolve_config_path(base_dir, str(config["netmhc_ii_file"]))
    config["output_dir"] = _resolve_config_path(base_dir, str(config.get("output_dir", "results")))

    if args.netmhc_i:
        config["netmhc_i_file"] = _resolve_config_path(base_dir, args.netmhc_i)
    if args.netmhc_ii:
        config["netmhc_ii_file"] = _resolve_config_path(base_dir, args.netmhc_ii)
    if args.output_dir:
        config["output_dir"] = _resolve_config_path(base_dir, args.output_dir)

    if args.fasta:
        config["input_fasta"] = _resolve_config_path(base_dir, args.fasta)
    elif args.uniprot_accession:
        output_dir = ensure_output_dir(config["output_dir"])
        fetched_path = output_dir / f"{args.uniprot_accession}.fasta"
        fasta_text = fetch_uniprot_fasta(args.uniprot_accession)
        save_fasta_text(fasta_text, fetched_path)
        config["input_fasta"] = str(fetched_path.resolve())

    sequence = parse_fasta(config["input_fasta"])
    preprocess_result = preprocess_sequence(
        sequence,
        segment_length=int(config.get("long_sequence_segment_length", 1200)),
        overlap=int(config.get("long_sequence_overlap", 50)),
    )

    mhc_i = load_netmhc_predictions(config["netmhc_i_file"], epitope_class="MHC_I")
    mhc_ii = load_netmhc_predictions(config["netmhc_ii_file"], epitope_class="MHC_II")
    bcell = predict_bcell_epitopes(
        preprocess_result.sequence,
        window_size=int(config["simulation"].get("bcell_window", 14)),
        min_score=float(config["simulation"].get("bcell_min_score", 0.58)),
    )

    raw_epitopes = merge_epitope_predictions(mhc_i, mhc_ii, bcell)
    filtered_epitopes = filter_epitopes(raw_epitopes, thresholds=config["thresholds"])
    ranked_epitopes = score_epitopes(filtered_epitopes, weights=config["weights"])
    selected_epitopes = select_top_epitopes(
        ranked_epitopes,
        top_n=int(config["selection"].get("top_n", 10)),
        max_bcell_epitopes=int(config["selection"].get("max_bcell_epitopes", 3)),
    )
    vaccine_construct = build_vaccine_construct(selected_epitopes, construct_config=config["construct"])

    output_dir = ensure_output_dir(config["output_dir"])
    save_records(raw_epitopes, output_dir / "raw_epitopes.csv")
    save_records(filtered_epitopes, output_dir / "filtered_epitopes.csv")
    save_records(ranked_epitopes, output_dir / "ranked_epitopes.csv")
    save_records(selected_epitopes, output_dir / "selected_epitopes.csv")
    plot_score_distribution(ranked_epitopes, output_dir / "score_distribution.svg")
    plot_epitope_positions(raw_epitopes, len(preprocess_result.sequence), output_dir / "epitope_positions.svg")
    write_markdown_report(
        output_path=output_dir / "report.md",
        preprocess_summary={
            "sequence_length": len(preprocess_result.sequence),
            "removed_residues": preprocess_result.removed_residues,
            "antigenicity_prefilter": preprocess_result.antigenicity_prefilter,
            "segment_count": len(preprocess_result.segments),
        },
        raw_epitopes=raw_epitopes,
        filtered_epitopes=filtered_epitopes,
        ranked_epitopes=ranked_epitopes,
        vaccine_construct=vaccine_construct,
    )
    (output_dir / "vaccine_construct.fasta").write_text(
        f">multi_epitope_vaccine_construct\n{vaccine_construct}\n",
        encoding="utf-8",
    )

    print("Pipeline completed successfully.")
    print(f"Input FASTA: {config['input_fasta']}")
    print(f"Raw epitopes: {len(raw_epitopes)}")
    print(f"Filtered epitopes: {len(filtered_epitopes)}")
    print(f"Ranked epitopes: {len(ranked_epitopes)}")
    print(f"Vaccine construct length: {len(vaccine_construct)} aa")


if __name__ == "__main__":
    main()
