# Immunoinformatics Epitope-to-Vaccine Pipeline

## Problem Statement
This project implements a modular, research-oriented Python pipeline for transforming a protein sequence into a ranked epitope shortlist and a multi-epitope vaccine construct. The workflow supports MHC class I, MHC class II, and linear B-cell epitope streams, then applies biologically motivated filtering and weighted prioritization before assembling a final construct.

## Pipeline Architecture
The execution order is intentionally strict and mirrors an epitope discovery workflow used in computational vaccine design:

1. `src/input.py`
   Parses FASTA input and validates that only canonical amino acids are present.
2. `src/preprocess.py`
   Cleans the sequence, segments long proteins when needed, and computes a basic antigenicity pre-filter score.
3. `src/epitope_prediction.py`
   Loads NetMHCpan and NetMHCIIpan outputs from CSV files and simulates BepiPred-style linear B-cell epitopes using hydrophilicity, flexibility, and antigenicity-aware heuristics.
4. `src/filtering.py`
   Integrates VaxiJen, AllerTOP, and ToxinPred through deterministic heuristic surrogates when direct services are unavailable.
5. `src/scoring.py`
   Normalizes biological features and ranks epitopes using configurable weights from `config.yaml`.
6. `src/construct_builder.py`
   Selects top candidates and assembles a vaccine construct using `AAY` for CTL epitopes and `GPGPG` for HTL epitopes, with optional N-terminal adjuvanting.
7. `src/utils.py`
   Saves structured outputs, generates plots, and writes a Markdown report.
8. `main.py`
   Orchestrates the full pipeline end-to-end.

## Tool Integration Strategy
This implementation is ready for both real prediction outputs and simulation mode:

- `NetMHCpan`
  Expected as a CSV with `peptide`, `position`, `allele`, `ic50`, and `rank`.
- `NetMHCIIpan`
  Expected as a CSV with the same schema for class II predictions.
- `BepiPred`
  Simulated by scoring overlapping sequence windows for hydrophilicity, flexibility, and exposure-prone residue composition.
- `VaxiJen`
  Simulated using sequence-derived antigenicity logic that favors exposed, diverse, and moderately hydrophilic peptides.
- `AllerTOP`
  Simulated using aromaticity, cysteine content, repetitive motifs, and charge imbalance as allergenicity proxies.
- `ToxinPred`
  Simulated using hydrophobicity, basic residue enrichment, proline content, and glycine/serine relief as toxicity proxies.

These heuristics are deterministic and biologically motivated rather than random, which makes the outputs reproducible and interpretable for development and benchmarking.

## Project Structure
```text
immunoinformatics_pipeline/
├── data/
├── src/
│   ├── input.py
│   ├── preprocess.py
│   ├── epitope_prediction.py
│   ├── filtering.py
│   ├── scoring.py
│   ├── construct_builder.py
│   └── utils.py
├── notebooks/
├── results/
├── config.yaml
├── main.py
└── README.md
```

## Installation
No third-party Python packages are required for the bundled implementation. The project is intentionally written against the Python standard library so it can run in constrained research or sandbox environments.

## Usage
From the project root:

```bash
python main.py
```

You can also provide the sequence source explicitly:

```bash
python main.py --fasta data/sample_protein.fasta
python main.py --uniprot-accession P05067
```

Optional prediction/result overrides:

```bash
python main.py --fasta myprotein.fasta --netmhc-i data/my_netmhcpan.csv --netmhc-ii data/my_netmhciipan.csv --output-dir results_run_01
```

## Example Inputs
- FASTA: `data/sample_protein.fasta`
- MHC-I predictions: `data/sample_netmhcpan.csv`
- MHC-II predictions: `data/sample_netmhciipan.csv`

## Example Outputs
The pipeline writes the following into `results/`:

- `raw_epitopes.csv`
- `filtered_epitopes.csv`
- `ranked_epitopes.csv`
- `selected_epitopes.csv`
- `vaccine_construct.fasta`
- `score_distribution.svg`
- `epitope_positions.svg`
- `report.md`

Example outcome on the bundled sample dataset:

- Raw epitope candidates from all streams: `67`
- Filtered shortlist after antigenicity, toxicity, allergenicity, and binding cutoffs: `60`
- Ranked epitope table with normalized component scores: `60`
- Final vaccine construct length: `228 aa`

Top sample candidates from `results/ranked_epitopes.csv`:

1. `PERNECFLSHKDDSP` (`MHC_II`) with composite score `0.5517`
2. `HRFKDLGEE` (`MHC_I`) with composite score `0.4737`
3. `NECFLSHKD` (`MHC_I`) with composite score `0.4650`

## Weighted Ranking Model
The ranking score follows:

```text
Score =
  w1 * antigenicity
  + w2 * binding_affinity
  - w3 * toxicity
  - w4 * allergenicity
```

The implementation first min-max normalizes each feature. Binding affinity is inverted during normalization so that stronger binders receive higher contributions.

## Notes on Biological Realism
- Stronger predicted MHC binding improves rank.
- Higher antigenicity improves rank.
- Higher toxicity and allergenicity reduce rank.
- Filtering retains intermediate biological annotations to support review and threshold tuning.
- B-cell epitope simulation is deterministic and based on interpretable sequence properties.

## Input Strategy Recommendation
- Use a local FASTA file as the canonical, reproducible pipeline input.
- Use a UniProt accession as an optional convenience mode when you want to fetch a reference sequence automatically.
- When a UniProt accession is used, the pipeline downloads the FASTA record first and stores it locally in the output directory before continuing.

## Expected Output Snapshot
After running the sample data, inspect:

- `results/ranked_epitopes.csv` for the highest-priority candidates
- `results/vaccine_construct.fasta` for the assembled sequence
- `results/report.md` for a compact run summary

## Future Improvements
- Replace heuristic surrogates with live wrappers around local NetMHCpan, NetMHCIIpan, VaxiJen, AllerTOP, and ToxinPred installations.
- Add IEDB population coverage analysis for allele set prioritization.
- Add structural filtering using solvent accessibility or disorder predictions.
- Introduce ensemble learning for antigenicity or immunogenicity estimation.
- Add unit tests and workflow packaging for reproducible deployment.
