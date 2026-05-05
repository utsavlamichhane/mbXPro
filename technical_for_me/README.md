# technical_for_me/

Internal-only **per-script algorithmic dissection** for mbX Pro. One `.docx`
per pipeline script, each containing every detail needed to reason about,
debug, or modify that script — without any reference to a specific run's
file paths.

## Files

| # | Document | Script | Purpose |
|---|---|---|---|
| 1  | `1_mbx_primer_identifier.docx`                  | `mbx_primer_identifier.sh`               | 16S primer auto-detection |
| 2  | `2_create_manifest.docx`                        | `create_manifest.sh`                     | QIIME2 manifest |
| 3  | `3_artifact_creator.docx`                       | `artifact_creator.sh`                    | qiime tools import |
| 4  | `4_mbx_dada2_parameter_finder.docx`             | `mbx_dada2_parameter_finder.sh`          | DADA2 parameter recommender |
| 5  | `5_mbx_dada2_run.docx`                          | `mbx_dada2_run.sh`                       | DADA2 denoising |
| 6  | `6_mbx_classifier_arranger.docx`                | `mbx_classifier_arranger.sh`             | GG2 prep |
| 7  | `7_mbx_classifier_run.docx`                     | `mbx_classifier_run.sh`                  | NB classifier |
| 8  | `8_mbx_taxonomy_run.docx`                       | `mbx_taxonomy_run.sh`                    | Tax CSVs + filter |
| 9  | `9_mbx_ezclean_all_levels.docx`                 | `mbx_ezclean_all_levels.sh`              | mbX::ezclean per level |
| 10 | `10_mbx_ezviz_all_levels_all_treatments.docx`   | `mbx_ezviz_all_levels_all_treatments.sh` | Stacked-bar viz |
| 11 | `11_mbx_ezstat_all_levels_all_treatments.docx`  | `mbx_ezstat_all_levels_all_treatments.sh`| KW/Dunn/CLD |
| 12 | `12_mbx_pre_diversity_parameters.docx`          | `mbx_pre_diversity_parameters.sh`        | Tree + depth + analytical rarefaction |
| 13 | `13_mbx_alpha_diversity_run.docx`               | `mbx_alpha_diversity_run.sh`             | 5 alpha metrics + stats |
| 14 | `14_mbx_beta_diversity_run.docx`                | `mbx_beta_diversity_run.sh`              | Beta + PERMANOVA + Adonis |
| 15 | `15_mbx_ancombc2_run.docx`                      | `mbx_ancombc2_run.sh`                    | ANCOMBC2 DA |
| 16 | `16_mbx_picrust_run.docx`                       | `mbx_picrust_run.sh`                     | PICRUSt2 functional |
| 17 | `17_mbx_ml_classifier_run.docx`                 | `mbx_ml_classifier_run.sh`               | Random Forest biomarkers |
| 18 | `18_mbx_network_run.docx`                       | `mbx_network_run.sh`                     | Co-occurrence networks |
| 19 | `19_mbx_final_report.docx`                      | `mbx_final_report.sh`                    | HTML + PDF report |
| 20 | `20_mbXPro_orchestrator.docx`                   | `mbXPro`                                 | End-to-end orchestrator |

## Each document contains

1. **Purpose & scope** — what this script does, why it exists.
2. **Inputs and outputs** — exact data structures, formats, and downstream consumers.
3. **Complete parameter reference** — every flag and env var, with the *decision rule* used.
4. **Algorithm walkthrough** — line-by-line algorithmic intent, in narrative form.
5. **Flow diagram** — auto-generated PNG, embedded.
6. **Edge cases and caveats** — every weird thing we have seen in the wild.
7. **Implementation notes** — Bash 3.2 tricks, conda env handling, idempotency, etc.
8. **Testing checklist** — what to verify before considering changes production-ready.

## Regenerate

```bash
cd mbXPro/technical_for_me
python3 build_all.py
```

The build:

- generates one PNG flow diagram per chapter into `_diagrams/`,
- writes every `.docx` in this directory.

Requires `python-docx` and `matplotlib` on the system Python.

## Files in this directory

- `_flowchart.py` — matplotlib helper for rendering flow-diagram PNGs.
- `_build_technical.py` — chapters 1-3 + reusable docx helpers.
- `_chapters_part2.py` — chapters 4-12.
- `_chapters_part3.py` — chapters 13-20.
- `build_all.py` — convenience top-level runner.
- `_diagrams/` — every flow-diagram PNG (generated; ok to delete and rebuild).
- `<N>_<script_name>.docx` — the 20 deliverables.

## Privacy / publishing note

These docs reference **abstract** pipeline behaviour only — no real run-time
file paths, no test data, no machine-specific info. Safe to commit to a
public repo or share with collaborators, but they were written for the
maintainer (you).
