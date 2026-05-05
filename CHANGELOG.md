# Changelog

All notable changes to mbX Pro are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
the project adheres to [Semantic Versioning](https://semver.org/).

## [1.2.0] -- 2026-05-04

### Added
- **Pre-trained classifier download from Zenodo**
  ([https://zenodo.org/records/20021035](https://zenodo.org/records/20021035))
  in step 5 (`mbx_classifier_arranger.sh`).  Eight pre-trained, sha256-verified
  full-length Greengenes2 Naive-Bayes classifiers cover QIIME2 2023.2, 2023.5,
  2023.7, 2023.9, 2024.2 (sklearn 0.24.x) and 2024.5, 2024.10, 2025.4
  (sklearn 1.4.x).  When `CLASSIFIER_MODE=full-length`, the script downloads
  a sklearn-compatible classifier instead of training one locally, saving
  30-90 minutes per run.
- New run-info field `CLASSIFIER_SOURCE` with values `zenodo`, `cached`,
  `local-training`, `local-training-fallback`.  Persisted in
  `mbx_classifier_run_info.txt` and surfaced verbatim in the final report.
- New run-info fields `SCIKIT_LEARN_VERSION`, `SCIKIT_LEARN_FAMILY`,
  `ZENODO_RECORD_URL`, `ZENODO_QIIME2_USED`, `ZENODO_GG2_USED`,
  `ZENODO_FILENAME`, `ZENODO_SHA256_EXPECTED`, `ZENODO_SHA256_ACTUAL`,
  `ZENODO_NOTE` for full provenance auditing.
- `--skip-zenodo` flag on `mbx_classifier_arranger.sh` to force local
  training even when a compatible Zenodo classifier exists.
- `--assume-primer-length N` and `--primer-info FILE` flags on
  `mbx_dada2_parameter_finder.sh` for explicit overrides.

### Changed
- **DADA2 trim-left logic now honours `DETECTION_STATUS`**.  Previously,
  whenever no primer flags were passed, both `--p-trim-left-f` and
  `--p-trim-left-r` defaulted to 20 -- which silently lost 20 bp of real
  16S sequence when primers had been trimmed by the sequencing facility.
  The new rules in `mbx_dada2_parameter_finder.sh` are:
    | DETECTION_STATUS    | --p-trim-left-f / -r          |
    |---------------------|-------------------------------|
    | `DETECTED`          | primer length (auto-loaded)   |
    | `USER_SUPPLIED`     | primer length (auto-loaded)   |
    | `TRIMMED`           | **0** (was incorrectly 20)    |
    | `UNKNOWN` or absent | 20 (defensive default)        |
  `mbx_dada2_parameter_finder.sh` now auto-discovers
  `0_primer_handling/mbx_primer_info.txt` and applies the rule above.
  Explicit `--forward-primer` / `--reverse-primer` /
  `--assume-primer-length` flags continue to override auto-detection.
- `mbx_classifier_arranger.sh` skips the GG2 backbone download (~2 GB) when
  a Zenodo pre-trained classifier was successfully obtained and verified.
- `mbx_classifier_run.sh` skips both extract-reads (Step 2) AND
  fit-classifier-naive-bayes (Step 3) when `CLASSIFIER_SOURCE=zenodo|cached`,
  jumping straight to classify-sklearn.
- `mbx_classifier_arranger.sh` respects `DETECTION_STATUS=TRIMMED` even when
  primer sequences happen to be present in `mbx_primer_info.txt` -- a
  TRIMMED status is a strict signal to use the FULL backbone (no
  extract-reads).
- Final report now displays the new CLASSIFIER_SOURCE, the Zenodo file
  used (if any) with its sha256, and explains the trim-left rule per
  DETECTION_STATUS in plain English.

### Fixed
- Self-healing fall-back path in `mbx_classifier_run.sh`: if a downloaded
  Zenodo classifier fails to load with classify-sklearn (e.g. an
  unanticipated sklearn pickle incompatibility), the file is deleted, the
  GG2 backbone is downloaded automatically (via a one-shot
  `mbx_classifier_arranger.sh --skip-zenodo`), the classifier is retrained
  in place, and classify-sklearn is retried.  The pipeline NEVER aborts
  because of a Zenodo problem.
- Cosmetic: `mbx_dada2_run.sh` no longer adds the misleading
  `(default - no primer provided)` annotation when trim-left equals 20,
  because 20 is also a perfectly normal primer length.

## [1.1.0] -- 2026-04-30

### Added
- **Three-tier primer detection** in step 0 (`mbx_primer_identifier.sh`):
  TIER 1 direct match (Cutadapt-style IUPAC sliding window, all orientations),
  TIER 2 V-region anchor motif match (recognises "primers were already trimmed"
  and infers the V-region), TIER 3 rich failure diagnostics (read-length
  distribution, per-position base composition, top-3 5'-prefixes).
- New output fields in `mbx_primer_info.txt`: `DETECTION_STATUS`
  (`DETECTED` / `TRIMMED` / `UNKNOWN` / `USER_SUPPLIED`), `CONFIDENCE_LEVEL`,
  `INFERRED_REGION`, `DETECTION_NOTE`.
- **Auto-fallback to FULL-LENGTH classifier mode** in steps 5+6 when primers
  are unavailable. The pipeline is now self-healing: if the sequencing
  facility trimmed primers before delivery (or detection fails), step 5
  sets `CLASSIFIER_MODE=full-length` and step 6 trains Naive-Bayes directly
  on the entire Greengenes2 backbone -- no user intervention required.
- `--extra-primers` flag for user-supplied primer TSV.
- `--debug-csv` flag dumps the full primer x orientation rate matrix.
- `--report-best` flag emits a tentative best-candidate even below threshold.

### Changed
- Default primer-detection thresholds relaxed: `--mismatches 2 -> 3`,
  `--offset 15 -> 25` (handles Fluidigm CS1/CS2 linkers), `--min-rate
  0.10 -> 0.05`.
- Treats `N` as wild on BOTH read and primer sides (canonical IUPAC).
- Primer database is sanitised at load time (NFKC normalise + Cyrillic ->
  Latin homoglyph map + strict ASCII IUPAC enforcement). Fixes the silent
  Cyrillic-N (`U+041D`) bug that was lurking in the 338F / 341F entries.
- Step 6 (`mbx_classifier_run.sh`) no longer hard-fails when primers are
  None -- it auto-detects `CLASSIFIER_MODE=full-length` and trains on the
  full backbone instead.
- Final-report step 5 section now clearly states which classifier mode was
  used and explains the trade-offs.

### Fixed
- Pre-existing bug: orchestrator passed `--forward-primer` /
  `--reverse-primer` to `mbx_primer_identifier.sh` but the script did not
  accept those flags. v2 now does, and emits a `USER_SUPPLIED` info file.
- R1/R2 detection regex tightened (no longer mis-classifies filenames
  like `myr1_R2.fastq.gz`).
- `.fastq` (uncompressed) FASTQ files now supported alongside `.fastq.gz`.

## [1.0.0] -- 2026-04-30

### Added
- 18-step end-to-end pipeline (steps 0 through 17 of analysis + step 18 final report).
- Single-command orchestrator `mbXPro` (step 19 of the build) that invokes every
  step in the correct order using one shared output directory.
- Idempotent `--resume` mode that skips already-completed steps.
- Manual primer override (`--forward-primer` / `--reverse-primer`) for datasets
  whose primers were already trimmed by the sequencing facility.
- Self-contained HTML + A4 PDF final report with:
  - Embedded base64 figures (no external file dependencies)
  - Cross-step convergence table (taxa flagged by 2+ independent methods)
  - Per-step "not run" placeholder boxes that show the exact command needed
    to add that step to a future re-run
  - Inline mbX Pro logo + how-to-cite block
  - Print-optimised CSS that respects A4 page boundaries
- Installer (`install/install_mbXPro.sh`) and uninstaller
  (`install/uninstall_mbXPro.sh`).

### Notes
- macOS Bash 3.2 compatibility is preserved (no `mapfile`, no `declare -A`,
  no `${var^^}`).
- All R-based steps strip conda environment variables to avoid Rcpp.so /
  methods.dylib loading conflicts when QIIME2 is active.
- PICRUSt2 step automatically installs PICRUSt2 into its own conda environment
  on first use (does not pollute the QIIME2 environment).
