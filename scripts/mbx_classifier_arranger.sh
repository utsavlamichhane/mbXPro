#!/usr/bin/env bash
# =============================================================================
#  mbx_classifier_arranger.sh
#  Pre-step for taxonomic classification: try a pre-trained classifier from
#  Zenodo first; if not available, download GG2 + arrange for local training.
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  INPUT  : path to mbX_pro_outputs_<timestamp>/ directory
#           (auto-discovers 0_primer_handling/ and 4_dada2_outputs/ within it)
#
#  STEPS (high level):
#    1  Detect QIIME2 version + sklearn version  -> pick GG2 version
#    2  Discover required input files (primers, rep-seqs, feature table)
#    3  Decide CLASSIFIER_MODE: region-specific (primers) vs full-length (no primers)
#    4  IF CLASSIFIER_MODE=full-length: try Zenodo pre-trained classifier first.
#         (Saves 30-90 minutes of training; ~180 MB download.)
#         If sha256 verifies -> CLASSIFIER_SOURCE=zenodo, skip Steps 5 + 6.
#         If anything goes wrong -> CLASSIFIER_SOURCE=local-training, continue.
#    5  Download GG2 backbone.full-length.fna.qza + backbone.tax.qza
#         (skipped if a Zenodo classifier was successfully downloaded)
#    6  Export representative_sequences.qza -> dna-sequences.fasta
#    7  Compute min/max ASV lengths -> length_summary.txt
#    8  Write mbx_classifier_run_info.txt with all paths + commands.
#
#  CLASSIFIER_SOURCE values written to mbx_classifier_run_info.txt:
#    zenodo          - downloaded pre-trained Naive-Bayes from Zenodo
#                      (https://zenodo.org/records/20021035), verified by sha256
#    local-training  - will be trained locally by mbx_classifier_run.sh
#
#  OUTPUT STRUCTURE:
#    mbX_pro_outputs_<timestamp>/
#    └── 5_classifier_working_dir/
#        ├── exported_rep_seqs/
#        │   └── dna-sequences.fasta
#        ├── <gg2_ver>.backbone.full-length.fna.qza   (downloaded if needed)
#        ├── <gg2_ver>.backbone.tax.qza               (downloaded if needed)
#        ├── gg2_full_length_trained_classifier.qza   (Zenodo OR local)
#        ├── length_summary.txt
#        ├── MANIFEST.zenodo.tsv                      (if Zenodo succeeded)
#        └── mbx_classifier_run_info.txt
# =============================================================================

set -euo pipefail

# ── Helpers ───────────────────────────────────────────────────────────────────
err() {
  echo "" >&2
  echo "╔══════════════════════════════════════════════════════════════╗" >&2
  echo "║  ERROR                                                       ║" >&2
  echo "╚══════════════════════════════════════════════════════════════╝" >&2
  echo "[ERROR] $*" >&2
  echo "" >&2
  exit 1
}
warn() { echo "[WARN]  $*" >&2; }
info() { echo "[INFO]  $*"; }
ok()   { echo "[OK]    $*"; }
step() {
  echo ""
  echo "┌─────────────────────────────────────────────────────────────────"
  echo "│  Step $*"
  echo "└─────────────────────────────────────────────────────────────────"
}
sep()  { echo "────────────────────────────────────────────────────────────────"; }

_abspath() {
  if [[ -d "$1" ]]; then cd "$1" && pwd
  elif [[ -f "$1" ]]; then echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  else return 1; fi
}

# ── Zenodo classifier registry ────────────────────────────────────────────────
# Pre-trained full-length Greengenes2 Naive-Bayes classifiers, one per QIIME2
# version.  Hosted at https://zenodo.org/records/20021035
#
# Each row is whitespace-separated:
#    qiime2_ver  gg2_ver  filename  sha256
# (the same SHA-256 the Zenodo MANIFEST.tsv reports).
ZENODO_RECORD_URL="https://zenodo.org/records/20021035"
ZENODO_FILES_URL="${ZENODO_RECORD_URL}/files"
ZENODO_REGISTRY=$(cat <<'__ZENODO__'
2023.2   2022.10  gg2-2022.10-full-length-naive-bayes-qiime2-2023.2.qza   7392e97dfaf68dadb6068f651b23dd9e5c81c92145bd5361c7657bc9f65613ae
2023.5   2022.10  gg2-2022.10-full-length-naive-bayes-qiime2-2023.5.qza   90c7ed7fc2c60977725552e5ec4cab9a0d8c9937fcd03be21f9db103dc4ceb64
2023.7   2022.10  gg2-2022.10-full-length-naive-bayes-qiime2-2023.7.qza   0d1ce4515f56844c9404860dcdf91be17b9db6c368ec7eda72fca7db6fd7cae6
2023.9   2022.10  gg2-2022.10-full-length-naive-bayes-qiime2-2023.9.qza   c8fc36fd5384e2f358407e225c735203330295f4aa26a4ae4f269a7fe5b8887a
2024.2   2022.10  gg2-2022.10-full-length-naive-bayes-qiime2-2024.2.qza   3e1eb9875fa1f2322688a99baaf518bfe7ea99f2c7fa7e75995759f929889249
2024.5   2024.09  gg2-2024.09-full-length-naive-bayes-qiime2-2024.5.qza   5e6d9b63145422678c87c6a36fa39e2e97c269069b7238c0c0cff951bfb388a2
2024.10  2024.09  gg2-2024.09-full-length-naive-bayes-qiime2-2024.10.qza  3441db798acb6e2d2a64d057d227a771861c23409b10a076fc5eb2321546ea26
2025.4   2024.09  gg2-2024.09-full-length-naive-bayes-qiime2-2025.4.qza   612075d9354fecfff7a2513e46891b3d9b0dc79bbcaf29f78de6b3e5d7bff3f8
__ZENODO__
)

# Map QIIME2 version -> sklearn family.
# These are the sklearn versions actually used by each QIIME2 release.
# Pickle compatibility is roughly tied to MAJOR.MINOR (0.24 vs 1.4 here).
_zenodo_sklearn_for_qiime2() {
  case "$1" in
    2023.2|2023.5|2023.7|2023.9|2024.2) echo "0.24" ;;
    2024.5|2024.10|2025.4)              echo "1.4"  ;;
    *)                                  echo ""     ;;
  esac
}

# Trim sklearn version "1.4.2" or "0.24.1" to "1.4" / "0.24"
_zenodo_sklearn_trim() {
  echo "$1" | awk -F. '{print $1"."$2}'
}

# Compute sha256 of a file using the most-portable tool available.
# Echoes the lower-case hex digest, or empty string on failure.
_zenodo_sha256() {
  local f="$1"
  if command -v shasum &>/dev/null; then
    shasum -a 256 "$f" 2>/dev/null | awk '{print $1}'
  elif command -v sha256sum &>/dev/null; then
    sha256sum "$f" 2>/dev/null | awk '{print $1}'
  else
    echo ""
  fi
}

# Download a URL to a destination, with up to 3 retries and resume support.
# Returns 0 on success, 1 on failure.  Uses curl by default (universally
# available on macOS), with wget fallback.
_zenodo_curl() {
  local url="$1" dest="$2"
  if command -v curl &>/dev/null; then
    curl -L --retry 3 --retry-delay 2 --connect-timeout 15 --max-time 600 \
         --continue-at - --silent --show-error \
         -o "$dest" "$url"
  elif command -v wget &>/dev/null; then
    wget -c --no-check-certificate -q -O "$dest" "$url"
  else
    return 1
  fi
}

# Try to fetch the live MANIFEST.tsv from Zenodo (so we automatically pick up
# any newly-uploaded versions in future).  Returns 0 + the manifest text on
# stdout if successful, 1 if not.  On failure callers must use the embedded
# ZENODO_REGISTRY fallback.
_zenodo_fetch_manifest() {
  local out="$1"
  _zenodo_curl "${ZENODO_FILES_URL}/MANIFEST.tsv" "$out" || return 1
  [[ -s "$out" ]] || return 1
  head -1 "$out" | grep -q '^qiime2_version' || return 1
  return 0
}

# Print a 4-column registry "qiime2_ver gg2_ver filename sha256" parsed from
# either a live MANIFEST.tsv (col1=qiime2_version, col3=gg2_version,
# col4=filename, col6=sha256) or the embedded fallback.
_zenodo_load_registry() {
  local manifest_tsv="$1"
  if [[ -s "$manifest_tsv" ]] && head -1 "$manifest_tsv" | grep -q '^qiime2_version'; then
    awk -F'\t' 'NR>1 && $4 ~ /\.qza$/ {printf "%s\t%s\t%s\t%s\n",$1,$3,$4,$6}' \
        "$manifest_tsv"
  else
    printf '%s\n' "$ZENODO_REGISTRY" | awk 'NF==4 {printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'
  fi
}

# Choose the best entry from registry for a given (qiime2_ver, sklearn_family).
# Search order:
#   1. exact qiime2 version match (best, guaranteed pickle compatibility)
#   2. same sklearn family, closest qiime2 version (good compat, may fail)
# Echoes a single line "qiime2_ver gg2_ver filename sha256" or empty string.
_zenodo_pick_entry() {
  local registry="$1" want_q="$2" want_sk="$3"
  awk -v wq="$want_q" -v wsk="$want_sk" '
    BEGIN { best_exact=""; best_compat=""; best_compat_dist=1e9 }
    {
      q=$1; gg=$2; fn=$3; sh=$4
      if (q == wq) { best_exact = $0; next }
      sk = ""
      if (q ~ /^2023\.[2579]$|^2024\.2$/) sk = "0.24"
      else if (q ~ /^2024\.5$|^2024\.10$|^2025\.4$/) sk = "1.4"
      if (sk != "" && sk == wsk) {
        # closeness: numeric distance on YYYY.M
        split(q, a, "."); split(wq, b, ".")
        d = (a[1]*100 + a[2]) - (b[1]*100 + b[2]); if (d<0) d=-d
        if (d < best_compat_dist) { best_compat = $0; best_compat_dist = d }
      }
    }
    END { if (best_exact != "") print best_exact; else print best_compat }
  ' <<< "$registry"
}

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_classifier_arranger.sh — Download GG2, compute lengths, write classifier run info

USAGE:
  mbx_classifier_arranger.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Auto-discovers all required inputs from the mbX_pro_outputs directory:
    • 0_primer_handling/mbx_primer_info.txt     → forward/reverse primers
    • 4_dada2_outputs/representative_sequences.qza → ASV sequences for lengths

  Downloads the correct GG2 version based on your QIIME2 version:
    • QIIME2 2024.5+  → GG2 2024.09  (current, recommended)
    • QIIME2 < 2024.5 → GG2 2022.10  (legacy compatibility)

  Writes mbx_classifier_run_info.txt with ready-to-copy commands for:
    1. qiime feature-classifier extract-reads   ← region-specific trimming
    2. qiime feature-classifier fit-classifier-naive-bayes
    3. qiime feature-classifier classify-sklearn

  WHY extract-reads matters:
    Training the classifier on your specific primer-amplified V-region
    (rather than full-length 16S) substantially improves classification
    accuracy, especially at genus/species level. Always do this step.

OPTIONS:
  --skip-download      Skip GG2 download (if files already exist)
  --skip-zenodo        Skip the Zenodo pre-trained classifier check
                       (forces local training even when a compatible
                        pre-trained classifier is available online).
  --gg2-version <VER>  Override GG2 version (e.g. 2024.09 or 2022.10)
  --dry-run            Print commands without executing them
  -h, --help           Show this help and exit

ZENODO CLASSIFIERS:
  When CLASSIFIER_MODE = full-length (primers were already trimmed by the
  sequencing facility, or detection failed), the script tries to download a
  pre-trained Naive-Bayes classifier from
      https://zenodo.org/records/20021035
  matching your active QIIME2 version (8 versions are pre-built: QIIME2
  2023.2, 2023.5, 2023.7, 2023.9, 2024.2, 2024.5, 2024.10, 2025.4).

  If a sklearn-compatible match is found AND its sha256 verifies, we use it
  directly -- no GG2 download, no 30-90 minute training.

  If anything goes wrong (no match, network failure, sha256 mismatch,
  corrupt download), the pipeline silently falls back to local training:
  the GG2 backbone is downloaded and the classifier is trained on the fly,
  exactly as in mbX Pro <= 1.1.x.

EXAMPLES:
  mbx_classifier_arranger.sh /path/to/mbX_pro_outputs_20250422_143022
  mbx_classifier_arranger.sh /path/to/mbX_pro_outputs_20250422_143022 --skip-download
  mbx_classifier_arranger.sh /path/to/mbX_pro_outputs_20250422_143022 --skip-zenodo
  mbx_classifier_arranger.sh /path/to/mbX_pro_outputs_20250422_143022 --gg2-version 2022.10

OUTPUT FILES:
  5_classifier_working_dir/
  ├── exported_rep_seqs/
  │   └── dna-sequences.fasta
  ├── <gg2>.backbone.full-length.fna.qza
  ├── <gg2>.backbone.tax.qza
  ├── length_summary.txt
  └── mbx_classifier_run_info.txt        ← copy commands from here

COMMON ERRORS:
  "Insufficient disk space"
    → GG2 files need ~4 GB free. Free up space and re-run.
  "primer_info.txt not found"
    → Run mbx_primer_identifier.sh first, or create the file manually.
  "representative_sequences.qza not found"
    → Run mbx_dada2_run.sh first to generate DADA2 outputs.
  "wget: command not found"
    → Install wget: brew install wget
      Or use --skip-download and download the GG2 files manually from:
      http://ftp.microbio.me/greengenes_release/current/

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
SKIP_DOWNLOAD=false
SKIP_ZENODO=false
GG2_VER_OVERRIDE=""
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)        usage ;;
    --dry-run)        DRY_RUN=true; shift ;;
    --skip-download)  SKIP_DOWNLOAD=true; shift ;;
    --skip-zenodo)    SKIP_ZENODO=true; shift ;;
    --gg2-version)
      GG2_VER_OVERRIDE="${2:?Missing value for --gg2-version}"; shift 2 ;;
    -*)  err "Unknown option: '${1}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Unexpected extra argument: '${1}'"; fi ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]]  && err "No mbX_pro_outputs directory provided."
[[ -d "$MBX_OUT_DIR" ]]  || err "Directory does not exist: '${MBX_OUT_DIR}'"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

# ── Verify QIIME2 ─────────────────────────────────────────────────────────────
command -v qiime &>/dev/null || err "qiime not found.
  → conda activate qiime2-amplicon-2025.4"

# ─────────────────────────────────────────────────────────────────────────────
step "1/8 — Detect QIIME2 version → select GG2 version"
# ─────────────────────────────────────────────────────────────────────────────

QIIME_RAW="$(qiime info 2>/dev/null | grep -i 'QIIME 2 release' | grep -oE '[0-9]{4}\.[0-9]+' | head -1 || true)"
[[ -z "$QIIME_RAW" ]] && QIIME_RAW="$(qiime info 2>/dev/null | grep -i 'version' | grep -oE '[0-9]{4}\.[0-9]+' | head -1 || true)"
[[ -z "$QIIME_RAW" ]] && { warn "Could not detect QIIME2 version. Defaulting to GG2 2024.09."; QIIME_RAW="2025.4"; }

# Also detect scikit-learn version (governs Zenodo classifier compatibility).
# This runs inside the active QIIME2 conda env, so Python's sklearn is the one
# actually used by qiime feature-classifier.
SKLEARN_RAW=""
if command -v python &>/dev/null; then
  SKLEARN_RAW="$(python -c 'import sklearn; print(sklearn.__version__)' 2>/dev/null || true)"
fi
[[ -z "$SKLEARN_RAW" ]] && SKLEARN_RAW="unknown"
SKLEARN_FAMILY="$(_zenodo_sklearn_trim "$SKLEARN_RAW")"
[[ -z "$SKLEARN_FAMILY" || "$SKLEARN_FAMILY" == "unknown" ]] && \
  SKLEARN_FAMILY="$(_zenodo_sklearn_for_qiime2 "$QIIME_RAW")"

info "QIIME2 version     : $QIIME_RAW"
info "scikit-learn       : $SKLEARN_RAW (family: ${SKLEARN_FAMILY:-unknown})"

# GG2 version selection:
# 2024.09 NB classifiers were built with QIIME2 2024.5 / sklearn 1.4.2
# 2022.10 NB classifiers were built with sklearn 0.24.1
# Use 2024.09 for QIIME2 2024.5+, 2022.10 for older
if [[ -n "$GG2_VER_OVERRIDE" ]]; then
  GG2_VER="$GG2_VER_OVERRIDE"
  info "GG2 version        : $GG2_VER (user override)"
else
  QIIME_YEAR="$(echo "$QIIME_RAW" | cut -d. -f1)"
  QIIME_MONTH="$(echo "$QIIME_RAW" | cut -d. -f2)"
  # Compare: >= 2024.5 → use 2024.09
  if [[ "$QIIME_YEAR" -gt 2024 ]] || \
     [[ "$QIIME_YEAR" -eq 2024 && "$QIIME_MONTH" -ge 5 ]]; then
    GG2_VER="2024.09"
  else
    GG2_VER="2022.10"
  fi
  info "GG2 version        : $GG2_VER (auto-selected for QIIME2 $QIIME_RAW)"
fi

GG2_BASE_URL="http://ftp.microbio.me/greengenes_release/${GG2_VER}"
GG2_FNA_FILE="${GG2_VER}.backbone.full-length.fna.qza"
GG2_TAX_FILE="${GG2_VER}.backbone.tax.qza"
GG2_FNA_URL="${GG2_BASE_URL}/${GG2_FNA_FILE}"
GG2_TAX_URL="${GG2_BASE_URL}/${GG2_TAX_FILE}"

# ─────────────────────────────────────────────────────────────────────────────
step "2/8 — Discover required input files"
# ─────────────────────────────────────────────────────────────────────────────

# ── Primer info ───────────────────────────────────────────────────────────────
PRIMER_TXT="${MBX_OUT_DIR}/0_primer_handling/mbx_primer_info.txt"
if [[ ! -f "$PRIMER_TXT" ]]; then
  warn "mbx_primer_info.txt not found at: $PRIMER_TXT"
  warn "Primer information will be set to 'None' in the run info file."
  warn "→ Run mbx_primer_identifier.sh first, or use --forward-primer / --reverse-primer"
  FWD_PRIMER="None"
  REV_PRIMER="None"
  FWD_LEN="None"
  REV_LEN="None"
  TARGET_REGION="Unknown"
else
  FWD_PRIMER="$(grep '^FORWARD_PRIMER_SEQUENCE=' "$PRIMER_TXT" | cut -d= -f2)"
  REV_PRIMER="$(grep '^REVERSE_PRIMER_SEQUENCE='  "$PRIMER_TXT" | cut -d= -f2)"
  FWD_LEN="$(grep   '^FORWARD_PRIMER_LENGTH='     "$PRIMER_TXT" | cut -d= -f2)"
  REV_LEN="$(grep   '^REVERSE_PRIMER_LENGTH='      "$PRIMER_TXT" | cut -d= -f2)"
  TARGET_REGION="$(grep '^FORWARD_PRIMER_TARGET_REGION=' "$PRIMER_TXT" | cut -d= -f2)"
  FWD_NAME="$(grep '^FORWARD_PRIMER_NAME='        "$PRIMER_TXT" | cut -d= -f2)"
  REV_NAME="$(grep '^REVERSE_PRIMER_NAME='         "$PRIMER_TXT" | cut -d= -f2)"
  [[ -z "$FWD_PRIMER" ]] && FWD_PRIMER="None"
  [[ -z "$REV_PRIMER" ]] && REV_PRIMER="None"
  [[ -z "$FWD_LEN"    ]] && FWD_LEN="None"
  [[ -z "$REV_LEN"    ]] && REV_LEN="None"
  [[ -z "$TARGET_REGION" ]] && TARGET_REGION="Unknown"
  ok "Primers loaded: fwd=${FWD_NAME:-None}  rev=${REV_NAME:-None}"
fi

info "Forward primer     : $FWD_PRIMER (${FWD_LEN} bp)"
info "Reverse primer     : $REV_PRIMER (${REV_LEN} bp)"
info "Target region      : $TARGET_REGION"

# ── Representative sequences ──────────────────────────────────────────────────
REP_SEQS_QZA="${MBX_OUT_DIR}/4_dada2_outputs/representative_sequences.qza"
[[ -f "$REP_SEQS_QZA" ]] || err "representative_sequences.qza not found at:
  $REP_SEQS_QZA
  → Run mbx_dada2_run.sh first to generate DADA2 outputs."

ok "Rep seqs found: $REP_SEQS_QZA"

# Feature table (needed for next step - optional, just log if present)
FEATURE_TABLE_QZA="${MBX_OUT_DIR}/4_dada2_outputs/feature_table.qza"
[[ -f "$FEATURE_TABLE_QZA" ]] \
  && ok "Feature table found: $FEATURE_TABLE_QZA" \
  || warn "feature_table.qza not found — you will need it for classify-sklearn downstream."

# ─────────────────────────────────────────────────────────────────────────────
step "3/8 — Check disk space (GG2 needs ~4 GB free)"
# ─────────────────────────────────────────────────────────────────────────────

TARGET_DIR="$(dirname "$MBX_OUT_DIR")"
if command -v df &>/dev/null; then
  # df on macOS uses 512-byte blocks, on Linux uses 1K blocks — use -k for portability
  AVAIL_KB="$(df -k "$TARGET_DIR" 2>/dev/null | awk 'NR==2 {print $4}' || echo 0)"
  AVAIL_GB=$(( AVAIL_KB / 1024 / 1024 ))
  info "Available disk space: ~${AVAIL_GB} GB at $(dirname "$MBX_OUT_DIR")"
  if [[ "$AVAIL_GB" -lt 4 ]]; then
    warn "Less than 4 GB free. GG2 files may not fit."
    warn "  backbone.full-length.fna.qza ~ 1.5–2.5 GB"
    warn "  backbone.tax.qza             ~ 0.5–1.0 GB"
    warn "  Continue? The download will fail if space runs out."
  fi
else
  warn "Could not check disk space — proceed with caution."
fi

# ─────────────────────────────────────────────────────────────────────────────
step "4/8 — Create 5_classifier_working_dir/"
# ─────────────────────────────────────────────────────────────────────────────

CLASSIFIER_DIR="${MBX_OUT_DIR}/5_classifier_working_dir"
EXPORT_DIR="${CLASSIFIER_DIR}/exported_rep_seqs"
LENGTH_TXT="${CLASSIFIER_DIR}/length_summary.txt"
RUN_INFO="${CLASSIFIER_DIR}/mbx_classifier_run_info.txt"

mkdir -p "$CLASSIFIER_DIR" "$EXPORT_DIR" \
  || err "Could not create output directories — check permissions."

ok "Output directory   : $CLASSIFIER_DIR"

GG2_FNA_PATH="${CLASSIFIER_DIR}/${GG2_FNA_FILE}"
GG2_TAX_PATH="${CLASSIFIER_DIR}/${GG2_TAX_FILE}"

# ─────────────────────────────────────────────────────────────────────────────
step "5/9 — Decide classifier mode (region-specific vs full-length)"
# ─────────────────────────────────────────────────────────────────────────────
# This is the SINGLE source of truth for which mode the classifier runs in.
# Decided here so that step 6 knows whether to even consider Zenodo.

PRIMERS_AVAILABLE=false
[[ "$FWD_PRIMER" != "None" && "$REV_PRIMER" != "None" \
   && -n "$FWD_PRIMER"   && -n "$REV_PRIMER" ]] && PRIMERS_AVAILABLE=true

# Honour DETECTION_STATUS too — TRIMMED is the explicit signal that the
# primers are NOT in the reads, which forces full-length even if a stale
# FORWARD_PRIMER_SEQUENCE leaked through somehow.
if [[ -f "$PRIMER_TXT" ]]; then
  DET_STATUS_NOW="$(grep '^DETECTION_STATUS=' "$PRIMER_TXT" | cut -d= -f2 | tr -d '[:space:]')"
  case "$DET_STATUS_NOW" in
    TRIMMED|UNKNOWN|"")
      PRIMERS_AVAILABLE=false
      info "DETECTION_STATUS=${DET_STATUS_NOW:-(empty)} -> forcing CLASSIFIER_MODE=full-length"
      ;;
    DETECTED|USER_SUPPLIED)
      :    # PRIMERS_AVAILABLE keeps whatever the FWD_PRIMER/REV_PRIMER check decided
      ;;
  esac
fi

if $PRIMERS_AVAILABLE; then
  CLASSIFIER_MODE="region-specific"
  TRAINED_CLASSIFIER_QZA="${CLASSIFIER_DIR}/gg2_trained_classifier.qza"
else
  CLASSIFIER_MODE="full-length"
  TRAINED_CLASSIFIER_QZA="${CLASSIFIER_DIR}/gg2_full_length_trained_classifier.qza"
fi

info "Classifier mode    : $CLASSIFIER_MODE"
info "Trained classifier : $TRAINED_CLASSIFIER_QZA"

# Filled in by Step 6 below.  Default = "local-training" (the previous behaviour).
CLASSIFIER_SOURCE="local-training"
ZENODO_QIIME2_USED=""
ZENODO_GG2_USED=""
ZENODO_SHA256_EXPECTED=""
ZENODO_SHA256_ACTUAL=""
ZENODO_FILENAME=""
ZENODO_NOTE="(not attempted)"

# ─────────────────────────────────────────────────────────────────────────────
step "6/9 — Try pre-trained Zenodo classifier (full-length mode only)"
# ─────────────────────────────────────────────────────────────────────────────
# When CLASSIFIER_MODE=full-length we can replace the local 30-90 min training
# step with a 1-3 min download from https://zenodo.org/records/20021035 .
#
# This block is intentionally defensive:
#   • on ANY failure (no match, network, sha256 mismatch, corrupt file)
#     we silently fall through to local-training mode.
#   • the pipeline NEVER aborts because of Zenodo issues.

if [[ "$CLASSIFIER_MODE" != "full-length" ]]; then
  info "CLASSIFIER_MODE=region-specific -> Zenodo classifiers do not apply."
  info "Will train locally on V-region-extracted reads (steps 7-9)."
  ZENODO_NOTE="skipped (region-specific mode)"
elif $SKIP_ZENODO; then
  warn "--skip-zenodo set: skipping Zenodo lookup, will train locally."
  ZENODO_NOTE="skipped (--skip-zenodo)"
elif [[ -f "$TRAINED_CLASSIFIER_QZA" ]]; then
  ok "A trained classifier already exists at: $TRAINED_CLASSIFIER_QZA"
  ok "Re-using it (delete the file to force a fresh download / training)."
  CLASSIFIER_SOURCE="cached"
  ZENODO_NOTE="not attempted (cached file already present)"
elif $DRY_RUN; then
  warn "[DRY-RUN] Would try Zenodo download for QIIME2 ${QIIME_RAW} (sklearn ${SKLEARN_FAMILY})."
  ZENODO_NOTE="dry-run"
else
  # 1) Try to fetch the live MANIFEST.tsv from Zenodo for up-to-date entries.
  ZENODO_MANIFEST="${CLASSIFIER_DIR}/MANIFEST.zenodo.tsv"
  info "Fetching Zenodo manifest: ${ZENODO_FILES_URL}/MANIFEST.tsv"
  if _zenodo_fetch_manifest "$ZENODO_MANIFEST"; then
    ok "Zenodo manifest downloaded -> $(basename "$ZENODO_MANIFEST")"
  else
    warn "Could not download Zenodo manifest (network issue?). Using built-in registry."
    rm -f "$ZENODO_MANIFEST" 2>/dev/null || true
  fi

  REGISTRY="$(_zenodo_load_registry "$ZENODO_MANIFEST")"
  if [[ -z "$REGISTRY" ]]; then
    warn "Zenodo registry empty after parse. Falling back to local training."
    ZENODO_NOTE="failed (empty registry)"
  else
    ENTRY="$(_zenodo_pick_entry "$REGISTRY" "$QIIME_RAW" "$SKLEARN_FAMILY")"
    if [[ -z "$ENTRY" ]]; then
      warn "No Zenodo classifier compatible with QIIME2 ${QIIME_RAW} / sklearn ${SKLEARN_FAMILY}."
      warn "Will fall back to local training."
      ZENODO_NOTE="no compatible match"
    else
      Z_QIIME="$(echo "$ENTRY" | awk -F'\t' '{print $1}')"
      Z_GG2="$(echo "$ENTRY" | awk -F'\t' '{print $2}')"
      Z_FILE="$(echo "$ENTRY" | awk -F'\t' '{print $3}')"
      Z_SHA="$(echo "$ENTRY" | awk -F'\t' '{print $4}')"
      Z_URL="${ZENODO_FILES_URL}/${Z_FILE}"

      if [[ "$Z_QIIME" == "$QIIME_RAW" ]]; then
        info "Exact QIIME2 match: ${Z_QIIME} -> ${Z_FILE}"
      else
        warn "No exact QIIME2 ${QIIME_RAW} match; closest sklearn-compatible:"
        warn "  picked ${Z_QIIME} (sklearn family ${SKLEARN_FAMILY}) -> ${Z_FILE}"
        warn "  This may still load; if not, we will fall back to local training."
      fi
      info "Expected SHA-256: ${Z_SHA}"
      info "URL             : ${Z_URL}"

      Z_TMP="${TRAINED_CLASSIFIER_QZA}.zenodo_download.part"
      rm -f "$Z_TMP"
      info "Downloading pre-trained classifier (~180 MB)..."
      if ! _zenodo_curl "$Z_URL" "$Z_TMP"; then
        warn "Download failed for ${Z_URL}. Falling back to local training."
        rm -f "$Z_TMP"
        ZENODO_NOTE="download failed"
      else
        # Verify size + sha256
        ACTUAL_SHA="$(_zenodo_sha256 "$Z_TMP")"
        if [[ -z "$ACTUAL_SHA" ]]; then
          warn "Could not compute SHA-256 (no shasum/sha256sum). Skipping integrity check."
          warn "Continuing on file-size check only."
        fi
        if [[ -n "$ACTUAL_SHA" && "$ACTUAL_SHA" != "$Z_SHA" ]]; then
          warn "SHA-256 mismatch!"
          warn "  expected: ${Z_SHA}"
          warn "  got     : ${ACTUAL_SHA}"
          warn "Falling back to local training."
          rm -f "$Z_TMP"
          ZENODO_NOTE="sha256 mismatch"
        elif ! qiime tools peek "$Z_TMP" 2>/dev/null | grep -q 'TaxonomicClassifier'; then
          warn "Downloaded file does not look like a TaxonomicClassifier .qza."
          warn "Falling back to local training."
          rm -f "$Z_TMP"
          ZENODO_NOTE="invalid artifact type"
        else
          mv -f "$Z_TMP" "$TRAINED_CLASSIFIER_QZA"
          ok "Pre-trained classifier verified -> $(basename "$TRAINED_CLASSIFIER_QZA")"
          CLASSIFIER_SOURCE="zenodo"
          ZENODO_QIIME2_USED="$Z_QIIME"
          ZENODO_GG2_USED="$Z_GG2"
          ZENODO_FILENAME="$Z_FILE"
          ZENODO_SHA256_EXPECTED="$Z_SHA"
          ZENODO_SHA256_ACTUAL="${ACTUAL_SHA:-(skipped)}"
          ZENODO_NOTE="OK (downloaded, sha256 verified)"
          # If we picked an inexact match, override GG2_VER so the rest of
          # the script reflects the actual training reference for the report.
          if [[ "$Z_GG2" != "$GG2_VER" ]]; then
            info "Adjusting GG2_VERSION to match downloaded classifier: $GG2_VER -> $Z_GG2"
            GG2_VER="$Z_GG2"
            GG2_FNA_FILE="${GG2_VER}.backbone.full-length.fna.qza"
            GG2_TAX_FILE="${GG2_VER}.backbone.tax.qza"
            GG2_FNA_PATH="${CLASSIFIER_DIR}/${GG2_FNA_FILE}"
            GG2_TAX_PATH="${CLASSIFIER_DIR}/${GG2_TAX_FILE}"
            GG2_BASE_URL="http://ftp.microbio.me/greengenes_release/${GG2_VER}"
            GG2_FNA_URL="${GG2_BASE_URL}/${GG2_FNA_FILE}"
            GG2_TAX_URL="${GG2_BASE_URL}/${GG2_TAX_FILE}"
          fi
        fi
      fi
    fi
  fi
fi

case "$CLASSIFIER_SOURCE" in
  zenodo)
    ok "CLASSIFIER_SOURCE=zenodo  (no local training will be needed)" ;;
  cached)
    ok "CLASSIFIER_SOURCE=cached  (existing classifier will be reused)" ;;
  *)
    info "CLASSIFIER_SOURCE=local-training  (the classifier will be trained on the fly)" ;;
esac

# ─────────────────────────────────────────────────────────────────────────────
step "7/9 — Download GG2 reference files (skipped if Zenodo classifier OK)"
# ─────────────────────────────────────────────────────────────────────────────

_download() {
  local url="$1" dest="$2" label="$3"
  if [[ -f "$dest" ]]; then
    ok "Already exists — skipping: $label"
    return 0
  fi
  info "Downloading $label ..."
  info "URL: $url"
  if $DRY_RUN; then
    warn "[DRY-RUN] Would run: wget -c --no-check-certificate -O \"$dest\" \"$url\""
    return 0
  fi
  # wget -c enables resume; falls back to curl if wget not available
  if command -v wget &>/dev/null; then
    wget -c --no-check-certificate --show-progress \
      -O "$dest" "$url" \
      || err "Download failed: $url
  → Check your internet connection.
  → If wget is missing: brew install wget
  → Or download manually and place at: $dest"
  elif command -v curl &>/dev/null; then
    curl -L --retry 3 --continue-at - \
      -o "$dest" "$url" \
      || err "Download failed: $url
  → Check your internet connection.
  → Or download manually and place at: $dest"
  else
    err "Neither wget nor curl found.
  → brew install wget
  → Or download manually:
      $url
    and place the file at: $dest"
  fi
  ok "Downloaded: $label"
}

# Skip GG2 download in two cases:
#   1. CLASSIFIER_SOURCE=zenodo|cached AND CLASSIFIER_MODE=full-length
#      -> we already have the trained classifier; the GG2 backbone is unused
#         downstream (the report still reads the GG2 version from the run info).
#   2. --skip-download (user explicit override)
if $SKIP_DOWNLOAD; then
  warn "--skip-download set. Skipping GG2 download."
  [[ -f "$GG2_FNA_PATH" ]] || warn "Expected file missing: $GG2_FNA_PATH"
  [[ -f "$GG2_TAX_PATH" ]] || warn "Expected file missing: $GG2_TAX_PATH"
elif [[ "$CLASSIFIER_MODE" == "full-length" ]] \
     && { [[ "$CLASSIFIER_SOURCE" == "zenodo" ]] || [[ "$CLASSIFIER_SOURCE" == "cached" ]]; }; then
  ok "Skipping GG2 download — pre-trained classifier already available"
  ok "  ($CLASSIFIER_SOURCE).  GG2 reference files are NOT needed."
  info "  This saves ~2 GB of disk space and ~5-15 minutes of download."
else
  _download "$GG2_FNA_URL" "$GG2_FNA_PATH" "${GG2_FNA_FILE}"
  _download "$GG2_TAX_URL" "$GG2_TAX_PATH" "${GG2_TAX_FILE}"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "8/9 — Export representative sequences from QZA"
# ─────────────────────────────────────────────────────────────────────────────

FASTA_PATH="${EXPORT_DIR}/dna-sequences.fasta"

if [[ -f "$FASTA_PATH" ]]; then
  ok "FASTA already exported — skipping."
elif $DRY_RUN; then
  warn "[DRY-RUN] Would export rep seqs to $EXPORT_DIR"
else
  info "Exporting representative_sequences.qza..."
  qiime tools export \
    --input-path  "$REP_SEQS_QZA" \
    --output-path "$EXPORT_DIR" \
    || err "qiime tools export failed on representative_sequences.qza.
  → Verify the artifact is valid: qiime tools peek $REP_SEQS_QZA"
  [[ -f "$FASTA_PATH" ]] || err "Export completed but dna-sequences.fasta not found at: $FASTA_PATH"
  ok "Exported: $FASTA_PATH"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "9a/9 — Compute min/max ASV lengths"
# ─────────────────────────────────────────────────────────────────────────────

if $DRY_RUN; then
  warn "[DRY-RUN] Would compute lengths from $FASTA_PATH"
  MIN_LEN="N/A"
  MAX_LEN="N/A"
  TOTAL_ASVS="N/A"
elif [[ ! -f "$FASTA_PATH" ]]; then
  warn "FASTA not found — skipping length computation."
  MIN_LEN="N/A"; MAX_LEN="N/A"; TOTAL_ASVS="N/A"
else
  info "Computing ASV lengths from dna-sequences.fasta..."

  # User's awk command (exact, as provided) — saves length_summary.txt
  awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq$0} END {print length(seq)}' \
    "$FASTA_PATH" | sort -n | \
    awk 'NR==1 {min=$1} {max=$1; count++} END {
      printf "Min Length: %d\nMax Length: %d\nTotal ASVs: %d\n", min, max, count
    }' > "$LENGTH_TXT"

  MIN_LEN="$(grep 'Min' "$LENGTH_TXT" | awk '{print $NF}')"
  MAX_LEN="$(grep 'Max' "$LENGTH_TXT" | awk '{print $NF}')"
  TOTAL_ASVS="$(grep 'Total' "$LENGTH_TXT" | awk '{print $NF}')"

  ok "Length summary:"
  ok "  Min length : $MIN_LEN bp"
  ok "  Max length : $MAX_LEN bp"
  ok "  Total ASVs : $TOTAL_ASVS"
fi

# ── Auto-detect CPU count for classify-sklearn ────────────────────────────────
if command -v nproc &>/dev/null; then
  N_JOBS="$(nproc)"
elif command -v sysctl &>/dev/null; then
  N_JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
else
  N_JOBS=1
fi

# ─────────────────────────────────────────────────────────────────────────────
step "9b/9 — Write mbx_classifier_run_info.txt"
# ─────────────────────────────────────────────────────────────────────────────

# CLASSIFIER_MODE, PRIMERS_AVAILABLE, TRAINED_CLASSIFIER_QZA, CLASSIFIER_SOURCE
# were all decided in Step 5 and (for SOURCE) Step 6 above.

if [[ "$CLASSIFIER_MODE" == "full-length" && "$CLASSIFIER_SOURCE" == "local-training" ]]; then
  warn "═══════════════════════════════════════════════════════════════════"
  warn "  No usable primers AND no Zenodo classifier available."
  warn "  → Falling back to LOCAL FULL-LENGTH training."
  warn "      • extract-reads will be SKIPPED (no primers to trim with)."
  warn "      • Naive-Bayes will be trained on the entire GG2 backbone."
  warn "      • Expect 30-90 min training the first time, 8-16 GB RAM."
  warn "      • The trained classifier is reusable for ANY future project."
  warn "      • Genus-level accuracy is virtually identical to region-specific."
  warn "═══════════════════════════════════════════════════════════════════"
fi

# Output paths for downstream commands
EXTRACTED_READS_QZA="${CLASSIFIER_DIR}/gg2_extracted_reads.qza"
TAXONOMY_QZA="${MBX_OUT_DIR}/5_classifier_working_dir/taxonomy.qza"
TAXONOMY_QZV="${MBX_OUT_DIR}/5_classifier_working_dir/taxonomy.qzv"

NOW="$(date '+%Y-%m-%d %H:%M:%S')"

cat > "$RUN_INFO" << RUNINFO
# ============================================================================
# mbx_classifier_run_info.txt
# Generated by mbx_classifier_arranger.sh
# Date              : $NOW
# ============================================================================
#
# This file contains all the information and ready-to-run QIIME2 commands
# for taxonomic classification using Greengenes2.
#
# WORKFLOW OVERVIEW:
#   Step A → extract-reads  (trim GG2 to your amplicon region using primers)
#   Step B → fit-classifier (train a Naive Bayes classifier on the trimmed ref)
#   Step C → classify-sklearn (assign taxonomy to your ASVs)
#   Step D → tabulate (make taxonomy viewable)
#
# WHY EXTRACT-READS IS ESSENTIAL:
#   Training on your specific V-region (e.g. V3-V4) rather than full-length 16S
#   substantially improves genus/species-level accuracy. Always run Step A+B
#   unless you are using the pre-built V4 classifier for EMP 515F/806R data.
# ============================================================================

# ── Run environment ──────────────────────────────────────────────────────────
QIIME2_VERSION=$QIIME_RAW
GG2_VERSION=$GG2_VER
TIMESTAMP=$NOW
N_CPU_CORES=$N_JOBS

# ── Classifier mode ──────────────────────────────────────────────────────────
# region-specific : primers were available → extract-reads → train on V-region
# full-length     : primers were missing/trimmed → train directly on the GG2
#                   backbone. extract-reads is skipped; classify-sklearn uses
#                   the resulting classifier as-is.
CLASSIFIER_MODE=$CLASSIFIER_MODE
PRIMERS_AVAILABLE=$PRIMERS_AVAILABLE

# ── Classifier source ────────────────────────────────────────────────────────
# zenodo          : pre-trained classifier downloaded from
#                   https://zenodo.org/records/20021035  (sha256-verified).
# cached          : a previous run already produced this trained classifier.
# local-training  : will be trained on the fly by mbx_classifier_run.sh.
CLASSIFIER_SOURCE=$CLASSIFIER_SOURCE
SCIKIT_LEARN_VERSION=$SKLEARN_RAW
SCIKIT_LEARN_FAMILY=$SKLEARN_FAMILY
ZENODO_RECORD_URL=$ZENODO_RECORD_URL
ZENODO_QIIME2_USED=${ZENODO_QIIME2_USED:-(none)}
ZENODO_GG2_USED=${ZENODO_GG2_USED:-(none)}
ZENODO_FILENAME=${ZENODO_FILENAME:-(none)}
ZENODO_SHA256_EXPECTED=${ZENODO_SHA256_EXPECTED:-(none)}
ZENODO_SHA256_ACTUAL=${ZENODO_SHA256_ACTUAL:-(none)}
ZENODO_NOTE=$ZENODO_NOTE

# ── Input files ──────────────────────────────────────────────────────────────
REPRESENTATIVE_SEQUENCES_QZA=$REP_SEQS_QZA
FEATURE_TABLE_QZA=$FEATURE_TABLE_QZA
PRIMER_INFO_TXT=$PRIMER_TXT

# ── GG2 reference files (downloaded — empty when CLASSIFIER_SOURCE=zenodo) ────
GG2_FULL_LENGTH_FNA_QZA=$GG2_FNA_PATH
GG2_TAXONOMY_QZA=$GG2_TAX_PATH
GG2_DOWNLOAD_URL_FNA=$GG2_FNA_URL
GG2_DOWNLOAD_URL_TAX=$GG2_TAX_URL

# ── Primer information ───────────────────────────────────────────────────────
TARGET_REGION=$TARGET_REGION
FORWARD_PRIMER_NAME=${FWD_NAME:-None}
FORWARD_PRIMER_SEQUENCE=$FWD_PRIMER
FORWARD_PRIMER_LENGTH=${FWD_LEN} bp
REVERSE_PRIMER_NAME=${REV_NAME:-None}
REVERSE_PRIMER_SEQUENCE=$REV_PRIMER
REVERSE_PRIMER_LENGTH=${REV_LEN} bp

# ── ASV length statistics ────────────────────────────────────────────────────
MIN_ASV_LENGTH=$MIN_LEN bp
MAX_ASV_LENGTH=$MAX_LEN bp
TOTAL_ASVS=$TOTAL_ASVS
FASTA_PATH=$FASTA_PATH
LENGTH_SUMMARY_TXT=$LENGTH_TXT

# ── Output paths for classifier commands ─────────────────────────────────────
EXTRACTED_READS_QZA=$EXTRACTED_READS_QZA
TRAINED_CLASSIFIER_QZA=$TRAINED_CLASSIFIER_QZA
TAXONOMY_QZA=$TAXONOMY_QZA
TAXONOMY_QZV=$TAXONOMY_QZV

RUNINFO

# ── Append ready-to-run commands ──────────────────────────────────────────────
cat >> "$RUN_INFO" << 'CMDS_HEADER'

# ============================================================================
# READY-TO-RUN COMMANDS
# Copy and paste these into your terminal (inside your QIIME2 conda env)
# ============================================================================

CMDS_HEADER

if $PRIMERS_AVAILABLE; then
  cat >> "$RUN_INFO" << STEP_A
# ── Step A: Extract reads from GG2 using your primers + ASV length range ─────
# This trims the full-length GG2 reference to your specific amplicon region.
# Using the actual min/max lengths from your data improves accuracy significantly.
# --p-n-jobs uses all available CPU cores for speed.

qiime feature-classifier extract-reads \\
  --i-sequences  $GG2_FNA_PATH \\
  --p-f-primer   $FWD_PRIMER \\
  --p-r-primer   $REV_PRIMER \\
  --p-min-length $MIN_LEN \\
  --p-max-length $MAX_LEN \\
  --p-n-jobs     $N_JOBS \\
  --o-reads      $EXTRACTED_READS_QZA

STEP_A
else
  cat >> "$RUN_INFO" << STEP_A_NONE
# ── Step A: Extract reads — AUTO-SKIPPED (FULL-LENGTH classifier mode) ────────
# Primers were not available (already trimmed by sequencing facility, or detection
# failed).  The pipeline has automatically switched to FULL-LENGTH classifier mode:
#   • Step A (extract-reads) is SKIPPED entirely.
#   • Step B trains Naive-Bayes directly on the FULL GG2 backbone.
#   • Step C classifies your ASVs against that classifier (same as usual).
#
# Trade-off:
#   • Pro:  fully automatic — no user intervention, no missing primers needed.
#   • Pro:  trained classifier (gg2_full_length_trained_classifier.qza) is reusable
#           across ANY future project, regardless of primer set.
#   • Con:  ~10-30% longer training time and ~2-3x larger classifier file.
#   • Con:  marginal precision loss at species level (V-region-trained is slightly
#           more discriminating for short reads).  Genus-level accuracy is virtually
#           identical.
#
# To force REGION-SPECIFIC mode instead, re-run:
#   mbx_primer_identifier.sh /path/to/FASTQ --forward-primer <SEQ> --reverse-primer <SEQ>
#   mbx_classifier_arranger.sh /path/to/mbX_pro_outputs_*

STEP_A_NONE
fi

if [[ "$CLASSIFIER_SOURCE" == "zenodo" || "$CLASSIFIER_SOURCE" == "cached" ]]; then
  cat >> "$RUN_INFO" << STEP_B_ZENODO
# ── Step B: Train Naive Bayes classifier — AUTO-SKIPPED (CLASSIFIER_SOURCE=$CLASSIFIER_SOURCE) ──
# A pre-trained Naive-Bayes classifier was obtained from
#   $ZENODO_RECORD_URL
# (file: $ZENODO_FILENAME, sha256-verified).
# This saves 30-90 minutes of local training time.
#
# The classifier was trained with QIIME2 ${ZENODO_QIIME2_USED} on the
# Greengenes2 ${ZENODO_GG2_USED} full-length backbone (sklearn family ${SKLEARN_FAMILY}).
#
# If classify-sklearn (Step C) fails to load this file (e.g. due to an
# unexpected sklearn pickle incompatibility), mbx_classifier_run.sh will
# DELETE this file and retry by training locally on the full GG2 backbone.

# Step B is intentionally NOT executed.

STEP_B_ZENODO
elif $PRIMERS_AVAILABLE; then
  cat >> "$RUN_INFO" << STEP_B_LOCAL
# ── Step B: Train Naive Bayes classifier (region-specific mode, local) ───────
# Trained on: extracted V-region reads
# This step can take 20–60 minutes.

qiime feature-classifier fit-classifier-naive-bayes \\
  --i-reference-reads     $EXTRACTED_READS_QZA \\
  --i-reference-taxonomy  $GG2_TAX_PATH \\
  --o-classifier          $TRAINED_CLASSIFIER_QZA

STEP_B_LOCAL
else
  cat >> "$RUN_INFO" << STEP_B_FALLBACK
# ── Step B: Train Naive Bayes classifier (full-length mode, local fallback) ──
# Trained on: full GG2 backbone (no primers available, no Zenodo classifier)
# This step can take 30–90 minutes.

qiime feature-classifier fit-classifier-naive-bayes \\
  --i-reference-reads     $GG2_FNA_PATH \\
  --i-reference-taxonomy  $GG2_TAX_PATH \\
  --o-classifier          $TRAINED_CLASSIFIER_QZA

STEP_B_FALLBACK
fi

cat >> "$RUN_INFO" << STEP_BC
# ── Step C: Classify your ASVs using the trained classifier ──────────────────
# --p-n-jobs uses all $N_JOBS available CPU cores.
# --p-confidence 0.7 is the QIIME2 default; lower to 0.6 for more assignments,
# raise to 0.8 for stricter, higher-confidence assignments.

qiime feature-classifier classify-sklearn \\
  --i-classifier   $TRAINED_CLASSIFIER_QZA \\
  --i-reads        $REP_SEQS_QZA \\
  --p-n-jobs       $N_JOBS \\
  --p-confidence   0.7 \\
  --o-classification $TAXONOMY_QZA

# ── Step D: Tabulate taxonomy → make it viewable ──────────────────────────────

qiime metadata tabulate \\
  --m-input-file  $TAXONOMY_QZA \\
  --o-visualization $TAXONOMY_QZV

# ── Step E: Taxonomy bar plots (requires metadata) ────────────────────────────

qiime taxa barplot \\
  --i-table          $FEATURE_TABLE_QZA \\
  --i-taxonomy       $TAXONOMY_QZA \\
  --m-metadata-file  /path/to/metadata.txt \\
  --o-visualization  ${CLASSIFIER_DIR}/taxa_barplot.qzv

# ── View results ──────────────────────────────────────────────────────────────
# Drag any .qzv file into: https://view.qiime2.org

# ── Expected output structure ─────────────────────────────────────────────────
# 5_classifier_working_dir/
# ├── exported_rep_seqs/
# │   └── dna-sequences.fasta
# ├── ${GG2_FNA_FILE}
# ├── ${GG2_TAX_FILE}
# ├── gg2_extracted_reads.qza        ← Step A output
# ├── gg2_trained_classifier.qza     ← Step B output (reusable!)
# ├── taxonomy.qza                   ← Step C output
# ├── taxonomy.qzv                   ← Step D output  (view at view.qiime2.org)
# ├── taxa_barplot.qzv               ← Step E output
# ├── length_summary.txt
# └── mbx_classifier_run_info.txt    ← this file

STEP_BC

ok "Run info written → $RUN_INFO"

# ── Final summary ─────────────────────────────────────────────────────────────
sep
ok "mbx_classifier_arranger.sh complete!"
sep
echo ""
echo "  ── Summary ──────────────────────────────────────────────────────"
echo "  QIIME2 version     : $QIIME_RAW"
echo "  scikit-learn       : $SKLEARN_RAW (family $SKLEARN_FAMILY)"
echo "  GG2 version        : $GG2_VER"
echo "  Classifier mode    : $CLASSIFIER_MODE"
echo "  Classifier source  : $CLASSIFIER_SOURCE"
case "$CLASSIFIER_SOURCE" in
  zenodo)
    echo "    Zenodo file      : $ZENODO_FILENAME"
    echo "    Trained QIIME2   : $ZENODO_QIIME2_USED"
    echo "    SHA-256          : ${ZENODO_SHA256_ACTUAL:0:16}... (verified)"
    ;;
  cached)
    echo "    (re-using a previously produced classifier)"
    ;;
esac
echo "  Forward primer     : $FWD_PRIMER"
echo "  Reverse primer     : $REV_PRIMER"
echo "  Target region      : $TARGET_REGION"
echo "  Min ASV length     : $MIN_LEN bp"
echo "  Max ASV length     : $MAX_LEN bp"
echo "  Total ASVs         : $TOTAL_ASVS"
echo "  CPU cores detected : $N_JOBS"
if [[ "$CLASSIFIER_MODE" == "full-length" ]]; then
  echo ""
  if [[ "$CLASSIFIER_SOURCE" == "zenodo" || "$CLASSIFIER_SOURCE" == "cached" ]]; then
    echo "  ✓ Full-length mode + pre-trained classifier (FAST PATH):"
    echo "      - extract-reads will be skipped"
    echo "      - fit-classifier-naive-bayes will be skipped (using Zenodo .qza)"
    echo "      - The pipeline will jump straight to classify-sklearn."
    echo "      - Estimated time saved: 30-90 minutes."
  else
    echo "  ⚠ Full-length mode active:"
    echo "      - extract-reads will be skipped"
    echo "      - NB will be trained directly on the full GG2 backbone"
    echo "      - This is automatic; the pipeline will continue to step 6 normally."
  fi
fi
echo ""
echo "  ── Output structure ─────────────────────────────────────────────"
echo "  $MBX_OUT_DIR/"
echo "  ├── 0_primer_handling/"
echo "  ├── 1_manifest_file/"
echo "  ├── 2_first_artifact_file/"
echo "  ├── 3_dada2_parameters/"
echo "  ├── 4_dada2_outputs/"
echo "  └── 5_classifier_working_dir/"
echo "      ├── exported_rep_seqs/dna-sequences.fasta"
echo "      ├── ${GG2_FNA_FILE}"
echo "      ├── ${GG2_TAX_FILE}"
echo "      ├── length_summary.txt"
echo "      └── mbx_classifier_run_info.txt    ← your next commands are here"
echo ""
echo "  ── Next step ────────────────────────────────────────────────────"
echo "  Open mbx_classifier_run_info.txt and run Steps A → B → C → D → E"
echo "  cat $RUN_INFO"
echo ""
