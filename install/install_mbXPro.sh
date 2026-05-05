#!/usr/bin/env bash
# =============================================================================
#  install_mbXPro.sh -- One-command installer for the mbX Pro pipeline
# =============================================================================
#
#  WHAT IT DOES
#    1. Verifies system requirements (macOS or Linux, Bash >= 3.2, conda, R)
#    2. Copies all 21 step scripts + the orchestrator to ~/bin/
#    3. Copies the mbX Pro logo (used by the final-report step) to ~/bin/
#    4. Adds ~/bin to your PATH in ~/.zshrc and / or ~/.bashrc (if missing)
#    5. Prints a one-line "ready to use" message
#
#  USAGE
#    bash install_mbXPro.sh
#    bash install_mbXPro.sh --prefix /custom/install/dir   (advanced)
#    bash install_mbXPro.sh --no-path                       (skip PATH edit)
#    bash install_mbXPro.sh -h | --help
#
#  WHAT IT DOES NOT DO
#    * It does NOT install QIIME2 (you must do that yourself; see
#      documentation/how_to_run_mbX_Pro.docx).
#    * It does NOT install R (you should `brew install r` on macOS or
#      `sudo apt install r-base` on Linux).
#    * It does NOT install PICRUSt2 (the picrust step auto-installs into
#      its own conda env on first run).
#
# =============================================================================

set -euo pipefail

# ── Pretty helpers ────────────────────────────────────────────────────────────
sep()  { printf '\n────────────────────────────────────────────────────────────────\n'; }
info() { printf '[INFO]  %s\n' "$*"; }
ok()   { printf '[OK]    %s\n' "$*"; }
warn() { printf '[WARN]  %s\n' "$*" >&2; }
err()  { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

# ── Args ──────────────────────────────────────────────────────────────────────
PREFIX="$HOME/bin"
EDIT_PATH=true

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)   PREFIX="$2"; shift 2 ;;
    --no-path)  EDIT_PATH=false; shift ;;
    -h|--help)
      awk '/^# ===/{n++;next} n==1 && /^#/{sub(/^# ?/,"");print} n==2{exit}' "$0"
      exit 0 ;;
    *) err "Unknown option: $1" ;;
  esac
done

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$THIS_DIR/.." && pwd)"
SCRIPTS_DIR="$ROOT_DIR/scripts"
ASSETS_DIR="$ROOT_DIR/assets"

[[ -d "$SCRIPTS_DIR" ]] || err "Cannot find scripts directory at: $SCRIPTS_DIR
  -> Run this installer from inside the cloned mbXPro/ directory:
        cd mbXPro/install && bash install_mbXPro.sh"

VERSION_FILE="$ROOT_DIR/VERSION"
[[ -f "$VERSION_FILE" ]] && MBX_VERSION="$(head -n 1 "$VERSION_FILE" | tr -d '[:space:]')" || MBX_VERSION="(unknown)"

cat <<BANNER

  ╔══════════════════════════════════════════════════════════════════╗
  ║              mbX Pro Installer  --  setup wizard                 ║
  ║                                                                  ║
  ║  Version : $(printf '%-54s' "$MBX_VERSION")║
  ║  Notes   : self-healing pipeline.  Step 0 decides DETECTION_     ║
  ║            STATUS and step 3 sets DADA2 trim-left accordingly    ║
  ║            (FOUND -> primer length, TRIMMED -> 0, UNKNOWN -> 20).║
  ║            Step 5 prefers a pre-trained classifier from Zenodo   ║
  ║            (https://zenodo.org/records/20021035) and falls back  ║
  ║            to local training automatically if anything fails.    ║
  ╚══════════════════════════════════════════════════════════════════╝

BANNER
info "Source repo : $ROOT_DIR"
info "Install to  : $PREFIX"
sep

# ── Pre-flight ────────────────────────────────────────────────────────────────
info "Checking your system..."

# OS
case "$(uname -s)" in
  Darwin) OS_NAME="macOS" ;;
  Linux)  OS_NAME="Linux" ;;
  *)      err "Unsupported OS: $(uname -s).  mbX Pro supports macOS and Linux only." ;;
esac
ok "OS              : $OS_NAME ($(uname -m))"

# Bash version
BV="${BASH_VERSION:-unknown}"
ok "Bash            : $BV"

# Conda + QIIME2 (warn-only; user can install later)
if command -v conda &>/dev/null; then
  CONDA_VER="$(conda --version 2>/dev/null | awk '{print $2}')"
  ok "conda           : $CONDA_VER  (found)"
  if conda env list 2>/dev/null | grep -q '^qiime2-amplicon-2025'; then
    QENV="$(conda env list | awk '/qiime2-amplicon-2025/{print $1; exit}')"
    ok "QIIME2 env      : $QENV  (found -- ready)"
  else
    warn "QIIME2 env not found.  Install with:"
    warn "  curl -LO https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2025.4-py310-osx-conda.yml"
    warn "  conda env create -n qiime2-amplicon-2025.4 --file qiime2-amplicon-2025.4-py310-osx-conda.yml"
    warn "  (Replace 'osx' with 'linux' on Linux.)"
  fi
else
  warn "conda not found in PATH."
  warn "  -> Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
  warn "  -> Then install the QIIME2 amplicon environment (see how_to_run_mbX_Pro.docx)."
fi

# Rscript (system-wide)
RSCRIPT_BIN=""
for _r in /usr/local/bin/Rscript /opt/homebrew/bin/Rscript /usr/bin/Rscript; do
  if [[ -x "$_r" ]]; then RSCRIPT_BIN="$_r"; break; fi
done
[[ -z "$RSCRIPT_BIN" ]] && command -v Rscript &>/dev/null && RSCRIPT_BIN="$(command -v Rscript)"

if [[ -n "$RSCRIPT_BIN" ]]; then
  R_VER="$("$RSCRIPT_BIN" --version 2>&1 | head -1 | awk '{print $4}')"
  ok "R (system)      : $RSCRIPT_BIN  ($R_VER)"
else
  warn "System-wide Rscript NOT found."
  if [[ "$OS_NAME" == "macOS" ]]; then
    warn "  -> Install with:  brew install r"
  else
    warn "  -> Install with:  sudo apt-get install r-base    (Debian/Ubuntu)"
    warn "                    sudo dnf  install R           (Fedora/RHEL)"
  fi
  warn "  -> R must NOT be installed inside the QIIME2 conda env (would clash)."
fi

sep

# ── Install scripts ───────────────────────────────────────────────────────────
mkdir -p "$PREFIX" || err "Cannot create install dir: $PREFIX"
info "Copying step scripts to $PREFIX ..."

N_INSTALLED=0
for f in "$SCRIPTS_DIR"/*.sh "$SCRIPTS_DIR"/mbXPro; do
  [[ -e "$f" ]] || continue
  cp -f "$f" "$PREFIX/"
  chmod +x "$PREFIX/$(basename "$f")"
  N_INSTALLED=$((N_INSTALLED + 1))
done
ok "Installed $N_INSTALLED scripts into $PREFIX/"

# Logo (used by mbx_final_report.sh for the report header)
if [[ -f "$ASSETS_DIR/mbX_Pro_icon.png" ]]; then
  cp -f "$ASSETS_DIR/mbX_Pro_icon.png" "$PREFIX/"
  ok "Installed mbX_Pro_icon.png   (used by the final-report step)"
fi

# ── PATH edit ─────────────────────────────────────────────────────────────────
if $EDIT_PATH; then
  PATH_LINE="export PATH=\"$PREFIX:\$PATH\"   # added by install_mbXPro.sh"
  RC_FILES=()
  [[ -f "$HOME/.zshrc"  ]] && RC_FILES+=( "$HOME/.zshrc"  )
  [[ -f "$HOME/.bashrc" ]] && RC_FILES+=( "$HOME/.bashrc" )
  [[ -f "$HOME/.bash_profile" ]] && RC_FILES+=( "$HOME/.bash_profile" )

  if [[ ${#RC_FILES[@]} -eq 0 ]]; then
    warn "No shell rc file found (~/.zshrc, ~/.bashrc, ~/.bash_profile)."
    warn "  -> Manually add this line to your shell startup file:"
    warn "       $PATH_LINE"
  else
    for rc in "${RC_FILES[@]}"; do
      if ! grep -q "install_mbXPro.sh" "$rc" 2>/dev/null; then
        printf '\n# mbX Pro pipeline\n%s\n' "$PATH_LINE" >> "$rc"
        ok "Added PATH entry to $rc"
      else
        info "PATH entry already present in $rc -- left untouched"
      fi
    done
  fi
else
  warn "--no-path: skipping PATH edit.  You must add this manually:"
  warn "  export PATH=\"$PREFIX:\$PATH\""
fi

sep

# ── Final summary ────────────────────────────────────────────────────────────
cat <<EOM

  ╔══════════════════════════════════════════════════════════════════╗
  ║                                                                  ║
  ║           mbX Pro Installation -- COMPLETE                       ║
  ║                                                                  ║
  ╚══════════════════════════════════════════════════════════════════╝

  Installed to    : $PREFIX
  Scripts copied  : $N_INSTALLED
  Logo copied     : yes
  PATH updated    : $($EDIT_PATH && echo yes || echo "no (--no-path)")

  Next steps:
    1) Open a NEW terminal (or run: source ~/.zshrc).
    2) Verify install:    mbXPro --help
    3) Activate QIIME2:   conda activate qiime2-amplicon-2025.4
    4) Run the pipeline:  mbXPro <fastq_dir> <metadata.txt>

  See documentation/how_to_run_mbX_Pro.docx for a full walk-through.

EOM
