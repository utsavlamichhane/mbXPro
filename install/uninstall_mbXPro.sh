#!/usr/bin/env bash
# =============================================================================
#  uninstall_mbXPro.sh -- Remove every file installed by install_mbXPro.sh
# =============================================================================
#
#  USAGE
#    bash uninstall_mbXPro.sh
#    bash uninstall_mbXPro.sh --prefix /custom/install/dir
#    bash uninstall_mbXPro.sh -h | --help
#
#  WHAT IT DOES
#    1. Deletes every script (and the logo) installed under <prefix>/
#       that came from this mbX Pro repo.
#    2. Removes the PATH-edit line from ~/.zshrc / ~/.bashrc / ~/.bash_profile.
#
#  WHAT IT DOES NOT DO
#    * Does not delete any of YOUR analysis output directories
#      (mbX_pro_outputs_*).  Your data is safe.
#    * Does not uninstall conda, QIIME2, R, or any third-party tool.
#
# =============================================================================

set -euo pipefail
sep()  { printf '\n────────────────────────────────────────────────────────────────\n'; }
info() { printf '[INFO]  %s\n' "$*"; }
ok()   { printf '[OK]    %s\n' "$*"; }
warn() { printf '[WARN]  %s\n' "$*" >&2; }
err()  { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

PREFIX="$HOME/bin"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)  PREFIX="$2"; shift 2 ;;
    -h|--help)
      awk '/^# ===/{n++;next} n==1 && /^#/{sub(/^# ?/,"");print} n==2{exit}' "$0"
      exit 0 ;;
    *) err "Unknown option: $1" ;;
  esac
done

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$THIS_DIR/.." && pwd)"
SCRIPTS_DIR="$ROOT_DIR/scripts"

[[ -d "$SCRIPTS_DIR" ]] || err "Cannot find scripts directory at $SCRIPTS_DIR"

cat <<BANNER

  ╔══════════════════════════════════════════════════════════════════╗
  ║          mbX Pro Uninstaller  --  removing installed files        ║
  ╚══════════════════════════════════════════════════════════════════╝

BANNER
info "Install prefix : $PREFIX"
sep

# Remove scripts + logo
N_REMOVED=0
for f in "$SCRIPTS_DIR"/*.sh "$SCRIPTS_DIR"/mbXPro; do
  [[ -e "$f" ]] || continue
  bn="$(basename "$f")"
  if [[ -f "$PREFIX/$bn" ]]; then
    rm -f "$PREFIX/$bn"
    N_REMOVED=$((N_REMOVED + 1))
  fi
done
[[ -f "$PREFIX/mbX_Pro_icon.png" ]] && rm -f "$PREFIX/mbX_Pro_icon.png" && N_REMOVED=$((N_REMOVED + 1))
ok "Removed $N_REMOVED files from $PREFIX/"

# Strip PATH-edit lines from rc files
for rc in "$HOME/.zshrc" "$HOME/.bashrc" "$HOME/.bash_profile"; do
  [[ -f "$rc" ]] || continue
  if grep -q "install_mbXPro.sh" "$rc"; then
    # Use a temp file (in-place sed differs across BSD/GNU)
    tmp="$(mktemp)"
    grep -v -E "(install_mbXPro\.sh|mbX Pro pipeline)" "$rc" > "$tmp" || true
    mv "$tmp" "$rc"
    ok "Stripped PATH-edit from $rc"
  fi
done

sep
cat <<EOM

  ╔══════════════════════════════════════════════════════════════════╗
  ║         mbX Pro Uninstallation -- COMPLETE                       ║
  ╚══════════════════════════════════════════════════════════════════╝

  Your analysis output directories (mbX_pro_outputs_*) were left untouched.
  You may also want to: conda env remove -n qiime2-amplicon-2025.4
                        (only if you no longer need QIIME2 at all)

EOM
