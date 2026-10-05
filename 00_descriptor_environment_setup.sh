#!/usr/bin/env bash
# =============================================================================
# 00_descriptor_environment_setup.sh
# =============================================================================
# Creates and validates a dedicated conda environment for Mordred + PaDEL
# descriptor benchmarking on Apple Silicon (M1/M2/M3 – osx-arm64).
#
# Usage:
#   bash 00_descriptor_environment_setup.sh
#
# What this script does:
#   1. Creates  conda env  "descriptor_env"  (Python 3.11, osx-arm64)
#   2. Installs rdkit, pandas, numpy, tqdm, requests, pathlib via conda-forge
#   3. Installs optional scientific libs: scipy, scikit-learn, joblib
#   4. Installs Mordred  via pip
#   5. Installs padelpy  via pip
#   6. Runs import + version checks and writes descriptor_environment_report.txt
#
# Requirements:
#   - Miniforge / Mambaforge (or Anaconda) installed for osx-arm64
#   - Internet connection
#   - No existing environment named "descriptor_env" is modified; the script
#     refuses to overwrite and asks the user to remove it first.
# =============================================================================

set -euo pipefail          # exit on error, undefined var, or pipe failure

# ─── colour helpers ──────────────────────────────────────────────────────────
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error()   { echo -e "${RED}[ERROR]${RESET} $*" >&2; }
section() { echo -e "\n${BOLD}${CYAN}══════════════════════════════════════════════${RESET}"; \
            echo -e "${BOLD}${CYAN}  $*${RESET}"; \
            echo -e "${BOLD}${CYAN}══════════════════════════════════════════════${RESET}"; }

# ─── configuration ────────────────────────────────────────────────────────────
ENV_NAME="descriptor_env"
PYTHON_VER="3.11"
REPORT_FILE="descriptor_environment_report.txt"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPORT_PATH="${SCRIPT_DIR}/${REPORT_FILE}"

# ─── locate conda ─────────────────────────────────────────────────────────────
section "STEP 0 – Locating conda"

# Source conda's shell integration so 'conda activate' works inside a script
CONDA_SH=""
for candidate in \
        "${HOME}/miniforge3/etc/profile.d/conda.sh" \
        "${HOME}/mambaforge/etc/profile.d/conda.sh" \
        "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh" \
        "/opt/anaconda3/etc/profile.d/conda.sh" \
        "/opt/miniconda3/etc/profile.d/conda.sh" \
        "${CONDA_PREFIX:-}/etc/profile.d/conda.sh"; do
    if [[ -f "$candidate" ]]; then
        CONDA_SH="$candidate"
        break
    fi
done

if [[ -z "$CONDA_SH" ]]; then
    error "conda shell initialisation script not found."
    error "Please ensure Miniforge/Anaconda is installed for osx-arm64."
    exit 1
fi

# shellcheck source=/dev/null
source "$CONDA_SH"
info "conda shell integration sourced from: ${CONDA_SH}"

# Confirm we are on Apple Silicon
ARCH="$(uname -m)"
info "Detected architecture: ${ARCH}"
if [[ "$ARCH" != "arm64" ]]; then
    warn "Expected arm64 (Apple Silicon) but found ${ARCH}."
    warn "The script will continue, but packages may not be ARM-native."
fi

# Conda version
CONDA_VER="$(conda --version 2>&1)"
info "Using ${CONDA_VER}"

# ─── guard: refuse to overwrite existing env ───────────────────────────────────
section "STEP 1 – Checking for existing environment"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    error "Environment '${ENV_NAME}' already exists."
    error "To rebuild it, first remove it:"
    error "    conda env remove -n ${ENV_NAME} -y"
    error "Then re-run this script."
    exit 1
fi
success "No existing '${ENV_NAME}' environment found – proceeding."

# ─── create env ───────────────────────────────────────────────────────────────
section "STEP 2 – Creating conda environment (Python ${PYTHON_VER})"

info "This may take a few minutes while conda resolves dependencies …"
conda create -n "${ENV_NAME}" \
    python="${PYTHON_VER}" \
    --channel conda-forge \
    --override-channels \
    --yes \
    --quiet

success "Environment '${ENV_NAME}' created."

# Activate the new env for all subsequent commands
conda activate "${ENV_NAME}"
info "Active environment: ${CONDA_DEFAULT_ENV}"

# ─── conda-forge packages ─────────────────────────────────────────────────────
section "STEP 3 – Installing core packages via conda-forge"

CONDA_PKGS=(
    rdkit
    pandas
    numpy
    tqdm
    requests
    pathlib2        # backport shim; pathlib is stdlib in 3.11 but harmless
    scipy
    scikit-learn
    joblib
)

info "Installing: ${CONDA_PKGS[*]}"
conda install -n "${ENV_NAME}" \
    "${CONDA_PKGS[@]}" \
    --channel conda-forge \
    --override-channels \
    --yes \
    --quiet \
    && success "conda-forge packages installed." \
    || { warn "One or more conda-forge packages failed."; \
         warn "Troubleshooting:"; \
         warn "  1. Check your internet connection."; \
         warn "  2. Try: conda clean --all -y  then re-run."; \
         warn "  3. Run: conda install -n ${ENV_NAME} rdkit -c conda-forge --debug"; }

# ─── pip packages ─────────────────────────────────────────────────────────────
section "STEP 4 – Installing Mordred via pip"

# Point to the pip inside the new env (conda activate may not export PATH yet)
ENV_PIP="$(conda run -n "${ENV_NAME}" which pip)"
info "pip binary: ${ENV_PIP}"

install_pip_pkg() {
    local pkg="$1"
    info "pip install ${pkg} …"
    if conda run -n "${ENV_NAME}" pip install --quiet "${pkg}"; then
        success "${pkg} installed via pip."
    else
        warn "pip install of '${pkg}' failed."
        warn "Troubleshooting suggestions:"
        warn "  • Ensure pip is up-to-date: conda run -n ${ENV_NAME} pip install --upgrade pip"
        warn "  • Check network / proxy settings."
        case "${pkg}" in
            mordred)
                warn "  • Mordred alternative: pip install mordred --no-deps"
                warn "  • GitHub source: pip install git+https://github.com/mordred-descriptor/mordred.git"
                warn "  • Ensure rdkit is installed first (conda-forge step above)."
                ;;
            padelpy)
                warn "  • padelpy needs Java ≥ 8 at runtime. Install via:"
                warn "    brew install openjdk  (then follow its PATH instructions)"
                warn "  • Verify: java -version"
                warn "  • padelpy GitHub: https://github.com/ecrl/padelpy"
                ;;
        esac
    fi
}

install_pip_pkg "mordred"

section "STEP 5 – Installing padelpy via pip"
install_pip_pkg "padelpy"

# ─── Java check (needed by padelpy / PaDEL-Descriptor) ────────────────────────
section "STEP 5a – Verifying Java runtime (required by padelpy)"

if command -v java &>/dev/null; then
    JAVA_VER="$(java -version 2>&1 | head -1)"
    success "Java found: ${JAVA_VER}"
else
    warn "Java not found on PATH."
    warn "padelpy uses PaDEL-Descriptor (a Java application) at runtime."
    warn "Install Java via Homebrew:"
    warn "    brew install openjdk"
    warn "    echo 'export PATH=\"/opt/homebrew/opt/openjdk/bin:\$PATH\"' >> ~/.zshrc"
    warn "    source ~/.zshrc"
    warn "Descriptor calculations with padelpy will fail without Java."
fi

# ─── verification ─────────────────────────────────────────────────────────────
section "STEP 6 – Running import verification"

VERIFY_SCRIPT=$(cat <<'PYEOF'
import sys, importlib, traceback, datetime, platform

sep = "=" * 60
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

lines = []
lines.append(sep)
lines.append("  DESCRIPTOR ENVIRONMENT VERIFICATION REPORT")
lines.append(f"  Generated : {now}")
lines.append(f"  Python    : {sys.version}")
lines.append(f"  Platform  : {platform.platform()}")
lines.append(f"  Machine   : {platform.machine()}")
lines.append(sep)
lines.append("")

PACKAGES = {
    "rdkit"       : ("rdkit",        "rdBase",         "__version__"),
    "mordred"     : ("mordred",      None,             "__version__"),
    "padelpy"     : ("padelpy",      None,             "__version__"),
    "pandas"      : ("pandas",       None,             "__version__"),
    "numpy"       : ("numpy",        None,             "__version__"),
    "scipy"       : ("scipy",        None,             "__version__"),
    "sklearn"     : ("sklearn",      None,             "__version__"),
    "joblib"      : ("joblib",       None,             "__version__"),
    "tqdm"        : ("tqdm",         None,             "__version__"),
    "requests"    : ("requests",     None,             "__version__"),
}

success_count = 0
fail_count    = 0
results       = []

for label, (mod, submod, ver_attr) in PACKAGES.items():
    try:
        m = importlib.import_module(mod)
        if submod:
            m = importlib.import_module(f"{mod}.{submod}")
        ver = getattr(m, ver_attr, "unknown")
        results.append(("OK", label, str(ver), ""))
        success_count += 1
    except Exception as e:
        tb = traceback.format_exc().strip().splitlines()[-1]
        results.append(("FAIL", label, "—", tb))
        fail_count += 1

lines.append(f"{'Status':<8} {'Package':<12} {'Version':<20} Notes")
lines.append(f"{'-'*8} {'-'*12} {'-'*20} {'-'*30}")
for status, pkg, ver, note in results:
    flag = "✓" if status == "OK" else "✗"
    lines.append(f"  {flag} {status:<6} {pkg:<12} {ver:<20} {note}")

lines.append("")
lines.append(f"Summary: {success_count} passed, {fail_count} failed")
lines.append("")

# ── Quick functional test for rdkit ──────────────────────────────────────────
lines.append(sep)
lines.append("  RDKIT FUNCTIONAL TEST")
lines.append(sep)
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")   # aspirin
    mw  = Descriptors.MolWt(mol)
    lines.append(f"  Aspirin SMILES parsed OK | MolWt = {mw:.4f} g/mol")
    lines.append("  RDKit functional test: PASSED")
except Exception as e:
    lines.append(f"  RDKit functional test: FAILED – {e}")

# ── Quick functional test for Mordred ────────────────────────────────────────
lines.append("")
lines.append(sep)
lines.append("  MORDRED FUNCTIONAL TEST")
lines.append(sep)
try:
    from rdkit import Chem
    from mordred import Calculator, descriptors
    mol  = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    calc = Calculator(descriptors, ignore_3D=True)
    result = calc(mol)
    n_desc = len(result)
    lines.append(f"  Aspirin: {n_desc} Mordred descriptors computed")
    lines.append("  Mordred functional test: PASSED")
except Exception as e:
    lines.append(f"  Mordred functional test: FAILED – {e}")

# ── padelpy import test (Java runtime check) ─────────────────────────────────
lines.append("")
lines.append(sep)
lines.append("  PADELPY IMPORT TEST")
lines.append(sep)
try:
    import padelpy
    lines.append(f"  padelpy imported OK (version: {getattr(padelpy,'__version__','unknown')})")
    lines.append("  NOTE: Full padelpy functionality requires Java ≥ 8 at runtime.")
    lines.append("  To verify: run  java -version  in your shell.")
except Exception as e:
    lines.append(f"  padelpy import FAILED – {e}")
    lines.append("  Ensure padelpy is installed: pip install padelpy")

lines.append("")
lines.append(sep)
lines.append("  END OF REPORT")
lines.append(sep)

report_text = "\n".join(lines)
print(report_text)

# Write report to file (passed via env var)
import os
report_path = os.environ.get("REPORT_PATH", "descriptor_environment_report.txt")
with open(report_path, "w", encoding="utf-8") as fh:
    fh.write(report_text + "\n")

print(f"\n[REPORT SAVED] → {report_path}")
sys.exit(0 if fail_count == 0 else 1)
PYEOF
)

# Run the verification inside the new env
export REPORT_PATH="${REPORT_PATH}"
if conda run -n "${ENV_NAME}" python -c "${VERIFY_SCRIPT}"; then
    success "All packages verified successfully."
else
    warn "Some packages failed verification. See ${REPORT_FILE} for details."
    warn ""
    warn "Common fixes:"
    warn "  • Mordred needs rdkit ≥ 2022. Confirm: conda run -n ${ENV_NAME} python -c 'from rdkit import rdBase; print(rdBase.rdkitVersion)'"
    warn "  • padelpy needs Java:  brew install openjdk"
    warn "  • Re-install a broken package:"
    warn "      conda run -n ${ENV_NAME} pip install --force-reinstall <package>"
    warn "  • For SSL/certificate errors: pip install --trusted-host pypi.org <package>"
fi

# ─── activation hint ──────────────────────────────────────────────────────────
section "SETUP COMPLETE"

echo ""
echo -e "${BOLD}To activate the environment in a new shell:${RESET}"
echo -e "    ${GREEN}conda activate ${ENV_NAME}${RESET}"
echo ""
echo -e "${BOLD}To run a script inside it without activating:${RESET}"
echo -e "    ${GREEN}conda run -n ${ENV_NAME} python your_script.py${RESET}"
echo ""
echo -e "${BOLD}Report saved to:${RESET}"
echo -e "    ${GREEN}${REPORT_PATH}${RESET}"
echo ""
echo -e "${BOLD}To remove this environment later:${RESET}"
echo -e "    ${YELLOW}conda env remove -n ${ENV_NAME} -y${RESET}"
echo ""
