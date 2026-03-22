#!/bin/bash
# setup.sh — MCCE4-Alpha Installer (Linux)
#
# Builds an Apptainer container with MCCE4, conda env, and all solvers.
# No sudo required — installs Apptainer to user-space if not available.
#
# Usage:
#   bash setup.sh                — smart re-run, skips completed steps
#   bash setup.sh --rebuild      — force full container rebuild
#   bash setup.sh --build-ngpb   — build NGPB from source instead of downloading
#
# macOS users: run setup_mac.sh instead.

trap "echo ''; echo '❌ Script interrupted by user'; exit 1" SIGINT

# ── Check OS ───────────────────────────────────────────────────────────────────
if [[ "$(uname -s)" != "Linux" ]]; then
    echo "❌ This script is for Linux only."
    echo "   Detected: $(uname -s)"
    echo ""
    echo "   macOS users: run setup_mac.sh instead."
    exit 1
fi

# ── Timing ─────────────────────────────────────────────────────────────────────
SETUP_START_TIME=$(date +%s)
TIMING_NAMES=()
TIMING_DURATIONS=()

_fmt_duration() {
    local secs=$1
    if [[ "$secs" -ge 3600 ]]; then
        printf "%dh %dm %ds" $((secs/3600)) $(((secs%3600)/60)) $((secs%60))
    elif [[ "$secs" -ge 60 ]]; then
        printf "%dm %ds" $((secs/60)) $((secs%60))
    else
        printf "%ds" "$secs"
    fi
}

_step_start_time=0
step_start() { _step_start_time=$(date +%s); echo ""; echo "⏱️  $1..."; }
step_end() {
    local d=$(( $(date +%s) - _step_start_time ))
    TIMING_NAMES+=("$1"); TIMING_DURATIONS+=("$d")
    echo "   ⏱️  $1: $(_fmt_duration $d)"
}

print_timing_summary() {
    local total=$(( $(date +%s) - SETUP_START_TIME ))
    echo ""; echo "⏱️  Timing Breakdown"
    echo "   ─────────────────────────────────────"
    local ml=0
    for n in "${TIMING_NAMES[@]}"; do [[ ${#n} -gt $ml ]] && ml=${#n}; done
    for i in "${!TIMING_NAMES[@]}"; do
        printf "   %-${ml}s  %s\n" "${TIMING_NAMES[$i]}" "$(_fmt_duration ${TIMING_DURATIONS[$i]})"
    done
    echo "   ─────────────────────────────────────"
    printf "   %-${ml}s  %s\n" "Total" "$(_fmt_duration $total)"
}

# ── Arguments ──────────────────────────────────────────────────────────────────
FORCE_REBUILD=0
BUILD_NGPB_SOURCE=0
for arg in "$@"; do
    [[ "$arg" == "--rebuild" ]]    && FORCE_REBUILD=1
    [[ "$arg" == "--build-ngpb" ]] && BUILD_NGPB_SOURCE=1
done

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT_DIR=$(dirname "$(readlink -f "$0")")
cd "$PROJECT_ROOT_DIR" || exit 1

BIN_DIR="$PROJECT_ROOT_DIR/bin"
CONDA_YML="$PROJECT_ROOT_DIR/mc4.yml"
CONDA_ENV_NAME="mc4"
MCCE4_ALPHA_DEF="$BIN_DIR/mcce4-alpha.def"
MCCE4_ALPHA_SIF="$BIN_DIR/mcce4-alpha.sif"

SHELL_CONFIG="$HOME/.bashrc"
[[ "${SHELL:-}" == *"zsh"* ]] && SHELL_CONFIG="$HOME/.zshrc"

# ── Shell Config Helper ───────────────────────────────────────────────────────
update_shell_config() {
    local dir_path="$1" tool_name="$2"
    local export_line="export PATH=\"$dir_path:\$PATH\""
    local marker="# $tool_name"

    if grep -qF "$marker" "$SHELL_CONFIG" 2>/dev/null; then
        sed -i "/$marker/,+1d" "$SHELL_CONFIG"
    fi

    if grep -qF "$dir_path" "$SHELL_CONFIG" 2>/dev/null; then
        echo "   ✅ PATH already contains $dir_path"
    else
        echo "" >> "$SHELL_CONFIG"
        echo "$marker" >> "$SHELL_CONFIG"
        echo "$export_line" >> "$SHELL_CONFIG"
        export PATH="$dir_path:$PATH"
        echo "   ✅ Added $dir_path to PATH"
    fi
}

# ── Conda Environment ─────────────────────────────────────────────────────────
setup_conda_env() {
    local force_apptainer="${1:-0}"

    echo "🐍 Setting up Conda environment '$CONDA_ENV_NAME'..."

    if ! command -v conda >/dev/null 2>&1; then
        echo "   ❌ 'conda' not found in PATH."
        echo "   👉 Install Miniconda (no sudo): https://docs.anaconda.com/miniconda/"
        echo "      Then: source ~/miniconda3/bin/activate && conda init"
        return 1
    fi

    if conda env list | grep -q "^${CONDA_ENV_NAME}\s\+"; then
        echo "   Updating existing environment..."
        yes "a" | conda env update -n "$CONDA_ENV_NAME" -f "$CONDA_YML"
    else
        echo "   Creating new environment..."
        yes "a" | conda env create -f "$CONDA_YML"
    fi

    if [[ "$force_apptainer" -eq 1 ]]; then
        echo "   ⬇️  Installing Apptainer via conda-forge..."
        conda install -n "$CONDA_ENV_NAME" -c conda-forge apptainer -y
        echo "   ✅ Apptainer installed in '$CONDA_ENV_NAME'."
    fi
}

# ── Apptainer Install ─────────────────────────────────────────────────────────
install_apptainer() {
    echo "🔍 Checking for Apptainer..."

    if command -v apptainer >/dev/null 2>&1; then
        echo "   ✅ Apptainer found: $(which apptainer)"
        return 0
    fi

    echo "   ⚠️  Apptainer not found."

    # Attempt 1: Unprivileged local install
    if command -v curl >/dev/null 2>&1; then
        echo "   ⬇️  Installing Apptainer locally (no sudo)..."
        local apt_dir="$HOME/apptainer"
        mkdir -p "$apt_dir"
        curl -fsSL \
            https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh \
            | bash -s - "$apt_dir"

        if "$apt_dir/bin/apptainer" version >/dev/null 2>&1; then
            export PATH="$apt_dir/bin:$PATH"
            update_shell_config "$apt_dir/bin" "Apptainer (Unprivileged)"
            echo "   ✅ Apptainer installed at $apt_dir"
            return 0
        else
            echo "   ⚠️  Unprivileged install failed."
            rm -rf "$apt_dir"
        fi
    fi

    # Attempt 2: conda-forge
    echo "   ⬇️  Falling back to conda-forge..."
    setup_conda_env 1
    eval "$(conda shell.bash hook)" 2>/dev/null
    conda activate "$CONDA_ENV_NAME" 2>/dev/null
    if command -v apptainer >/dev/null 2>&1; then
        echo "   ✅ Apptainer available via conda."
        return 0
    fi

    echo "   ❌ Could not install Apptainer via any method."
    return 1
}

# ── Container Build ────────────────────────────────────────────────────────────
should_rebuild_container() {
    [[ ! -f "$MCCE4_ALPHA_SIF" ]]              && echo "   No container found — building..." && return 0
    [[ "$FORCE_REBUILD" -eq 1 ]]               && echo "   --rebuild — forcing rebuild..." && return 0
    [[ "$CONDA_YML" -nt "$MCCE4_ALPHA_SIF" ]]  && echo "   ⚠️  mc4.yml newer — rebuilding..." && return 0
    [[ "$MCCE4_ALPHA_DEF" -nt "$MCCE4_ALPHA_SIF" ]] && echo "   ⚠️  .def newer — rebuilding..." && return 0
    echo "   ✅ Container up to date. Use --rebuild to force."
    return 1
}

build_container() {
    echo "🔧 Building Container Image: $MCCE4_ALPHA_SIF"

    mkdir -p "$(dirname "$MCCE4_ALPHA_SIF")"

    # Ensure apptainer is in PATH (may be in conda env)
    if ! command -v apptainer >/dev/null 2>&1; then
        if command -v conda >/dev/null 2>&1; then
            eval "$(conda shell.bash hook)" 2>/dev/null
            conda activate "$CONDA_ENV_NAME" 2>/dev/null
        fi
    fi

    if ! command -v apptainer >/dev/null 2>&1; then
        echo "   ❌ Apptainer not found."
        return 1
    fi

    # --force suppresses "overwrite?" prompt
    if apptainer build --force --fakeroot "$MCCE4_ALPHA_SIF" "$MCCE4_ALPHA_DEF"; then
        echo "   ✅ Build successful: $MCCE4_ALPHA_SIF"
        return 0
    fi

    echo "   ⚠️  --fakeroot failed. Retrying without it..."
    if apptainer build --force "$MCCE4_ALPHA_SIF" "$MCCE4_ALPHA_DEF"; then
        echo "   ✅ Build successful: $MCCE4_ALPHA_SIF"
        return 0
    fi

    echo "   ❌ Container build failed."
    return 1
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NGPB SOLVER (download pre-built or build from source)
# ═══════════════════════════════════════════════════════════════════════════════

NGPB_RELEASE="NextGenPB_v1.0.0"
NGPB_URL="https://github.com/concept-lab/NextGenPB/releases/download/$NGPB_RELEASE/NextGenPB.sif"
NGPB_GENERIC="$BIN_DIR/NextGenPB.sif"
NGPB_MC4="$BIN_DIR/NextGenPB_MCCE4.sif"
NGPB_RECIPE="$BIN_DIR/recipe_MCCE.def"

setup_ngpb() {
    echo "🧪 Setting up NextGenPB solver..."

    # Already present and not forcing
    if [[ -f "$NGPB_MC4" && "$FORCE_REBUILD" -eq 0 && "$BUILD_NGPB_SOURCE" -eq 0 ]]; then
        echo "   ✅ NextGenPB_MCCE4.sif exists. Use --rebuild or --build-ngpb to replace."
        return 0
    fi

    if [[ "$BUILD_NGPB_SOURCE" -eq 1 ]]; then
        # ── Build from source ────────────────────────────────────────────────
        if [[ ! -f "$NGPB_RECIPE" ]]; then
            echo "   ❌ bin/recipe_MCCE.def not found. Cannot build from source."
            echo "      Use default download instead (omit --build-ngpb)."
            return 1
        fi

        if ! command -v apptainer >/dev/null 2>&1; then
            echo "   ❌ apptainer required to build NGPB from source."
            return 1
        fi

        echo "   🔨 Building NextGenPB from source (~20-30 min)..."
        if apptainer build --force --fakeroot "$NGPB_MC4" "$NGPB_RECIPE"; then
            echo "   ✅ NextGenPB_MCCE4.sif built from source."
            return 0
        fi

        echo "   ❌ NGPB source build failed."
        return 1
    else
        # ── Download pre-built (~1.6GB) ──────────────────────────────────────
        if [[ ! -f "$NGPB_GENERIC" || "$FORCE_REBUILD" -eq 1 ]]; then
            echo "   ⬇️  Downloading NextGenPB container (~1.6GB)..."
            curl -L -o "$NGPB_GENERIC" "$NGPB_URL" || {
                echo "   ⚠️  Download failed. DELPHI/APBS will still work."
                echo "      Retry later or build from source: bash setup.sh --build-ngpb"
                return 0
            }
        fi

        if [[ -f "$NGPB_GENERIC" ]]; then
            chmod +x "$NGPB_GENERIC"
            ln -sf "$(basename "$NGPB_GENERIC")" "$NGPB_MC4"
            echo "   ✅ NextGenPB installed: NextGenPB_MCCE4.sif → NextGenPB.sif"
        fi
        return 0
    fi
}


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🔧 MCCE4-Alpha Setup (Linux)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "   Clone root : $PROJECT_ROOT_DIR"
echo "   Platform   : Linux / $(uname -m)"
[[ "$FORCE_REBUILD" -eq 1 ]]    && echo "   Mode       : full rebuild (--rebuild)"
[[ "$BUILD_NGPB_SOURCE" -eq 1 ]] && echo "   NGPB       : build from source (--build-ngpb)"
echo ""

step_start "Apptainer install"
if ! install_apptainer; then
    echo ""
    echo "⛔ Could not install Apptainer. Cannot proceed."
    echo "   Ask your sysadmin to install it, or install Miniconda and re-run."
    exit 1
fi
step_end "Apptainer install"

update_shell_config "$BIN_DIR" "MCCE4-Alpha CLI"

# Temporarily move any existing .sif files out of bin/ so they don't get
# baked into the container (they'd make mksquashfs OOM on low-memory servers).
# They'll be bind-mounted at runtime by mc4 instead.
_sif_tmp=""
if ls "$BIN_DIR"/*.sif >/dev/null 2>&1; then
    _sif_tmp=$(mktemp -d)
    echo "   Moving .sif files out of bin/ temporarily (will be restored)..."
    mv "$BIN_DIR"/*.sif "$_sif_tmp/" 2>/dev/null || true
fi

step_start "Container build"
if should_rebuild_container; then
    if ! build_container; then
        # Restore .sif files before exiting
        [[ -n "$_sif_tmp" ]] && mv "$_sif_tmp"/*.sif "$BIN_DIR/" 2>/dev/null && rm -rf "$_sif_tmp"
        echo ""
        echo "⛔ Container build failed. Cannot proceed."
        exit 1
    fi
fi
step_end "Container build"

# Restore .sif files
[[ -n "$_sif_tmp" ]] && mv "$_sif_tmp"/*.sif "$BIN_DIR/" 2>/dev/null && rm -rf "$_sif_tmp"

# Download/build NGPB AFTER container build (kept on host, bind-mounted at runtime)
step_start "NGPB solver"
setup_ngpb
step_end "NGPB solver"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🚀 MCCE4-Alpha Setup Complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "   Restart your shell or run: source $SHELL_CONFIG"
echo ""
echo "   Mode    : Apptainer container"
echo "   Image   : $MCCE4_ALPHA_SIF"
echo "   Solvers : DELPHI ✅  APBS ✅  NGPB ✅"
echo ""
echo "   Usage:"
echo "     mc4 step1.py <pdbfile>"
echo "     mc4 step2.py"
echo "     mc4 step3.py"
echo "     mc4 step4.py"
echo "     mc4 --shell"

print_timing_summary

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
