#!/bin/bash
# setup.sh — MCCE4-Alpha Installer
#
# Unified installer for Linux and macOS:
#
#   Linux  → Apptainer container (mcce4-alpha.sif)
#            All solvers including NGPB work directly.
#            Fallback: native conda + compiled binaries.
#
#   macOS  → Lima Linux VM with native build inside the VM.
#            Your home directory is mounted writable — code edits on Mac are
#            instantly visible. MCCE4 compiles from source inside the VM (~30s).
#            delphi_precompiled works (Linux ELF in a Linux VM).
#            Apptainer is available for NGPB.
#            No 20-30 min container rebuilds during development.
#
# No sudo required. All tools install to user-space.
#
# Usage:
#   bash setup.sh            — smart re-run: skips completed steps
#   bash setup.sh --rebuild  — force full rebuild

set -euo pipefail
trap "echo ''; echo '❌ Script interrupted by user'; exit 1" SIGINT

# ── Timing utilities ──────────────────────────────────────────────────────────
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
step_start() {
    _step_start_time=$(date +%s)
    echo ""
    echo "⏱️  $1..."
}

step_end() {
    local end_time=$(date +%s)
    local duration=$(( end_time - _step_start_time ))
    TIMING_NAMES+=("$1")
    TIMING_DURATIONS+=("$duration")
    echo "   ⏱️  $1: $(_fmt_duration $duration)"
}

print_timing_summary() {
    local total_time=$(( $(date +%s) - SETUP_START_TIME ))
    echo ""
    echo "⏱️  Timing Breakdown"
    echo "   ─────────────────────────────────────"
    local max_len=0
    for name in "${TIMING_NAMES[@]}"; do
        [[ ${#name} -gt $max_len ]] && max_len=${#name}
    done
    for i in "${!TIMING_NAMES[@]}"; do
        printf "   %-${max_len}s  %s\n" "${TIMING_NAMES[$i]}" "$(_fmt_duration ${TIMING_DURATIONS[$i]})"
    done
    echo "   ─────────────────────────────────────"
    printf "   %-${max_len}s  %s\n" "Total" "$(_fmt_duration $total_time)"
}

# ── Parse arguments ────────────────────────────────────────────────────────────
FORCE_REBUILD=0
for arg in "$@"; do
    [[ "$arg" == "--rebuild" ]] && FORCE_REBUILD=1
done

# ── Portable realpath ──────────────────────────────────────────────────────────
_realpath() {
    if command -v python3 >/dev/null 2>&1; then
        python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$1"
    elif command -v realpath >/dev/null 2>&1; then
        realpath "$1"
    else
        cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"
    fi
}

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT_DIR=$(dirname "$(_realpath "$0")")
cd "$PROJECT_ROOT_DIR" || exit 1

BIN_DIR="$PROJECT_ROOT_DIR/bin"
CONDA_YML="$PROJECT_ROOT_DIR/mc4.yml"
CONDA_ENV_NAME="mc4"
MCCE4_ALPHA_DEF="$BIN_DIR/mcce4-alpha.def"

INSTALL_DIR="$(dirname "$PROJECT_ROOT_DIR")"
MCCE4_ALPHA_SIF="$INSTALL_DIR/mcce4-alpha.sif"
APPTAINER_LOCAL_DIR="$INSTALL_DIR/apptainer"

# Lima VM settings (macOS)
LIMA_INSTANCE="mcce4"
LIMA_CONFIG="$BIN_DIR/mcce4-lima.yaml"

SHELL_CONFIG="$HOME/.bashrc"
[[ "${SHELL:-}" == *"zsh"* ]] && SHELL_CONFIG="$HOME/.zshrc"

# ── Platform Detection ─────────────────────────────────────────────────────────
OS=""
ARCH=""

detect_platform() {
    case "$(uname -s)" in
        Linux*)  OS="linux"  ;;
        Darwin*) OS="macos"  ;;
        *)       echo "❌ Unsupported OS: $(uname -s)"; exit 1 ;;
    esac
    case "$(uname -m)" in
        x86_64)        ARCH="amd64"  ;;
        aarch64|arm64) ARCH="arm64"  ;;
        *)             echo "❌ Unsupported architecture: $(uname -m)"; exit 1 ;;
    esac
    echo "🖥️  Detected platform: $OS / $ARCH"
}

# ── Shell Config Helper ────────────────────────────────────────────────────────
update_shell_config() {
    local dir_path="$1"
    local tool_name="$2"
    local export_line="export PATH=\"$dir_path:\$PATH\""
    local marker="# $tool_name"
    if grep -qF "$marker" "$SHELL_CONFIG" 2>/dev/null; then
        sed -i.bak "/$marker/,+1d" "$SHELL_CONFIG"
        rm -f "${SHELL_CONFIG}.bak"
    fi
    if grep -qF "$dir_path" "$SHELL_CONFIG" 2>/dev/null; then
        echo "   ✅ PATH already contains $dir_path"
    else
        printf "\n%s\n%s\n" "$marker" "$export_line" >> "$SHELL_CONFIG"
        export PATH="$dir_path:$PATH"
        echo "   ✅ Added $dir_path to PATH in $SHELL_CONFIG"
    fi
}


# ═══════════════════════════════════════════════════════════════════════════════
#  macOS: LIMA VM WITH NATIVE BUILD (no container needed)
# ═══════════════════════════════════════════════════════════════════════════════

generate_lima_config() {
    echo "   Generating Lima VM configuration..."
    cat > "$LIMA_CONFIG" << 'YAML'
# mcce4-lima.yaml — Lima VM for MCCE4-Alpha
# Provides a full Linux dev environment: gcc, conda, apptainer.
# Home directory is mounted writable — Mac edits are instantly visible.

images:
  - location: "https://cloud-images.ubuntu.com/releases/24.04/release/ubuntu-24.04-server-cloudimg-amd64.img"
    arch: "x86_64"
  - location: "https://cloud-images.ubuntu.com/releases/24.04/release/ubuntu-24.04-server-cloudimg-arm64.img"
    arch: "aarch64"

mounts:
  - location: "~"
    writable: true
  - location: "/tmp/lima"
    writable: true

provision:
  - mode: system
    script: |
      #!/bin/bash
      set -eux -o pipefail

      # Skip full provision if already done
      [ -f /etc/mcce4-provisioned ] && exit 0

      apt-get update
      apt-get install -y --no-install-recommends \
          build-essential gcc g++ gfortran make cmake \
          curl wget git ca-certificates file \
          software-properties-common

      # Apptainer (for NGPB solver)
      add-apt-repository -y ppa:apptainer/ppa
      apt-get update
      apt-get install -y apptainer

      touch /etc/mcce4-provisioned

  - mode: user
    script: |
      #!/bin/bash
      set -eux -o pipefail

      # Install Miniforge if conda is not available
      if ! command -v conda >/dev/null 2>&1; then
          if [ ! -d "$HOME/miniforge3" ]; then
              echo "Installing Miniforge..."
              ARCH=$(uname -m)
              curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh" \
                  -o /tmp/miniforge.sh
              bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
              rm /tmp/miniforge.sh
          fi
          eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
          conda init bash
      fi

cpus: 4
memory: "8GiB"
disk: "30GiB"
YAML
    echo "   ✅ Config written to $LIMA_CONFIG"
}

check_homebrew() {
    if command -v brew >/dev/null 2>&1; then
        echo "   ✅ Homebrew found"
        return 0
    fi
    echo "   ❌ Homebrew not found."
    echo '      Install: /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
    return 1
}

setup_lima() {
    echo ""
    echo "🍎 Setting up Lima Linux VM for macOS..."

    check_homebrew || return 1

    if ! command -v limactl >/dev/null 2>&1; then
        echo "   ⬇️  Installing Lima..."
        brew install lima
    else
        echo "   ✅ Lima found: $(limactl --version 2>/dev/null || echo 'installed')"
    fi

    generate_lima_config

    if limactl list 2>/dev/null | grep -q "$LIMA_INSTANCE"; then
        if ! limactl list 2>/dev/null | grep "$LIMA_INSTANCE" | grep -q "Running"; then
            echo "   ▶️  Starting existing Lima VM '$LIMA_INSTANCE'..."
            limactl start "$LIMA_INSTANCE"
        else
            echo "   ✅ Lima VM '$LIMA_INSTANCE' is already running."
        fi
    else
        echo "   🆕 Creating Lima VM '$LIMA_INSTANCE' (first time — may take a few minutes)..."
        limactl start --name="$LIMA_INSTANCE" "$LIMA_CONFIG"
    fi

    # Verify
    if limactl shell "$LIMA_INSTANCE" -- bash -c "command -v gcc && command -v apptainer" >/dev/null 2>&1; then
        echo "   ✅ Lima VM ready (gcc + apptainer available)."
        return 0
    else
        echo "   ❌ Lima VM missing required tools."
        echo "      Try: limactl delete $LIMA_INSTANCE && bash setup.sh --rebuild"
        return 1
    fi
}

build_mcce_in_lima() {
    echo ""
    echo "🔨 Building MCCE4 natively inside Lima VM..."

    # Create conda env if needed, then compile MCCE4
    limactl shell "$LIMA_INSTANCE" -- bash -c "
        set -e
        cd '$PROJECT_ROOT_DIR'

        # Initialize conda
        if [ -f \"\$HOME/miniforge3/etc/profile.d/conda.sh\" ]; then
            . \"\$HOME/miniforge3/etc/profile.d/conda.sh\"
        elif command -v conda >/dev/null 2>&1; then
            eval \"\$(conda shell.bash hook)\"
        else
            echo '❌ conda not found inside Lima VM'
            exit 1
        fi

        # Create or update conda env
        if conda env list | grep -q '^mc4 '; then
            echo '   Updating conda env mc4...'
            conda env update -n mc4 -f mc4.yml --prune -q
        else
            echo '   Creating conda env mc4...'
            conda env create -f mc4.yml -q
        fi

        conda activate mc4

        # Clean and rebuild
        if [ '$FORCE_REBUILD' -eq 1 ] || [ ! -f bin/mcce ]; then
            echo '   Compiling MCCE4...'
            make clean 2>/dev/null || true
            rm -f bin/mcce bin/delphi lib/*.o lib/*.a 2>/dev/null || true
            make
        else
            echo '   bin/mcce already exists. Use --rebuild to force.'
        fi

        # Verify
        if [ -f bin/mcce ]; then
            echo \"   ✅ mcce built: \$(file bin/mcce)\"
        else
            echo '   ❌ mcce was not produced by the build.'
            exit 1
        fi

        # Verify delphi_precompiled works
        if [ -f bin/delphi_precompiled ]; then
            echo \"   ✅ delphi_precompiled: \$(file bin/delphi_precompiled)\"
        fi
    "
}

build_ngpb_in_lima() {
    echo ""
    echo "🧪 Checking for NextGenPB solver..."

    local ngpb_sif="$BIN_DIR/NextGenPB_MCCE4.sif"

    if [[ -f "$ngpb_sif" && "$FORCE_REBUILD" -eq 0 ]]; then
        echo "   ✅ NextGenPB_MCCE4.sif already exists."
        echo "      Use --rebuild to force a fresh build."
        return 0
    fi

    if [[ ! -f "$BIN_DIR/recipe_MCCE.def" ]]; then
        echo "   ⚠️  bin/recipe_MCCE.def not found. Skipping NGPB build."
        echo "      DELPHI and APBS will still work."
        return 0
    fi

    echo "   Building NextGenPB_MCCE4.sif inside Lima VM..."
    echo "   (This is a one-time build — ~20-30 min)"

    if limactl shell "$LIMA_INSTANCE" -- bash -c "
        cd '$PROJECT_ROOT_DIR'
        apptainer build --force --fakeroot bin/NextGenPB_MCCE4.sif bin/recipe_MCCE.def
    "; then
        echo "   ✅ NextGenPB_MCCE4.sif built."
        return 0
    fi

    echo "   ⚠️  NGPB build failed. DELPHI and APBS will still work."
    echo "      You can build it later: mc4 --shell then run the apptainer build manually."
    return 0
}


# ═══════════════════════════════════════════════════════════════════════════════
#  LINUX: APPTAINER CONTAINER
# ═══════════════════════════════════════════════════════════════════════════════

install_apptainer() {
    echo ""
    echo "🔍 Checking for Apptainer..."

    if command -v apptainer >/dev/null 2>&1; then
        echo "   ✅ Apptainer found: $(which apptainer)"
        return 0
    fi

    echo "   ⚠️  Apptainer not found. Attempting user-local install..."

    if command -v curl >/dev/null 2>&1; then
        echo "   ⬇️  Downloading Apptainer to $APPTAINER_LOCAL_DIR ..."
        mkdir -p "$APPTAINER_LOCAL_DIR"
        curl -fsSL \
            https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh \
            | bash -s - "$APPTAINER_LOCAL_DIR"

        if "$APPTAINER_LOCAL_DIR/bin/apptainer" version >/dev/null 2>&1; then
            export PATH="$APPTAINER_LOCAL_DIR/bin:$PATH"
            update_shell_config "$APPTAINER_LOCAL_DIR/bin" "Apptainer (Unprivileged)"
            echo "   ✅ Apptainer installed at $APPTAINER_LOCAL_DIR"
            return 0
        else
            rm -rf "$APPTAINER_LOCAL_DIR"
        fi
    fi

    echo "   ❌ Could not install Apptainer."
    return 1
}

setup_conda_env() {
    echo ""
    echo "🐍 Setting up Conda environment '$CONDA_ENV_NAME'..."

    if ! command -v conda >/dev/null 2>&1; then
        echo "   ❌ 'conda' not found. Install Miniconda: https://docs.anaconda.com/miniconda/"
        return 1
    fi

    if ! conda list -n base 2>/dev/null | grep -q "conda-libmamba-solver"; then
        conda install -n base -c conda-forge conda-libmamba-solver -y --quiet
    fi
    conda config --set solver libmamba 2>/dev/null || true

    if conda env list | grep -q "^${CONDA_ENV_NAME}\s\+"; then
        echo "   Updating environment '$CONDA_ENV_NAME'..."
        yes "a" | conda env update -n "$CONDA_ENV_NAME" -f "$CONDA_YML" --prune
    else
        echo "   Creating environment '$CONDA_ENV_NAME'..."
        yes "a" | conda env create -f "$CONDA_YML"
    fi
}

build_container_linux() {
    echo ""
    echo "🔧 Building Apptainer container..."
    echo "   Definition : $MCCE4_ALPHA_DEF"
    echo "   Output     : $MCCE4_ALPHA_SIF"

    mkdir -p "$(dirname "$MCCE4_ALPHA_SIF")"

    if ! command -v apptainer >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)" 2>/dev/null
        conda activate "$CONDA_ENV_NAME" 2>/dev/null
    fi

    if apptainer build --force --fakeroot "$MCCE4_ALPHA_SIF" "$MCCE4_ALPHA_DEF" 2>&1; then
        echo "   ✅ Container built: $MCCE4_ALPHA_SIF"
        return 0
    fi

    echo "   ⚠️  --fakeroot failed. Retrying..."
    if apptainer build --force "$MCCE4_ALPHA_SIF" "$MCCE4_ALPHA_DEF" 2>&1; then
        echo "   ✅ Container built: $MCCE4_ALPHA_SIF"
        return 0
    fi

    echo "   ❌ Container build failed."
    return 1
}

should_rebuild_container() {
    [[ ! -f "$MCCE4_ALPHA_SIF" ]] && echo "   No existing container — building..." && return 0
    [[ "$FORCE_REBUILD" -eq 1 ]] && echo "   --rebuild flag — forcing rebuild..." && return 0
    [[ "$CONDA_YML" -nt "$MCCE4_ALPHA_SIF" ]] && echo "   ⚠️  mc4.yml newer — rebuilding..." && return 0
    [[ "$MCCE4_ALPHA_DEF" -nt "$MCCE4_ALPHA_SIF" ]] && echo "   ⚠️  .def newer — rebuilding..." && return 0
    echo "   ✅ Container up to date. Use --rebuild to force."
    return 1
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NATIVE FALLBACK (Linux without container, or macOS without Lima)
# ═══════════════════════════════════════════════════════════════════════════════

configure_compiler() {
    [[ "$OS" != "macos" ]] && return 0

    if ! command -v xcrun >/dev/null 2>&1; then
        echo "   ❌ xcrun not found. Run: xcode-select --install"
        return 1
    fi

    export CC CXX
    CC=$(xcrun -find cc 2>/dev/null)
    CXX=$(xcrun -find c++ 2>/dev/null)
    [[ -z "$CC" ]] && return 1

    SDK_PATH=$(xcrun --show-sdk-path 2>/dev/null)
    [[ -z "$SDK_PATH" ]] && return 1

    export CFLAGS="${CFLAGS:+$CFLAGS }-isysroot $SDK_PATH"
    export CXXFLAGS="${CXXFLAGS:+$CXXFLAGS }-isysroot $SDK_PATH"
    [[ "$ARCH" == "arm64" ]] && export CFLAGS="$CFLAGS -arch arm64" && export CXXFLAGS="$CXXFLAGS -arch arm64"

    CXX_INCLUDE_DIR=""
    if [[ -f "$SDK_PATH/usr/include/c++/v1/cstdlib" ]]; then
        CXX_INCLUDE_DIR="$SDK_PATH/usr/include/c++/v1"
    elif [[ -f "/Library/Developer/CommandLineTools/usr/include/c++/v1/cstdlib" ]]; then
        CXX_INCLUDE_DIR="/Library/Developer/CommandLineTools/usr/include/c++/v1"
    else
        local found
        found=$(find /Library/Developer/CommandLineTools/SDKs -path "*/include/c++/v1/cstdlib" 2>/dev/null | head -1)
        [[ -n "$found" ]] && CXX_INCLUDE_DIR="$(dirname "$found")"
    fi
    [[ -n "$CXX_INCLUDE_DIR" ]] && export CXXFLAGS="$CXXFLAGS -isystem $CXX_INCLUDE_DIR"

    export CFLAGS="$CFLAGS -Wno-implicit-function-declaration -Wno-format -Wno-error"
    export CXXFLAGS="$CXXFLAGS -Wno-implicit-function-declaration -Wno-format -Wno-error -stdlib=libc++"
}

build_native() {
    echo ""
    echo "🔨 Building MCCE binaries natively..."
    echo "   ⚠️  NGPB and delphi_precompiled will NOT work (requires Linux)."

    if ! command -v make >/dev/null 2>&1; then
        echo "   ❌ 'make' not found."
        return 1
    fi

    cd "$PROJECT_ROOT_DIR" || return 1

    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)" 2>/dev/null
        conda activate "$CONDA_ENV_NAME" 2>/dev/null
    fi

    configure_compiler || return 1

    local makefile_patched=0
    if [[ "$OS" == "macos" ]]; then
        cp Makefile Makefile._setup_bak
        local sysroot_flags="-isysroot $SDK_PATH"
        [[ "$ARCH" == "arm64" ]] && sysroot_flags="$sysroot_flags -arch arm64"
        local cxx_extra=""
        [[ -n "${CXX_INCLUDE_DIR:-}" ]] && cxx_extra="-isystem $CXX_INCLUDE_DIR"
        local gcc_rep="$CC $sysroot_flags -Wno-implicit-function-declaration -Wno-format -Wno-error"
        local gxx_rep="$CXX $sysroot_flags -stdlib=libc++ $cxx_extra -Wno-implicit-function-declaration -Wno-format -Wno-error"
        sed -i.tmp "s|g++ |${gxx_rep} |g" Makefile
        sed -i.tmp "s|gcc |${gcc_rep} |g" Makefile
        rm -f Makefile.tmp
        makefile_patched=1
    fi

    [[ "$FORCE_REBUILD" -eq 1 ]] && make clean 2>/dev/null || true

    make ${CC:+CC="$CC"} ${CXX:+CXX="$CXX"} ${CFLAGS:+CFLAGS="$CFLAGS"} ${CXXFLAGS:+CXXFLAGS="$CXXFLAGS"} || {
        [[ "$makefile_patched" -eq 1 ]] && mv Makefile._setup_bak Makefile
        return 1
    }

    [[ "$makefile_patched" -eq 1 ]] && mv Makefile._setup_bak Makefile
    echo "   ✅ Native build complete."
}


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🔧 MCCE4-Alpha Setup"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "   Clone root : $PROJECT_ROOT_DIR"
echo "   Install dir: $INSTALL_DIR"
[[ "$FORCE_REBUILD" -eq 1 ]] && echo "   Mode       : full rebuild (--rebuild)"
echo ""

detect_platform

INSTALL_MODE=""

# ── macOS: Lima VM with native build ──────────────────────────────────────────
if [[ "$OS" == "macos" ]]; then
    echo ""
    echo "🍎 macOS detected — setting up Lima Linux VM."
    echo "   Your code edits on Mac are instantly visible inside the VM."
    echo "   delphi_precompiled + NGPB work because you're running real Linux."
    echo ""

    step_start "Lima VM setup"
    if setup_lima; then
        step_end "Lima VM setup"
        INSTALL_MODE="lima"

        step_start "MCCE4 build (inside VM)"
        build_mcce_in_lima
        step_end "MCCE4 build (inside VM)"

        step_start "NGPB solver check"
        build_ngpb_in_lima
        step_end "NGPB solver check"
    else
        step_end "Lima VM setup"
        echo ""
        echo "⚠️  Lima setup failed. Falling back to native macOS build."
        echo "   NGPB and delphi_precompiled will NOT work."
        INSTALL_MODE="native"

        step_start "Conda environment"
        setup_conda_env || exit 1
        step_end "Conda environment"

        step_start "Native build"
        build_native || exit 1
        step_end "Native build"
    fi

# ── Linux: Apptainer container ───────────────────────────────────────────────
elif [[ "$OS" == "linux" ]]; then
    step_start "Apptainer install"
    if install_apptainer; then
        step_end "Apptainer install"
        INSTALL_MODE="container"

        step_start "Container build"
        if should_rebuild_container; then
            build_container_linux || {
                echo "⚠️  Container failed. Falling back to native."
                INSTALL_MODE="native"
            }
        fi
        step_end "Container build"
    else
        step_end "Apptainer install"
        INSTALL_MODE="native"
    fi

    if [[ "$INSTALL_MODE" == "native" ]]; then
        step_start "Conda environment"
        setup_conda_env || exit 1
        step_end "Conda environment"

        step_start "Native build"
        build_native || exit 1
        step_end "Native build"
    fi
fi

# ── Finalize ──────────────────────────────────────────────────────────────────
update_shell_config "$BIN_DIR" "MCCE4-Alpha CLI"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🚀 MCCE4-Alpha Setup Complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "   To activate, restart your shell or run:"
echo "     source $SHELL_CONFIG"
echo ""
case "$INSTALL_MODE" in
    lima)
        echo "   Mode       : Lima VM (macOS → Linux)"
        echo "   Lima VM    : $LIMA_INSTANCE"
        echo "   Solvers    : DELPHI ✅  APBS ✅  NGPB ✅"
        echo ""
        echo "   Code edits on Mac are instantly reflected — no rebuild needed."
        echo "   Only run setup.sh --rebuild if Makefile or mc4.yml changes."
        ;;
    container)
        echo "   Mode       : Apptainer container"
        echo "   Image      : $MCCE4_ALPHA_SIF"
        echo "   Solvers    : DELPHI ✅  APBS ✅  NGPB ✅"
        ;;
    native)
        echo "   Mode       : Native build"
        echo "   Solvers    : DELPHI ⚠️   APBS ⚠️   NGPB ❌"
        echo "   (Precompiled Linux binaries won't run natively on macOS)"
        ;;
esac
echo ""
echo "   Usage:"
echo "     mc4 step1.py prot.pdb"
echo "     mc4 step2.py"
echo "     mc4 step3.py"
echo "     mc4 step4.py"
echo "     mc4 --shell              # Interactive Linux shell"

print_timing_summary

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
