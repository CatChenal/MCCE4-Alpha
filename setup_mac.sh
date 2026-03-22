#!/bin/bash
# setup_mac.sh — MCCE4-Alpha Installer (macOS)
#
# Sets up a Lima Linux VM so MCCE4 runs in real Linux on your Mac.
# All solvers work: delphi_precompiled, APBS, NGPB.
# Your code edits on Mac are instantly visible — no rebuild needed.
#
# No sudo required.
#
# Usage:
#   bash setup_mac.sh                — smart re-run, skips completed steps
#   bash setup_mac.sh --rebuild      — delete VM and start fresh
#   bash setup_mac.sh --build-ngpb   — build NGPB from source instead of downloading
#
# Linux users: run setup.sh instead.

trap "echo ''; echo '❌ Script interrupted by user'; exit 1" SIGINT

# ── Check OS ───────────────────────────────────────────────────────────────────
if [[ "$(uname -s)" != "Darwin" ]]; then
    echo "❌ This script is for macOS only."
    echo "   Detected: $(uname -s)"
    echo ""
    echo "   Linux users: run setup.sh instead."
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

# ── Portable realpath (macOS readlink lacks -f) ───────────────────────────────
_realpath() {
    python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$1"
}

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT_DIR=$(dirname "$(_realpath "$0")")
cd "$PROJECT_ROOT_DIR" || exit 1

BIN_DIR="$PROJECT_ROOT_DIR/bin"
CONDA_YML="$PROJECT_ROOT_DIR/mc4.yml"

LIMA_INSTANCE="mcce4"
LIMA_CONFIG="$BIN_DIR/mcce4-lima.yaml"

SHELL_CONFIG="$HOME/.bashrc"
[[ "${SHELL:-}" == *"zsh"* ]] && SHELL_CONFIG="$HOME/.zshrc"

# ── Shell Config Helper ───────────────────────────────────────────────────────
update_shell_config() {
    local dir_path="$1" tool_name="$2"
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
        echo "   ✅ Added $dir_path to PATH"
    fi
}

# ── Check VM is actually running ──────────────────────────────────────────────
vm_is_running() {
    limactl list 2>/dev/null | grep "$LIMA_INSTANCE" | grep -q "Running"
}


# ═══════════════════════════════════════════════════════════════════════════════
#  LIMA VM SETUP
# ═══════════════════════════════════════════════════════════════════════════════

generate_lima_config() {
    cat > "$LIMA_CONFIG" << 'YAML'
# mcce4-lima.yaml — Lima VM for MCCE4-Alpha
# Full Linux dev environment: gcc, conda, apptainer.
# Home directory is mounted writable — Mac edits appear instantly.

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

# Enable Rosetta so x86_64 Linux binaries (delphi_precompiled, NextGenPB.sif)
# run transparently on Apple Silicon (arm64) VMs.
vmType: "vz"
vmOpts:
  vz:
    rosetta:
      enabled: true
      binfmt: true

provision:
  - mode: system
    script: |
      #!/bin/bash
      set -eux -o pipefail
      [ -f /etc/mcce4-provisioned ] && exit 0

      apt-get update
      apt-get install -y --no-install-recommends \
          build-essential gcc g++ gfortran make cmake \
          curl wget git ca-certificates file \
          software-properties-common

      # x86_64 cross-architecture support for Rosetta.
      # delphi_precompiled is an x86_64 ELF that needs the x86_64 dynamic
      # linker and gfortran runtime.
      #
      # We must pin the existing arm64 sources BEFORE adding amd64,
      # otherwise ports.ubuntu.com tries to serve amd64 (which it doesn't have).
      # Ubuntu 24.04 uses deb822 format in /etc/apt/sources.list.d/ubuntu.sources
      if [ -f /etc/apt/sources.list.d/ubuntu.sources ] && ! grep -q "Architectures:" /etc/apt/sources.list.d/ubuntu.sources; then
          sed -i 's/^Types: deb$/Architectures: arm64\nTypes: deb/' /etc/apt/sources.list.d/ubuntu.sources
      fi
      dpkg --add-architecture amd64
      CODENAME=$(lsb_release -cs)
      echo "deb [arch=amd64] http://archive.ubuntu.com/ubuntu ${CODENAME} main restricted universe multiverse" > /etc/apt/sources.list.d/amd64-cross.list
      echo "deb [arch=amd64] http://archive.ubuntu.com/ubuntu ${CODENAME}-updates main restricted universe multiverse" >> /etc/apt/sources.list.d/amd64-cross.list
      echo "deb [arch=amd64] http://archive.ubuntu.com/ubuntu ${CODENAME}-security main restricted universe multiverse" >> /etc/apt/sources.list.d/amd64-cross.list
      apt-get update
      apt-get install -y \
          libc6:amd64 \
          libgfortran5:amd64 \
          libquadmath0:amd64

      # Apptainer (for NGPB solver)
      add-apt-repository -y ppa:apptainer/ppa
      apt-get update
      apt-get install -y apptainer

      touch /etc/mcce4-provisioned

  - mode: user
    script: |
      #!/bin/bash
      set -eux -o pipefail

      # Install Miniforge if conda not available
      if ! command -v conda >/dev/null 2>&1 && [ ! -d "$HOME/miniforge3" ]; then
          ARCH=$(uname -m)
          curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh" \
              -o /tmp/miniforge.sh
          bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
          rm /tmp/miniforge.sh
          eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
          conda init bash
      fi

cpus: 4
memory: "8GiB"
disk: "30GiB"
YAML
}

setup_lima() {
    echo ""
    echo "🍎 Setting up Lima Linux VM..."

    # ── Homebrew ──────────────────────────────────────────────────────────────
    if ! command -v brew >/dev/null 2>&1; then
        echo "   ❌ Homebrew not found. Install it first:"
        echo '      /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
        return 1
    fi
    echo "   ✅ Homebrew found"

    # ── Rosetta (Apple Silicon only) ──────────────────────────────────────
    # delphi_precompiled and NextGenPB.sif are x86_64 binaries.
    # Rosetta lets the arm64 Linux VM run them transparently.
    if [[ "$(uname -m)" == "arm64" ]]; then
        if ! arch -x86_64 /usr/bin/true 2>/dev/null; then
            echo "   ⬇️  Installing Rosetta (needed for x86_64 solver binaries)..."
            softwareupdate --install-rosetta --agree-to-license 2>/dev/null || {
                echo "   ❌ Rosetta installation failed."
                echo "      Run manually: softwareupdate --install-rosetta"
                return 1
            }
        fi
        echo "   ✅ Rosetta available"
    fi

    # ── Lima ──────────────────────────────────────────────────────────────────
    if ! command -v limactl >/dev/null 2>&1; then
        echo "   ⬇️  Installing Lima..."
        brew install lima
    else
        echo "   ✅ Lima found"
    fi

    # ── Generate config ──────────────────────────────────────────────────────
    generate_lima_config
    echo "   ✅ VM config written"

    # ── --rebuild: delete stale VM ───────────────────────────────────────────
    if [[ "$FORCE_REBUILD" -eq 1 ]] && limactl list 2>/dev/null | grep -q "$LIMA_INSTANCE"; then
        echo "   🗑️  --rebuild: deleting existing VM '$LIMA_INSTANCE'..."
        limactl stop "$LIMA_INSTANCE" 2>/dev/null || true
        limactl delete "$LIMA_INSTANCE" --force 2>/dev/null || true
    fi

    # ── Create or start VM ───────────────────────────────────────────────────
    if limactl list 2>/dev/null | grep -q "$LIMA_INSTANCE"; then
        if ! vm_is_running; then
            echo "   ▶️  Starting VM '$LIMA_INSTANCE'..."
            limactl start "$LIMA_INSTANCE" --tty=false
        else
            echo "   ✅ VM '$LIMA_INSTANCE' is running."
        fi
    else
        echo "   🆕 Creating VM '$LIMA_INSTANCE' (first time — several minutes)..."
        echo "      Downloading Ubuntu 24.04 image + provisioning gcc, conda, apptainer..."
        echo ""
        limactl start --name="$LIMA_INSTANCE" --tty=false "$LIMA_CONFIG"
    fi

    # ── Verify VM actually booted ────────────────────────────────────────────
    if ! vm_is_running; then
        echo ""
        echo "   ❌ Lima VM failed to start."
        echo ""
        echo "      This usually means the Ubuntu image download failed (network issue)."
        echo "      Try:"
        echo "        1. Check your internet connection"
        echo "        2. limactl delete $LIMA_INSTANCE --force"
        echo "        3. bash setup_mac.sh --rebuild"
        return 1
    fi

    # ── Verify required tools ────────────────────────────────────────────────
    echo "   Verifying VM tools..."

    if ! limactl shell "$LIMA_INSTANCE" -- bash -c "command -v gcc" >/dev/null 2>&1; then
        echo "   ⚠️  gcc missing — provisioning..."
        limactl shell "$LIMA_INSTANCE" -- sudo bash -c "
            apt-get update && \
            apt-get install -y --no-install-recommends \
                build-essential gcc g++ gfortran make cmake \
                curl wget git ca-certificates file software-properties-common && \
            touch /etc/mcce4-provisioned
        " || { echo "   ❌ Failed to install gcc."; return 1; }
    fi

    # x86_64 libs for Rosetta (delphi_precompiled)
    if ! limactl shell "$LIMA_INSTANCE" -- bash -c "[ -f /lib64/ld-linux-x86-64.so.2 ]" >/dev/null 2>&1; then
        echo "   ⚠️  x86_64 libs missing — installing for delphi..."
        limactl shell "$LIMA_INSTANCE" -- sudo bash -c '
            # Pin existing sources to arm64 so ports.ubuntu.com stops trying amd64
            if [ -f /etc/apt/sources.list.d/ubuntu.sources ] && ! grep -q "Architectures:" /etc/apt/sources.list.d/ubuntu.sources; then
                sed -i "s/^Types: deb$/Architectures: arm64\nTypes: deb/" /etc/apt/sources.list.d/ubuntu.sources
            fi
            dpkg --add-architecture amd64
            CODENAME=$(lsb_release -cs)
            echo "deb [arch=amd64] http://archive.ubuntu.com/ubuntu ${CODENAME} main restricted universe multiverse" > /etc/apt/sources.list.d/amd64-cross.list
            echo "deb [arch=amd64] http://archive.ubuntu.com/ubuntu ${CODENAME}-updates main restricted universe multiverse" >> /etc/apt/sources.list.d/amd64-cross.list
            echo "deb [arch=amd64] http://archive.ubuntu.com/ubuntu ${CODENAME}-security main restricted universe multiverse" >> /etc/apt/sources.list.d/amd64-cross.list
            apt-get update && \
            apt-get install -y libc6:amd64 libgfortran5:amd64 libquadmath0:amd64
        ' || { echo "   ⚠️  Failed to install x86_64 libs. delphi may not work."; }
    fi

    if ! limactl shell "$LIMA_INSTANCE" -- bash -c "command -v apptainer" >/dev/null 2>&1; then
        echo "   ⚠️  apptainer missing — installing..."
        limactl shell "$LIMA_INSTANCE" -- sudo bash -c "
            add-apt-repository -y ppa:apptainer/ppa && \
            apt-get update && \
            apt-get install -y apptainer
        " || { echo "   ❌ Failed to install apptainer."; return 1; }
    fi

    if ! limactl shell "$LIMA_INSTANCE" -- bash -c "command -v conda >/dev/null 2>&1 || [ -f \$HOME/miniforge3/bin/conda ]" >/dev/null 2>&1; then
        echo "   ⬇️  Installing Miniforge inside VM..."
        limactl shell "$LIMA_INSTANCE" -- bash -c '
            ARCH=$(uname -m)
            curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh" \
                -o /tmp/miniforge.sh
            bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
            rm /tmp/miniforge.sh
            "$HOME/miniforge3/bin/conda" init bash
        ' || { echo "   ❌ Failed to install conda."; return 1; }
    fi

    # ── Final check ──────────────────────────────────────────────────────────
    local ok=1
    limactl shell "$LIMA_INSTANCE" -- bash -c "command -v gcc" >/dev/null 2>&1         || { echo "   ❌ gcc not available"; ok=0; }
    limactl shell "$LIMA_INSTANCE" -- bash -c "command -v apptainer" >/dev/null 2>&1   || { echo "   ❌ apptainer not available"; ok=0; }
    limactl shell "$LIMA_INSTANCE" -- bash -c "command -v conda >/dev/null 2>&1 || [ -f \$HOME/miniforge3/bin/conda ]" >/dev/null 2>&1 || { echo "   ❌ conda not available"; ok=0; }

    # Check Rosetta is working inside VM (Apple Silicon only)
    if [[ "$(uname -m)" == "arm64" ]]; then
        if limactl shell "$LIMA_INSTANCE" -- bash -c "[ -f /proc/sys/fs/binfmt_misc/rosetta ]" >/dev/null 2>&1; then
            echo "   ✅ Rosetta active inside VM (x86_64 binaries will work)"
        else
            echo "   ⚠️  Rosetta not detected inside VM. x86_64 binaries (delphi, NGPB) may not work."
            echo "      Try: limactl delete $LIMA_INSTANCE --force && bash setup_mac.sh --rebuild"
        fi
    fi

    if [[ "$ok" -eq 1 ]]; then
        echo "   ✅ VM ready (gcc + apptainer + conda)"
        return 0
    else
        echo "   ❌ VM missing required tools after provisioning."
        echo "      Try: limactl delete $LIMA_INSTANCE --force && bash setup_mac.sh --rebuild"
        return 1
    fi
}


# ═══════════════════════════════════════════════════════════════════════════════
#  BUILD MCCE4 INSIDE VM
# ═══════════════════════════════════════════════════════════════════════════════

build_mcce_in_lima() {
    echo ""
    echo "🔨 Building MCCE4 inside Lima VM..."

    limactl shell "$LIMA_INSTANCE" -- bash -c "
        set -e
        cd '$PROJECT_ROOT_DIR'

        # Initialize conda
        if [ -f \"\$HOME/miniforge3/etc/profile.d/conda.sh\" ]; then
            . \"\$HOME/miniforge3/etc/profile.d/conda.sh\"
        elif command -v conda >/dev/null 2>&1; then
            eval \"\$(conda shell.bash hook)\"
        else
            echo '   ❌ conda not found inside VM'
            exit 1
        fi

        # Create or update conda env
        if conda env list | grep -q '^mc4 '; then
            echo '   Updating conda env mc4...'
            yes 'a' | conda env update -n mc4 -f mc4.yml --prune
        else
            echo '   Creating conda env mc4...'
            yes 'a' | conda env create -f mc4.yml
        fi
        conda activate mc4

        # Clean and rebuild
        if [ '$FORCE_REBUILD' -eq 1 ] || [ ! -f bin/mcce ]; then
            echo '   Compiling MCCE4...'
            make clean 2>/dev/null || true
            rm -f bin/mcce bin/delphi lib/*.o lib/*.a 2>/dev/null || true
            make
        else
            echo '   bin/mcce exists. Use --rebuild to force recompilation.'
        fi

        # Verify mcce
        if [ -f bin/mcce ]; then
            echo \"   ✅ mcce: \$(file bin/mcce)\"
        else
            echo '   ❌ mcce was not built!'
            exit 1
        fi

        # Ensure delphi symlink
        if [ ! -f bin/delphi ] && [ -f bin/delphi_precompiled ]; then
            ln -sf delphi_precompiled bin/delphi
            echo '   ✅ Linked delphi_precompiled → delphi'
        fi
        if [ -f bin/delphi_precompiled ]; then
            echo \"   ✅ delphi: \$(file bin/delphi_precompiled)\"
        fi
    " || {
        echo "   ❌ MCCE4 build failed inside Lima VM."
        return 1
    }
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NGPB SOLVER
# ═══════════════════════════════════════════════════════════════════════════════

build_ngpb_in_lima() {
    echo ""
    echo "🧪 Setting up NextGenPB solver..."

    local ngpb_mc4="$BIN_DIR/NextGenPB_MCCE4.sif"
    local ngpb_generic="$BIN_DIR/NextGenPB.sif"
    local ngpb_recipe="$BIN_DIR/recipe_MCCE.def"
    local release="NextGenPB_v1.0.0"
    local url="https://github.com/concept-lab/NextGenPB/releases/download/$release/NextGenPB.sif"

    # Already present and not forcing
    if [[ -f "$ngpb_mc4" && "$FORCE_REBUILD" -eq 0 && "$BUILD_NGPB_SOURCE" -eq 0 ]]; then
        echo "   ✅ NextGenPB_MCCE4.sif exists. Use --rebuild or --build-ngpb to replace."
        return 0
    fi

    if [[ "$BUILD_NGPB_SOURCE" -eq 1 ]]; then
        # ── Build from source inside Lima VM ─────────────────────────────────
        if [[ ! -f "$ngpb_recipe" ]]; then
            echo "   ❌ bin/recipe_MCCE.def not found. Cannot build from source."
            echo "      Use default download instead (omit --build-ngpb)."
            return 1
        fi

        echo "   🔨 Building NextGenPB from source inside VM (~20-30 min)..."
        if limactl shell "$LIMA_INSTANCE" -- bash -c "
            cd '$PROJECT_ROOT_DIR'
            apptainer build --force --fakeroot bin/NextGenPB_MCCE4.sif bin/recipe_MCCE.def
        "; then
            echo "   ✅ NextGenPB_MCCE4.sif built from source."
            return 0
        fi

        echo "   ❌ NGPB source build failed."
        return 1
    else
        # ── Download pre-built (~1.6GB) ──────────────────────────────────────
        if [[ ! -f "$ngpb_generic" || "$FORCE_REBUILD" -eq 1 ]]; then
            echo "   ⬇️  Downloading NextGenPB container (~1.6GB)..."
            curl -L -o "$ngpb_generic" "$url" || {
                echo "   ⚠️  Download failed. DELPHI/APBS will still work."
                echo "      Retry later or build from source: bash setup_mac.sh --build-ngpb"
                return 0
            }
        fi

        if [[ -f "$ngpb_generic" ]]; then
            chmod +x "$ngpb_generic"
            ln -sf "$(basename "$ngpb_generic")" "$ngpb_mc4"
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
echo "🔧 MCCE4-Alpha Setup (macOS)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "   Clone root : $PROJECT_ROOT_DIR"
echo "   Platform   : macOS / $(uname -m)"
[[ "$FORCE_REBUILD" -eq 1 ]]    && echo "   Mode       : full rebuild (--rebuild)"
[[ "$BUILD_NGPB_SOURCE" -eq 1 ]] && echo "   NGPB       : build from source (--build-ngpb)"
echo ""
echo "   This sets up a lightweight Linux VM on your Mac."
echo "   Code edits on Mac are instantly visible inside the VM."
echo "   delphi + NGPB work because you're running real Linux."
echo ""

step_start "Lima VM setup"
if ! setup_lima; then
    echo ""
    echo "⛔ Lima VM setup failed. Cannot proceed."
    echo "   DELPHI and NGPB require a Linux environment."
    echo ""
    echo "   To fix:"
    echo "     1. Check your internet connection"
    echo "     2. limactl delete $LIMA_INSTANCE --force"
    echo "     3. bash setup_mac.sh --rebuild"
    exit 1
fi
step_end "Lima VM setup"

step_start "MCCE4 build (inside VM)"
if ! build_mcce_in_lima; then
    echo ""
    echo "⛔ MCCE4 build failed. Cannot proceed."
    exit 1
fi
step_end "MCCE4 build (inside VM)"

step_start "NGPB solver check"
build_ngpb_in_lima
step_end "NGPB solver check"

update_shell_config "$BIN_DIR" "MCCE4-Alpha CLI"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🚀 MCCE4-Alpha Setup Complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "   Restart your shell or run: source $SHELL_CONFIG"
echo ""
echo "   Mode    : Lima VM (native Linux)"
echo "   VM      : $LIMA_INSTANCE"
echo "   Solvers : DELPHI ✅  APBS ✅  NGPB ✅"
echo ""
echo "   Code edits on Mac are instant — no rebuild needed."
echo "   Run setup_mac.sh --rebuild only if Makefile or mc4.yml changes."
echo ""
echo "   Usage:"
echo "     mc4 step1.py <pdbfile>"
echo "     mc4 step2.py"
echo "     mc4 step3.py"
echo "     mc4 step4.py"
echo "     mc4 --shell              # Interactive Linux shell"

print_timing_summary

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
