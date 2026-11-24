#!/bin/bash -l
# Build GROMACS 2026.4 with CUDA on a Lem (WCSS) H100 node.
# Lem GPU nodes: Intel Xeon Platinum 8462Y+ (x86_64, AVX-512), 64 cores, 4x H100-94GB.
# Installs two flavours into the same prefix:
#   gmx      (thread-MPI) - single-node runs, grompp, analysis tools
#   gmx_mpi  (real MPI)   - multi-node runs and -multidir replica ensembles
#
# BEFORE submitting, on ui.wcss.pl (don't rely on compute nodes having internet):
#   mkdir -p ~/src && cd ~/src
#   wget https://ftp.gromacs.org/gromacs/gromacs-2026.4.tar.gz
#   wget https://ftp.gromacs.org/regressiontests/regressiontests-2026.4.tar.gz
#   wget https://www.fftw.org/fftw-3.3.10.tar.gz          # keep as tarball, don't unpack
#   tar xzf gromacs-2026.4.tar.gz && tar xzf regressiontests-2026.4.tar.gz
#
# Pick modules (all from the same toolchain generation, e.g. all *-GCCcore-13.x):
#   source /usr/local/sbin/modules.sh
#   module avail -i gromacs     # check first - a recent GROMACS may already exist
#   module avail -i gcc ; module avail -i cuda ; module avail -i cmake ; module avail -i openmpi
#
# Then edit MODULES below and:  sbatch build_gromacs_2026.4_lem.sh
# (or build interactively: sub-interactive-lem-gpu, then: bash build_gromacs_2026.4_lem.sh)

#SBATCH --job-name=build-gmx-2026.4
#SBATCH --partition=lem-gpu          # routes to lem-gpu-short / lem-gpu-normal by time limit
#SBATCH --gres=gpu:hopper:1          # Lem requires an explicit GPU request
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # 64 cores / 4 GPUs per node
#SBATCH --time=03:00:00
#SBATCH --output=build-gmx-%j.log

set -eo pipefail   # no -u: module scripts often touch unset variables

### ---- edit these ----------------------------------------------------------
VER=2026.4
SRC=$HOME/src
PREFIX=$HOME/software/gromacs-$VER
CUDA_ARCH=90        # H100 = sm_90
SIMD=""             # empty = auto-detect (AVX_512 on Sapphire Rapids); AVX2_256 is worth a benchmark
# EasyBuild 2024a-generation toolchain + CUDA 12.8 (SYSTEM module, works with any GCC; needs driver >= 570)
MODULES="GCC/13.3.0 CMake/3.29.3-GCCcore-13.3.0 OpenMPI/5.0.3-GCC-13.3.0 CUDA/12.8.0"
# Alternative (GCC 14 host compiler is supported from CUDA 12.8 on):
# MODULES="GCC/14.3.0 CMake/3.31.8-GCCcore-14.3.0 OpenMPI/5.0.8-GCC-14.3.0 CUDA/12.8.0"
### ---------------------------------------------------------------------------

if [ "$MODULES" = "CHANGE_ME" ]; then
  echo "ERROR: set MODULES at the top of the script first"; exit 1
fi

# 1. Clean environment: conda envs inject their own compilers, libstdc++ and cmake
conda deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true
unset PYTHONPATH CPATH LIBRARY_PATH
source /usr/local/sbin/modules.sh
module load $MODULES

case "$(command -v g++) $(command -v cmake)" in
  *conda*|*mamba*|*miniforge*) echo "ERROR: conda toolchain in PATH - fix before building"; exit 1 ;;
esac

echo "== toolchain on $(hostname) ($(uname -m)) =="
command -v gcc g++ nvcc cmake mpicc || true
g++ --version | head -1
nvcc --version | tail -1
cmake --version | head -1
nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv || true

# 2. Configure, build, test, install (thread-MPI first, then MPI)
NPROC=${SLURM_CPUS_PER_TASK:-16}
SIMD_FLAG=${SIMD:+-DGMX_SIMD=$SIMD}

for MPI in OFF ON; do
  BUILD=$SRC/gromacs-$VER/build-cuda-mpi$MPI
  rm -rf "$BUILD" && mkdir -p "$BUILD" && cd "$BUILD"

  cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DGMX_GPU=CUDA \
    -DCMAKE_CUDA_ARCHITECTURES="$CUDA_ARCH" \
    -DGMX_MPI=$MPI \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_BUILD_OWN_FFTW_URL="file://$SRC/fftw-3.3.10.tar.gz" \
    -DREGRESSIONTEST_PATH="$SRC/regressiontests-$VER" \
    $SIMD_FLAG

  make -j "$NPROC"
  if [ "$MPI" = OFF ]; then
    make check -j "$NPROC"      # unit + regression tests, exercises the GPU
  fi
  make install
done

# 3. Sanity check
source "$PREFIX/bin/GMXRC"
gmx --version | grep -E "GROMACS version|SIMD instructions|GPU support|CUDA compiler|CUDA driver|CUDA runtime" || true
echo
echo "Done. In job scripts use:"
echo "  source /usr/local/sbin/modules.sh && module load $MODULES"
echo "  source $PREFIX/bin/GMXRC"