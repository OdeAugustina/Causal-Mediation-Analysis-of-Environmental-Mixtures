#!/bin/bash
#SBATCH --job-name=cma_deepmed
#SBATCH --array=1-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120g
#SBATCH --time=10-00:00:00
#SBATCH --partition=l40-gpu
#SBATCH --qos=gpu_access
#SBATCH --gres=gpu:1
#SBATCH --output=logs/deepmed_%A_%a.out
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=aoodediran@aggies.ncat.edu
set -euo pipefail
module purge
module load r/4.5.0
module load python/3.12.4
module load cuda/12.4
cd "${SLURM_SUBMIT_DIR}"

# --- Activate the DeepMed Python venv ---
DEEPMED_VENV="/work/users/o/d/odeaug/envs/deepmed_py"
if [ -d "${DEEPMED_VENV}" ]; then
  source "${DEEPMED_VENV}/bin/activate"
  export RETICULATE_PYTHON="${DEEPMED_VENV}/bin/python"
else
  echo "ERROR: DeepMed venv not found at ${DEEPMED_VENV}" >&2
  exit 1
fi

# --- CRITICAL: put the venv's pip-installed CUDA libraries on the loader path ---
# TensorFlow 2.16's GPU libs (cudnn, cublas, cufft, ...) are installed as
# nvidia-*-cu12 pip packages under the venv. Without this, TF cannot dlopen
# them, silently falls back to CPU, and every rep takes ~2 hours instead of
# seconds. This export is what makes DeepMed actually use the GPU.
NVIDIA_LIBS=$("${DEEPMED_VENV}/bin/python" -c "import os, nvidia; b=os.path.dirname(nvidia.__file__); print(':'.join(os.path.join(b,d,'lib') for d in os.listdir(b) if os.path.isdir(os.path.join(b,d,'lib'))))")
export LD_LIBRARY_PATH="${NVIDIA_LIBS}:${LD_LIBRARY_PATH:-}"

echo "===== DeepMed job ====="
echo "Host:           $(hostname)"
echo "Job ID:         ${SLURM_JOB_ID}"
echo "Array Task:     ${SLURM_ARRAY_TASK_ID}"
echo "GPU:            ${CUDA_VISIBLE_DEVICES:-unset}"
echo "Python:         ${RETICULATE_PYTHON:-default}"
echo "Working Dir:    $(pwd)"
echo "Start:          $(date)"
echo

# --- Sanity check: confirm TF sees the GPU BEFORE running the scenario ---
echo "--- GPU visibility check (through the venv python) ---"
python -c "import tensorflow as tf; g=tf.config.list_physical_devices('GPU'); print('GPUs seen:', g); import sys; sys.exit(0 if g else 1)" \
  && echo "GPU OK" || echo "WARNING: TF sees NO GPU — job will run on CPU (slow). Check LD_LIBRARY_PATH."
echo

Rscript R/run_deepmed.R "${SLURM_ARRAY_TASK_ID}"
echo
echo "End: $(date)"
