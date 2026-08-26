#!/bin/bash
# =============================================================================
#  transRNA pipeline v1.2 — login-node launcher
#
#  Usage:
#    bash run.sh                                    # fresh run, pipeline defaults (params.config)
#    bash run.sh --resume                            # resume a previous run (reuse cached steps)
#    bash run.sh --params-file config/params/X.yml   # project-specific overrides layered on top
#                                                      # of params.config (e.g. alvi_rnaome.yml)
#    bash run.sh --resume --params-file config/params/X.yml   # both together
# =============================================================================

set -euo pipefail

# ── Edit these two if needed ──────────────────────────────────────────────────
PARTITION="icelake"
PROJECT="MAORI-SL2-CPU"
# ─────────────────────────────────────────────────────────────────────────────

# Parse optional --resume and --params-file <path> flags
RESUME_FLAG=""
PARAMS_FILE_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --resume)
            RESUME_FLAG="-resume"
            shift
            ;;
        --params-file)
            PARAMS_FILE_ARGS=(-params-file "$2")
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# Suppress tput warnings from Nextflow in non-interactive sessions
export TERM=xterm

PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p "${PIPELINE_DIR}/logs" "${PIPELINE_DIR}/pipeline_reports"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG="${PIPELINE_DIR}/logs/nextflow_${TIMESTAMP}.log"

echo "============================================================"
echo "  transRNA pipeline v1.2"
echo "  $(date)"
echo "  Partition   : ${PARTITION}"
echo "  Project     : ${PROJECT}"
echo "  Mode        : ${RESUME_FLAG:-(fresh run)}"
echo "  Params file : ${PARAMS_FILE_ARGS[1]:-(pipeline defaults only — params.config)}"
echo "  Nextflow    : $(nextflow -version 2>&1 | head -1)"
echo "  Log         : ${LOG}"
echo "============================================================"

nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile cambridge \
    "${PARAMS_FILE_ARGS[@]}" \
    --partition "${PARTITION}" \
    --project   "${PROJECT}" \
    ${RESUME_FLAG} \
    2>&1 | tee "${LOG}"

echo "============================================================"
echo "  Finished: $(date)"
echo "  Log: ${LOG}"
echo "============================================================"
