#!/usr/bin/env bash
# =============================================================================
# Retry failed SRR downloads with better settings
# Usage: bash retry_failed.sh
# =============================================================================

OUTDIR="${HOME}/datasets"          # Must match your original OUTDIR
THREADS=6
LOG="${OUTDIR}/retry.log"
MAX_RETRIES=3

exec > >(tee -a "${LOG}") 2>&1
echo "========================================"
echo " Retry started: $(date)"
echo "========================================"

# Extract failed SRRs from original log automatically
FAILED_SRRS=($(grep "\[FAIL\]" "${OUTDIR}/download.log" | grep -oP 'SRR\d+'))

if [[ ${#FAILED_SRRS[@]} -eq 0 ]]; then
    echo "No failed SRRs found in log. Exiting."
    exit 0
fi

echo "Found ${#FAILED_SRRS[@]} failed SRRs:"
printf '  %s\n' "${FAILED_SRRS[@]}"
echo ""

for SRR in "${FAILED_SRRS[@]}"; do
    # Determine output directory based on SRR prefix
    if [[ "${SRR}" == SRR183* ]]; then
        DEST_DIR="${OUTDIR}/SRR_spinalcord"
    else
        DEST_DIR="${OUTDIR}/SRR_cortex"
    fi

    # Skip if already successfully downloaded
    if ls "${DEST_DIR}/${SRR}"*.fastq.gz 2>/dev/null | grep -q .; then
        echo "[SKIP] ${SRR} already exists"
        continue
    fi

    echo "[RETRY] ${SRR} -> ${DEST_DIR}"
    SUCCESS=false

    for ATTEMPT in $(seq 1 ${MAX_RETRIES}); do
        echo "  Attempt ${ATTEMPT}/${MAX_RETRIES}..."
        sleep $((ATTEMPT * 10))   # backoff: 10s, 20s, 30s between retries

        fasterq-dump "${SRR}" \
            --outdir "${DEST_DIR}" \
            --threads "${THREADS}" \
            --progress \
            --split-files \
            --location NCBI \
        && pigz -p "${THREADS}" "${DEST_DIR}/${SRR}"*.fastq \
        && SUCCESS=true && break \
        || echo "  Attempt ${ATTEMPT} failed, retrying..."
    done

    if $SUCCESS; then
        echo "[OK]   ${SRR}"
    else
        echo "[FAIL] ${SRR} — gave up after ${MAX_RETRIES} attempts"
    fi
done

echo ""
echo "========================================"
echo " Retry finished: $(date)"
echo "========================================"

# Summary
echo ""
echo "--- Summary ---"
echo "Still failing:"
grep "\[FAIL\]" "${LOG}" | grep -oP 'SRR\d+' || echo "  None!"
