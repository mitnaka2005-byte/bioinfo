#!/usr/bin/env bash
# =============================================================================
# Download all datasets: CellxGene .h5ad files + SRA FASTQs
# Usage: bash download_all_datasets.sh
# Then: Ctrl+Z -> bg -> disown
# =============================================================================

set -euo pipefail

# --------------------------------------------------------------------------
# CONFIGURATION — edit these before running
# --------------------------------------------------------------------------
OUTDIR="${HOME}/datasets"          # Change to your lab server path
THREADS=6                          # fasterq-dump parallel threads
LOG="${OUTDIR}/download.log"

mkdir -p "${OUTDIR}"/{h5ad,SRR_cortex,SRR_cortex2,SRR_spinalcord}

exec > >(tee -a "${LOG}") 2>&1
echo "========================================"
echo " Download started: $(date)"
echo "========================================"

# =============================================================================
# PART 1: CellxGene .h5ad direct downloads (wget)
# =============================================================================
echo ""
echo "--- Downloading CellxGene .h5ad files ---"

declare -A H5AD_FILES=(
    ["SMARTer_cells_MOp"]="https://datasets.cellxgene.cziscience.com/8964ca96-7583-4ce3-976f-512297c67be9.h5ad"
    ["DNA_Methylation_CHN"]="https://datasets.cellxgene.cziscience.com/294a2745-7102-448b-a55c-b6f01d774c27.h5ad"
    ["DNA_Methylation_CGN"]="https://datasets.cellxgene.cziscience.com/17e16125-ad0c-4ff1-96a4-c63175bf8de3.h5ad"
    ["10X_cells_v3_AIBS"]="https://datasets.cellxgene.cziscience.com/9c3db9b7-f0dc-4dc4-9e7d-e825b969e272.h5ad"
)

for NAME in "${!H5AD_FILES[@]}"; do
    URL="${H5AD_FILES[$NAME]}"
    DEST="${OUTDIR}/h5ad/${NAME}.h5ad"
    if [[ -f "${DEST}" ]]; then
        echo "[SKIP] ${NAME}.h5ad already exists"
    else
        echo "[DL]   ${NAME}.h5ad"
        wget -c --show-progress -O "${DEST}" "${URL}" \
            && echo "[OK]   ${NAME}.h5ad" \
            || echo "[FAIL] ${NAME}.h5ad — check log"
    fi
done

# =============================================================================
# PART 2: SRA — Mouse Primary Motor Cortex (SraRunTable.csv)
# RNA-Seq: SRR14407605-06-07, SRR12082755-764
# ATAC-seq: SRR12082772-774
# =============================================================================
echo ""
echo "--- fasterq-dump: Motor Cortex SRR runs ---"

CORTEX_SRRS=(
    # RNA-Seq
    SRR14407605 SRR14407606 SRR14407607
    SRR12082755 SRR12082756 SRR12082757 SRR12082758
    SRR12082759 SRR12082760 SRR12082761 SRR12082762
    SRR12082763 SRR12082764
    # ATAC-seq
    SRR12082772 SRR12082773 SRR12082774
)

for SRR in "${CORTEX_SRRS[@]}"; do
    DEST="${OUTDIR}/SRR_cortex/${SRR}"
    if ls "${OUTDIR}/SRR_cortex/${SRR}"*.fastq.gz 2>/dev/null | grep -q .; then
        echo "[SKIP] ${SRR} already downloaded"
    else
        echo "[DL]   ${SRR}"
        fasterq-dump "${SRR}" \
            --outdir "${OUTDIR}/SRR_cortex" \
            --threads "${THREADS}" \
            --progress \
            --split-files \
        && pigz -p "${THREADS}" "${OUTDIR}/SRR_cortex/${SRR}"*.fastq \
        && echo "[OK]   ${SRR}" \
        || echo "[FAIL] ${SRR} — check log"
    fi
done

# =============================================================================
# PART 3: SRA — Mouse Spinal Cord (SraRunTable_(1).csv)
# RNA-Seq: SRR18338966-83
# ATAC-seq: SRR18338984-97
# =============================================================================
echo ""
echo "--- fasterq-dump: Spinal Cord SRR runs ---"

SPINAL_SRRS=(
    # RNA-Seq
    SRR18338966 SRR18338967 SRR18338968 SRR18338969
    SRR18338970 SRR18338971 SRR18338972 SRR18338973
    SRR18338974 SRR18338975 SRR18338976 SRR18338977
    SRR18338978 SRR18338979 SRR18338980 SRR18338981
    SRR18338982 SRR18338983
    # ATAC-seq
    SRR18338984 SRR18338985 SRR18338986 SRR18338987
    SRR18338988 SRR18338989 SRR18338990 SRR18338991
    SRR18338992 SRR18338993 SRR18338994 SRR18338995
    SRR18338996 SRR18338997
)

for SRR in "${SPINAL_SRRS[@]}"; do
    if ls "${OUTDIR}/SRR_spinalcord/${SRR}"*.fastq.gz 2>/dev/null | grep -q .; then
        echo "[SKIP] ${SRR} already downloaded"
    else
        echo "[DL]   ${SRR}"
        fasterq-dump "${SRR}" \
            --outdir "${OUTDIR}/SRR_spinalcord" \
            --threads "${THREADS}" \
            --progress \
            --split-files \
        && pigz -p "${THREADS}" "${OUTDIR}/SRR_spinalcord/${SRR}"*.fastq \
        && echo "[OK]   ${SRR}" \
        || echo "[FAIL] ${SRR} — check log"
    fi
done

echo ""
echo "========================================"
echo " All downloads finished: $(date)"
echo "========================================"
