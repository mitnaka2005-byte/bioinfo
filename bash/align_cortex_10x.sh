#!/usr/bin/env bash
# =============================================================================
# Pipeline 1: Cortex 10X scRNA-seq alignment (E10.5 → P4)
# SraRunTable.csv — 10X Chromium v2 (E10.5-E18.5) and v3 (P4)
# Uses STARsolo for barcode demultiplexing + alignment in one step
# =============================================================================
# Usage: nohup bash align_cortex_10x.sh &
# =============================================================================

set -euo pipefail

# ── Load modules ──────────────────────────────────────────────────────────────
module load star/2.7.10a
module load samtools/1.14

# ── CONFIGURATION — edit these ────────────────────────────────────────────────
FASTQ_DIR="/home/nakagawa/datasets/SRR_cortex"
OUT_DIR="/home/nakagawa/datasets/aligned_cortex"
LOG_DIR="${OUT_DIR}/logs"
THREADS=32                        # conservative — leaves cores for others

# Genome — update paths after checking ls /home/yukik/gtf/GRCm38/
GENOME_DIR="/home/nakagawa/datasets/genome/star_mm10_index"   # will be built below
GENOME_FA="/home/yukik/gtf/GRCm38/genome.fa"                 # ← CHECK & EDIT
GTF="/home/yukik/gtf/GRCm38/genes.gtf"                       # ← CHECK & EDIT

# 10X barcode whitelists — download if missing (see note at bottom)
WL_V2="/home/nakagawa/datasets/genome/whitelist_10x_v2.txt"
WL_V3="/home/nakagawa/datasets/genome/whitelist_10x_v3.txt"

LOGFILE="${LOG_DIR}/pipeline.log"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${OUT_DIR}" "${LOG_DIR}" "${GENOME_DIR}" \
         "/home/nakagawa/datasets/genome"

exec > >(tee -a "${LOGFILE}") 2>&1
echo "========================================"
echo " Cortex 10X pipeline started: $(date)"
echo "========================================"


# =============================================================================
# STEP 1: Build STAR genome index (only runs if index doesn't exist)
# =============================================================================
if [[ ! -f "${GENOME_DIR}/SA" ]]; then
    echo ""
    echo "--- Building STAR genome index (mm10 / GRCm38) ---"
    echo "    This takes ~30 min and ~30GB RAM. Run once only."

    STAR \
        --runMode genomeGenerate \
        --runThreadN "${THREADS}" \
        --genomeDir "${GENOME_DIR}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GTF}" \
        --sjdbOverhang 130 \
        && echo "[OK] Genome index built" \
        || { echo "[FAIL] Genome index build failed"; exit 1; }
else
    echo "[SKIP] STAR genome index already exists"
fi


# =============================================================================
# STEP 2: Download 10X barcode whitelists if missing
# =============================================================================
if [[ ! -f "${WL_V2}" ]]; then
    echo ""
    echo "--- Downloading 10X v2 barcode whitelist ---"
    wget -q -O "${WL_V2}.gz" \
        "https://github.com/10XGenomics/cellranger/raw/main/lib/python/cellranger/barcodes/737K-august-2016.txt.gz" \
    && gunzip "${WL_V2}.gz" \
    && echo "[OK] v2 whitelist downloaded" \
    || echo "[FAIL] v2 whitelist — download manually"
fi

if [[ ! -f "${WL_V3}" ]]; then
    echo ""
    echo "--- Downloading 10X v3 barcode whitelist ---"
    wget -q -O "${WL_V3}.gz" \
        "https://github.com/10XGenomics/cellranger/raw/main/lib/python/cellranger/barcodes/3M-february-2018.txt.gz" \
    && gunzip "${WL_V3}.gz" \
    && echo "[OK] v3 whitelist downloaded" \
    || echo "[FAIL] v3 whitelist — download manually"
fi


# =============================================================================
# STEP 3: STARsolo alignment per sample
#
# Layout for 10X Chromium (both v2 and v3):
#   _1.fastq.gz = R1 = barcode + UMI (28bp for v2, 28bp for v3)
#   _2.fastq.gz = R2 = cDNA read (biological)
#
# v2: 16bp barcode + 10bp UMI = 26bp total (CB len=16, UMI len=10)
# v3: 16bp barcode + 12bp UMI = 28bp total (CB len=16, UMI len=12)
# =============================================================================

# Sample metadata: SRR -> (timepoint, chemistry)
# E10.5-E18.5 = v2 (older HiSeq runs)
# P4 = v3 (newer NovaSeq runs, SRR14407605-607)
declare -A TIMEPOINT=(
    [SRR14407605]="P4"
    [SRR14407606]="E17.5"
    [SRR14407607]="P4"
    [SRR12082755]="E11.5"
    [SRR12082756]="E12.5"
    [SRR12082757]="E13.5"
    [SRR12082758]="E14.5"
    [SRR12082759]="E15.5"
    [SRR12082760]="E16.5"
    [SRR12082761]="E18.5"
    [SRR12082762]="E18.5b"
    [SRR12082763]="P1"
    [SRR12082764]="P1b"
)

declare -A CHEMISTRY=(
    [SRR14407605]="v3"
    [SRR14407606]="v3"
    [SRR14407607]="v3"
    [SRR12082755]="v2"
    [SRR12082756]="v2"
    [SRR12082757]="v2"
    [SRR12082758]="v2"
    [SRR12082759]="v2"
    [SRR12082760]="v2"
    [SRR12082761]="v2"
    [SRR12082762]="v2"
    [SRR12082763]="v2"
    [SRR12082764]="v2"
)

echo ""
echo "--- Running STARsolo alignment ---"

for SRR in "${!TIMEPOINT[@]}"; do
    TP="${TIMEPOINT[$SRR]}"
    CHEM="${CHEMISTRY[$SRR]}"
    SAMPLE="${SRR}_${TP}"
    SAMPLE_OUT="${OUT_DIR}/${SAMPLE}"

    # Skip if already done
    if [[ -f "${SAMPLE_OUT}/Solo.out/GeneFull/filtered/matrix.mtx.gz" ]]; then
        echo "[SKIP] ${SAMPLE} already aligned"
        continue
    fi

    # Check FASTQs exist
    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "[WAIT] ${SAMPLE} — FASTQs not yet downloaded, skipping"
        continue
    fi

    echo "[ALIGN] ${SAMPLE} (10X ${CHEM}) ..."

    # Set chemistry-specific parameters
    if [[ "${CHEM}" == "v2" ]]; then
        WHITELIST="${WL_V2}"
        CB_LEN=16
        UMI_LEN=10
    else
        WHITELIST="${WL_V3}"
        CB_LEN=16
        UMI_LEN=12
    fi

    mkdir -p "${SAMPLE_OUT}"

    STAR \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "${WHITELIST}" \
        --soloCBstart 1 --soloCBlen ${CB_LEN} \
        --soloUMIstart $((CB_LEN + 1)) --soloUMIlen ${UMI_LEN} \
        --runThreadN "${THREADS}" \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${R2}" "${R1}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --outFileNamePrefix "${SAMPLE_OUT}/" \
        --soloFeatures GeneFull \
        --soloOutDir "${SAMPLE_OUT}/Solo.out" \
        --soloCellFilter EmptyDrops_CR \
        --soloOutQuant GeneFull \
        --outBAMsortingThreadN 8 \
        --limitBAMsortRAM 60000000000 \
        2>"${LOG_DIR}/${SAMPLE}_STAR.log" \
    && echo "[OK]    ${SAMPLE}" \
    || { echo "[FAIL]  ${SAMPLE} — check ${LOG_DIR}/${SAMPLE}_STAR.log"; continue; }

    # Index BAM
    samtools index -@ 8 "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam"

    # Compress STARsolo output matrices
    gzip -f "${SAMPLE_OUT}/Solo.out/GeneFull/filtered/matrix.mtx"
    gzip -f "${SAMPLE_OUT}/Solo.out/GeneFull/filtered/features.tsv"
    gzip -f "${SAMPLE_OUT}/Solo.out/GeneFull/filtered/barcodes.tsv"

done


# =============================================================================
# STEP 4: Generate summary table
# =============================================================================
echo ""
echo "--- Alignment summary ---"
echo "Sample,Timepoint,Chemistry,CellsDetected,MappingRate" \
    > "${OUT_DIR}/alignment_summary.csv"

for SRR in "${!TIMEPOINT[@]}"; do
    TP="${TIMEPOINT[$SRR]}"
    CHEM="${CHEMISTRY[$SRR]}"
    SAMPLE="${SRR}_${TP}"
    LOG="${LOG_DIR}/${SAMPLE}_STAR.log"

    if [[ -f "${LOG}" ]]; then
        CELLS=$(grep "Estimated Number of Cells" "${LOG}" \
                | awk '{print $NF}' || echo "NA")
        MAPRATE=$(grep "Uniquely mapped reads %" "${LOG}" \
                  | awk '{print $NF}' || echo "NA")
        echo "${SAMPLE},${TP},${CHEM},${CELLS},${MAPRATE}" \
            >> "${OUT_DIR}/alignment_summary.csv"
    fi
done

cat "${OUT_DIR}/alignment_summary.csv"

echo ""
echo "========================================"
echo " Cortex 10X pipeline finished: $(date)"
echo "========================================"

# =============================================================================
# NOTES:
# - Output count matrices: ${OUT_DIR}/<SAMPLE>/Solo.out/GeneFull/filtered/
#   matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
# - Load into scanpy: sc.read_10x_mtx(path)
# - If whitelist wget fails (network restriction), get from:
#   https://kb.10xgenomics.com/hc/en-us/articles/115004506263
# =============================================================================
