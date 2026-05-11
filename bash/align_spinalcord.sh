#!/usr/bin/env bash
# =============================================================================
# Pipeline 2: Spinal Cord snRNA-seq + ATAC-seq alignment (E13.5 → 2yr)
# SraRunTable_(1).csv — SMARTer Stranded Total RNA + INTACT nuclei protocol
# RNA  → STAR paired-end alignment
# ATAC → Bowtie2 alignment + MACS2 peak calling
# =============================================================================
# Usage: nohup bash align_spinalcord_bulk.sh &
# =============================================================================

set -euo pipefail

# ── Load modules ──────────────────────────────────────────────────────────────
module load star/2.7.10a
module load bowtie2/2.4.5
module load samtools/1.14

# ── CONFIGURATION — edit these ────────────────────────────────────────────────
FASTQ_DIR="/home/nakagawa/datasets/SRR_spinalcord"
OUT_DIR="/home/nakagawa/datasets/aligned_spinalcord"
LOG_DIR="${OUT_DIR}/logs"
THREADS=32

# Genome — same mm10 index as cortex pipeline
GENOME_DIR="/home/nakagawa/datasets/genome/star_mm10_index"
GENOME_FA="/home/yukik/gtf/GRCm38/genome.fa"
GTF="/home/yukik/gtf/GRCm38/genes.gtf"

# Bowtie2 index for mm10 (built in Step 1b below)
BT2_INDEX="/home/nakagawa/datasets/genome/bowtie2_mm10/mm10"

LOGFILE="${LOG_DIR}/pipeline.log"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${OUT_DIR}"/{rna,atac,logs} \
         "/home/nakagawa/datasets/genome/bowtie2_mm10"

exec > >(tee -a "${LOGFILE}") 2>&1
echo "========================================"
echo " Spinal cord pipeline started: $(date)"
echo "========================================"


# =============================================================================
# STEP 1a: Build STAR genome index (shared with cortex pipeline — skip if done)
# =============================================================================
if [[ ! -f "${GENOME_DIR}/SA" ]]; then
    echo ""
    echo "--- Building STAR genome index ---"
    STAR \
        --runMode genomeGenerate \
        --runThreadN "${THREADS}" \
        --genomeDir "${GENOME_DIR}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GTF}" \
        --sjdbOverhang 151 \
        && echo "[OK] STAR index built" \
        || { echo "[FAIL] STAR index"; exit 1; }
else
    echo "[SKIP] STAR genome index already exists"
fi


# =============================================================================
# STEP 1b: Build Bowtie2 genome index for ATAC (skip if done)
# =============================================================================
if [[ ! -f "${BT2_INDEX}.1.bt2" ]]; then
    echo ""
    echo "--- Building Bowtie2 genome index ---"
    bowtie2-build \
        --threads "${THREADS}" \
        "${GENOME_FA}" \
        "${BT2_INDEX}" \
        && echo "[OK] Bowtie2 index built" \
        || { echo "[FAIL] Bowtie2 index"; exit 1; }
else
    echo "[SKIP] Bowtie2 index already exists"
fi


# =============================================================================
# STEP 2: RNA-seq alignment (snRNA — SMARTer Stranded Total RNA)
#
# This is paired-end bulk/pseudo-bulk nuclear RNA, NOT 10X droplet.
# No barcode demultiplexing needed — standard STAR paired-end alignment.
# Each SRR = one sample (one age/replicate).
#
# Samples and ages from SraRunTable_(1).csv:
# =============================================================================

declare -A RNA_SAMPLES=(
    # SRR          AGE         GENOTYPE
    [SRR18338966]="2yr_ChatCre_r1"
    [SRR18338967]="2yr_ChatCre_r2"
    [SRR18338968]="P56_ChatCre_r1"
    [SRR18338969]="P56_ChatCre_r2"
    [SRR18338970]="P56_ChatCre_r3"
    [SRR18338971]="P21_ChatCre_r1"
    [SRR18338972]="P21_ChatCre_r2"
    [SRR18338973]="P21_ChatCre_r3"
    [SRR18338974]="P13_ChatCre_r1"
    [SRR18338975]="P13_ChatCre_r2"
    [SRR18338976]="P4_ChatCre_r1"
    [SRR18338977]="P4_ChatCre_r2"
    [SRR18338978]="P4_ChatCre_r3"
    [SRR18338979]="E13.5_ChatCre_r1"
    [SRR18338980]="E13.5_ChatCre_r2"
    [SRR18338981]="E10.5_HB9GFP_r1"
    [SRR18338982]="E10.5_HB9GFP_r2"
    [SRR18338983]="E10.5_HB9GFP_r3"
)

echo ""
echo "--- RNA-seq alignment (snRNA SMARTer) ---"

for SRR in "${!RNA_SAMPLES[@]}"; do
    LABEL="${RNA_SAMPLES[$SRR]}"
    SAMPLE="${SRR}_${LABEL}"
    SAMPLE_OUT="${OUT_DIR}/rna/${SAMPLE}"

    if [[ -f "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam.bai" ]]; then
        echo "[SKIP] ${SAMPLE}"
        continue
    fi

    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "[WAIT] ${SAMPLE} — FASTQs not yet downloaded"
        continue
    fi

    echo "[ALIGN RNA] ${SAMPLE} ..."
    mkdir -p "${SAMPLE_OUT}"

    STAR \
        --runThreadN "${THREADS}" \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${R1}" "${R2}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM \
        --outFileNamePrefix "${SAMPLE_OUT}/" \
        --quantMode GeneCounts \
        --outBAMsortingThreadN 8 \
        --limitBAMsortRAM 60000000000 \
        --outFilterMultimapNmax 10 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        2>"${LOG_DIR}/${SAMPLE}_STAR.log" \
    && samtools index -@ 8 "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam" \
    && echo "[OK]   ${SAMPLE}" \
    || echo "[FAIL] ${SAMPLE} — check ${LOG_DIR}/${SAMPLE}_STAR.log"

done


# =============================================================================
# STEP 3: ATAC-seq alignment (Bowtie2)
#
# ATAC reads: short paired-end reads (~84bp from CSV AvgSpotLen)
# Goal: open chromatin regions per timepoint
# =============================================================================

declare -A ATAC_SAMPLES=(
    [SRR18338984]="2yr_ChatCre_r1"
    [SRR18338985]="2yr_ChatCre_r2"
    [SRR18338986]="P56_ChatCre_r1"
    [SRR18338987]="P56_ChatCre_r2"
    [SRR18338988]="P21_ChatCre_r1"
    [SRR18338989]="P21_ChatCre_r2"
    [SRR18338990]="P13_ChatCre_r1"
    [SRR18338991]="P13_ChatCre_r2"
    [SRR18338992]="P4_ChatCre_r1"
    [SRR18338993]="P4_ChatCre_r2"
    [SRR18338994]="E13.5_ChatCre_r1"
    [SRR18338995]="E13.5_ChatCre_r2"
    [SRR18338996]="E10.5_HB9GFP_r1"
    [SRR18338997]="E10.5_HB9GFP_r2"
)

echo ""
echo "--- ATAC-seq alignment (Bowtie2) ---"

for SRR in "${!ATAC_SAMPLES[@]}"; do
    LABEL="${ATAC_SAMPLES[$SRR]}"
    SAMPLE="${SRR}_${LABEL}"
    SAMPLE_OUT="${OUT_DIR}/atac/${SAMPLE}"

    if [[ -f "${SAMPLE_OUT}/${SAMPLE}_peaks.narrowPeak" ]]; then
        echo "[SKIP] ${SAMPLE}"
        continue
    fi

    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "[WAIT] ${SAMPLE} — FASTQs not yet downloaded"
        continue
    fi

    echo "[ALIGN ATAC] ${SAMPLE} ..."
    mkdir -p "${SAMPLE_OUT}"

    # Bowtie2 alignment — ATAC uses very sensitive local for short reads
    bowtie2 \
        -p "${THREADS}" \
        -x "${BT2_INDEX}" \
        -1 "${R1}" \
        -2 "${R2}" \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        --phred33 \
        -X 2000 \
        2>"${LOG_DIR}/${SAMPLE}_bowtie2.log" \
    | samtools view -bS -q 30 - \
    | samtools sort -@ 8 -o "${SAMPLE_OUT}/${SAMPLE}_sorted.bam" \
    && samtools index -@ 8 "${SAMPLE_OUT}/${SAMPLE}_sorted.bam" \
    && echo "[OK bowtie2] ${SAMPLE}" \
    || { echo "[FAIL bowtie2] ${SAMPLE}"; continue; }

    # Remove mitochondrial reads and duplicates
    samtools view -@ 8 -b \
        "${SAMPLE_OUT}/${SAMPLE}_sorted.bam" \
        $(samtools view -H "${SAMPLE_OUT}/${SAMPLE}_sorted.bam" \
          | grep "^@SQ" | grep -v "chrM" \
          | sed 's/.*SN://;s/\t.*//') \
        > "${SAMPLE_OUT}/${SAMPLE}_nochrM.bam"

    samtools markdup -@ 8 -r \
        "${SAMPLE_OUT}/${SAMPLE}_nochrM.bam" \
        "${SAMPLE_OUT}/${SAMPLE}_dedup.bam" \
        2>"${LOG_DIR}/${SAMPLE}_markdup.log"

    samtools index -@ 8 "${SAMPLE_OUT}/${SAMPLE}_dedup.bam"
    rm "${SAMPLE_OUT}/${SAMPLE}_nochrM.bam"

    # Tn5 offset shift: +4bp on + strand, -5bp on - strand
    # Convert to BED for MACS2
    samtools view -@ 8 -f 0x2 \
        "${SAMPLE_OUT}/${SAMPLE}_dedup.bam" \
    | awk 'BEGIN{OFS="\t"} {
        if ($2==99 || $2==83)  { print $3, $4+4,  $4+4+1,  $1, ".", "+" }
        if ($2==147 || $2==163){ print $3, $4-5,  $4-5+1,  $1, ".", "-" }
    }' \
    | sort -k1,1 -k2,2n \
    > "${SAMPLE_OUT}/${SAMPLE}_tn5shift.bed" \
    && echo "[OK Tn5 shift] ${SAMPLE}" \
    || { echo "[FAIL Tn5 shift] ${SAMPLE}"; continue; }

    echo "[OK ATAC prep] ${SAMPLE}"
done


# =============================================================================
# STEP 4: Peak calling with MACS2 (run after all ATAC alignments done)
# Check if macs2 is available; install if not
# =============================================================================
echo ""
echo "--- ATAC peak calling (MACS2) ---"

if ! command -v macs2 &>/dev/null; then
    echo "macs2 not found — installing into scrna conda env..."
    conda run -n scrna pip install macs2 --quiet \
        && echo "[OK] macs2 installed" \
        || echo "[FAIL] Install macs2 manually: conda activate scrna && pip install macs2"
fi

for SRR in "${!ATAC_SAMPLES[@]}"; do
    LABEL="${ATAC_SAMPLES[$SRR]}"
    SAMPLE="${SRR}_${LABEL}"
    SAMPLE_OUT="${OUT_DIR}/atac/${SAMPLE}"
    BED="${SAMPLE_OUT}/${SAMPLE}_tn5shift.bed"

    if [[ ! -f "${BED}" ]]; then
        echo "[SKIP peaks] ${SAMPLE} — BED not ready"
        continue
    fi
    if [[ -f "${SAMPLE_OUT}/${SAMPLE}_peaks.narrowPeak" ]]; then
        echo "[SKIP peaks] ${SAMPLE} already called"
        continue
    fi

    echo "[PEAKS] ${SAMPLE} ..."
    macs2 callpeak \
        -t "${BED}" \
        -f BED \
        -g mm \
        -n "${SAMPLE}" \
        --outdir "${SAMPLE_OUT}" \
        --nomodel \
        --shift -75 \
        --extsize 150 \
        --nolambda \
        --keep-dup all \
        -q 0.05 \
        2>"${LOG_DIR}/${SAMPLE}_macs2.log" \
    && echo "[OK peaks] ${SAMPLE}" \
    || echo "[FAIL peaks] ${SAMPLE} — check ${LOG_DIR}/${SAMPLE}_macs2.log"

done


# =============================================================================
# STEP 5: Summary
# =============================================================================
echo ""
echo "--- Alignment summary ---"
echo ""
echo "RNA samples:"
printf "%-45s %s\n" "Sample" "MappingRate"
for SRR in "${!RNA_SAMPLES[@]}"; do
    LABEL="${RNA_SAMPLES[$SRR]}"
    SAMPLE="${SRR}_${LABEL}"
    LOG="${LOG_DIR}/${SAMPLE}_STAR.log"
    if [[ -f "${LOG}" ]]; then
        RATE=$(grep "Uniquely mapped reads %" "${LOG}" \
               | awk '{print $NF}' || echo "NA")
        printf "%-45s %s\n" "${SAMPLE}" "${RATE}"
    fi
done

echo ""
echo "ATAC samples:"
printf "%-45s %s\n" "Sample" "OverallAlignRate"
for SRR in "${!ATAC_SAMPLES[@]}"; do
    LABEL="${ATAC_SAMPLES[$SRR]}"
    SAMPLE="${SRR}_${LABEL}"
    LOG="${LOG_DIR}/${SAMPLE}_bowtie2.log"
    if [[ -f "${LOG}" ]]; then
        RATE=$(grep "overall alignment rate" "${LOG}" \
               | awk '{print $1}' || echo "NA")
        printf "%-45s %s\n" "${SAMPLE}" "${RATE}"
    fi
done

echo ""
echo "========================================"
echo " Spinal cord pipeline finished: $(date)"
echo "========================================"

# =============================================================================
# OUTPUT STRUCTURE:
# aligned_spinalcord/
#   rna/<SRR>_<age>/
#     Aligned.sortedByCoord.out.bam      ← for IGV / downstream
#     ReadsPerGene.out.tab               ← gene count table → DESeq2
#   atac/<SRR>_<age>/
#     <sample>_dedup.bam                 ← cleaned ATAC BAM
#     <sample>_peaks.narrowPeak          ← open chromatin peaks → trajectory
#     <sample>_tn5shift.bed              ← Tn5-corrected cut sites
# =============================================================================
