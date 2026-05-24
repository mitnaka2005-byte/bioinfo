#!/usr/bin/env python
# coding: utf-8

# # FASTQ to h5ad Pipeline (kallisto | bustools)
# This pipeline aligns raw single-cell RNA-seq FASTQ files to the mm10 mouse reference genome and generates `.h5ad` count matrices.
# 
# ### Steps:
# 1. **Reference Generation:** Build the transcriptome index from the genome `.fa` and `.gtf`.
# 2. **Quantification:** Align reads and count UMIs to generate cell x gene matrices.
# 

# In[2]:


import os
import glob

# 1. Define Base Directories
base_dir = '/home/nakagawa/datasets'
fastq_dir = os.path.join(base_dir, 'SRR_cortex/SRR_cortex_scRNA')
genome_dir = os.path.join(base_dir, 'genome', 'mm10')
output_dir = os.path.join(base_dir, 'SRR_cortex/SRR_cortex_scRNA/processed_h5ad')

os.makedirs(output_dir, exist_ok=True)

# 2. Define Reference Files
fasta_path = os.path.join(genome_dir, 'Mus_musculus.GRCm38.dna.primary_assembly.fa')
gtf_path = os.path.join(genome_dir, 'Mus_musculus.GRCm38.84.gtf')

# 3. Define Index Output Paths
index_path = os.path.join(genome_dir, 'transcriptome.idx')
t2g_path = os.path.join(genome_dir, 'transcripts_to_genes.txt')
t_fasta_path = os.path.join(genome_dir, 'transcriptome.fa')

print("Directories and paths configured.")


# In[ ]:


get_ipython().run_cell_magic('bash', '', '\nFASTQ_DIR="/home/nakagawa/datasets/SRR_cortex/SRR_cortex_scRNA"\nOUT_DIR="/home/nakagawa/datasets/SRR_cortex/SRR_cortex_scRNA/processed_h5ad"\n\nINDEX="/home/nakagawa/datasets/genome/mm10/transcriptome.idx"\nT2G="/home/nakagawa/datasets/genome/mm10/transcripts_to_genes.txt"\n\nmkdir -p "$OUT_DIR"\n\nfor R2 in ${FASTQ_DIR}/*_2.fastq.gz; do\n\n    SAMPLE=$(basename "$R2" _2.fastq.gz)\n\n    echo "---------------------------------------------------"\n    echo "Processing $SAMPLE..."\n\n    SAMPLE_OUT="$OUT_DIR/$SAMPLE"\n\n    R1="$FASTQ_DIR/${SAMPLE}_1.fastq.gz"\n\n    # Skip if R1 missing\n    if [[ ! -f "$R1" ]]; then\n        echo "Missing R1 for $SAMPLE"\n        continue\n    fi\n\n    # Remove old failed outputs if they exist\n    rm -rf "$SAMPLE_OUT"\n\n    echo "Trying 10xv3..."\n\n    kb count \\\n        -i "$INDEX" \\\n        -g "$T2G" \\\n        -x 10xv3 \\\n        -o "$SAMPLE_OUT" \\\n        --h5ad \\\n        "$R1" "$R2"\n\n    STATUS=$?\n\n    # If 10xv3 failed, retry with 10xv2\n    if [[ $STATUS -ne 0 ]]; then\n\n        echo "10xv3 failed for $SAMPLE"\n        echo "Retrying with 10xv2..."\n\n        rm -rf "$SAMPLE_OUT"\n\n        kb count \\\n            -i "$INDEX" \\\n            -g "$T2G" \\\n            -x 10xv2 \\\n            -o "$SAMPLE_OUT" \\\n            --h5ad \\\n            "$R1" "$R2"\n\n        STATUS=$?\n    fi\n\n    if [[ $STATUS -eq 0 ]]; then\n        echo "Finished $SAMPLE successfully"\n    else\n        echo "FAILED: $SAMPLE"\n    fi\n\ndone\n')


# In[ ]:




