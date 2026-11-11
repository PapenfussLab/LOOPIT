#!/bin/bash
# download_reference_data.sh
# Downloads LOOPIT reference data
# Usage: bash download_reference_data.sh [genome] [outdir]
# Supported genomes: GRCh38, GRCh37

set -euo pipefail

GENOME=${1:-GRCh38}
OUTDIR=${2:-./loopit_data/${GENOME}}
mkdir -p ${OUTDIR}

echo "=============================================="
echo "  LOOPIT reference data setup for ${GENOME}"
echo "  Output: ${OUTDIR}"
echo "=============================================="
echo ""

case ${GENOME} in

    GRCh38|hg38)
        FASTA_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        FAI_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        BLACKLIST_URL="https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"   
        ANNOTATION_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.basic.annotation.gtf.gz"
        REPEAT_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
        FASTA_NAME="Homo_sapiens_assembly38.fasta"
        ANNOTATION_NAME="gencode.v42.primary_assembly.basic.annotation.gtf"
        ;;

    GRCh37|hg19)
        FASTA_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
        FAI_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
        BLACKLIST_URL="https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz"
        ANNOTATION_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.basic.annotation.gtf.gz"
        REPEAT_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz"
        FASTA_NAME="Homo_sapiens_assembly19.fasta"
        ANNOTATION_NAME="gencode.v19.primary_assembly.basic.annotation.gtf"
        ;;

    *)
        echo "ERROR: Unsupported genome ${GENOME}"
        echo "Supported genomes: GRCh38, GRCh37"
        exit 1
        ;;
esac

# ── reference genome ──
if [ ! -f ${OUTDIR}/${FASTA_NAME} ]; then
    echo "Downloading ${GENOME} reference genome (~3GB)..."
    wget --show-progress -O ${OUTDIR}/${FASTA_NAME} ${FASTA_URL}
    echo "  ✓ ${FASTA_NAME}"
else
    echo "  ✓ ${FASTA_NAME} (already exists)"
fi

# ── FAI index ──
if [ ! -f ${OUTDIR}/${FASTA_NAME}.fai ]; then
    echo "Downloading FAI index..."
    wget --show-progress -O ${OUTDIR}/${FASTA_NAME}.fai ${FAI_URL}
    echo "  ✓ ${FASTA_NAME}.fai"
else
    echo "  ✓ ${FASTA_NAME}.fai (already exists)"
fi

# ── blacklist ──
if [ ! -f ${OUTDIR}/encode_blacklist.bed ]; then
    echo "Downloading blacklist..."
    wget --show-progress -O ${OUTDIR}/encode_blacklist.bed.gz ${BLACKLIST_URL}
    gunzip ${OUTDIR}/encode_blacklist.bed.gz
    echo "  ✓ encode_blacklist.bed"
else
    echo "  ✓ encode_blacklist.bed (already exists)"
fi

# ── annotation ──
if [ ! -f ${OUTDIR}/${ANNOTATION_NAME} ]; then
    echo "Downloading gene annotation..."
    wget --show-progress -O ${OUTDIR}/${ANNOTATION_NAME}.gz ${ANNOTATION_URL}
    gunzip ${OUTDIR}/${ANNOTATION_NAME}.gz
    echo "  ✓ ${ANNOTATION_NAME}"
else
    echo "  ✓ ${ANNOTATION_NAME} (already exists)"
fi
# ── repeat mask ──
if [ ! -f ${OUTDIR}/repeat_mask.tsv ]; then
    echo "Downloading UCSC RepeatMasker (~500MB)..."
    wget --show-progress -O ${OUTDIR}/repeat_mask.tsv.gz ${REPEAT_URL}
    gunzip ${OUTDIR}/repeat_mask.tsv.gz

    NCOLS=$(head -1 ${OUTDIR}/repeat_mask.tsv | awk -F'\t' '{print NF}')
    EXPECTED_NCOLS=17
    if [ "$NCOLS" -ne "$EXPECTED_NCOLS" ]; then
        echo "ERROR: repeat_mask.tsv has $NCOLS columns, expected $EXPECTED_NCOLS." >&2
        echo "UCSC rmsk schema may have changed - check ${REPEAT_URL%.txt.gz}.sql before proceeding." >&2
        exit 1
    fi

    sed -i '1i bin\tswScore\tmilliDiv\tmilliDel\tmilliIns\tgenoName\tgenoStart\tgenoEnd\tgenoLeft\tstrand\trepName\trepClass\trepFamily\trepStart\trepEnd\trepLeft\tid' ${OUTDIR}/repeat_mask.tsv
    echo "  ✓ repeat_mask.tsv"
else
    echo "  ✓ repeat_mask.tsv (already exists)"
fi

echo ""
echo "=============================================="
echo "  All reference files ready in ${OUTDIR}"
echo "=============================================="
echo ""
echo "Run LOOPIT with:"
echo ""
echo "  nextflow run main.nf \\"
echo "    --genome        ${GENOME} \\"
echo "    --reference_dir ${OUTDIR} \\"
echo "    --blacklist     ${OUTDIR}/encode_blacklist.bed \\"
echo "    --annotation    ${OUTDIR}/${ANNOTATION_NAME}"
