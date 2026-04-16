#!/bin/bash
set -euo pipefail

# ==================================================
# CONFIG
# ==================================================
PLINK_BIN="${PLINK_BIN:-plink_install/plink}"
PLINK2_BIN="${PLINK2_BIN:-plink2}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
VCFTOOLS_BIN="${VCFTOOLS_BIN:-vcftools}"
BGZIP_BIN="${BGZIP_BIN:-bgzip}"
TABIX_BIN="${TABIX_BIN:-tabix}"
JAVA_BIN="${JAVA_BIN:-java}"

BEAGLE_JAR="${BEAGLE_JAR:-beagle.jar}"
HAPIBD_JAR="${HAPIBD_JAR:-hap-ibd.jar}"
GENETIC_MAP_DIR="${GENETIC_MAP_DIR:-genetic_map_GRCh37}"

THREADS="${THREADS:-4}"
OUTDIR="${OUTDIR:-./behar_vcf_ibd_output}"

MAF="${MAF:-0.01}"
MAX_MISSING="${MAX_MISSING:-0.05}"
REL_CUTOFF="${REL_CUTOFF:-0.25}"

PRUNE_WINDOW="${PRUNE_WINDOW:-200}"
PRUNE_STEP="${PRUNE_STEP:-25}"
PRUNE_R2="${PRUNE_R2:-0.2}"

HAPIBD_MIN_OUTPUT="${HAPIBD_MIN_OUTPUT:-3}"

# ==================================================
# HELPERS
# ==================================================
log_step() {
    echo
    echo "============================================================"
    echo "$1"
    echo "============================================================"
}

require_file() {
    local f="$1"
    if [ ! -f "$f" ]; then
        echo "ERROR: required file not found: $f"
        exit 1
    fi
}

require_dir() {
    local d="$1"
    if [ ! -d "$d" ]; then
        echo "ERROR: required directory not found: $d"
        exit 1
    fi
}

require_cmd_or_file() {
    local x="$1"
    if command -v "$x" >/dev/null 2>&1; then
        return 0
    fi
    if [ -f "$x" ]; then
        return 0
    fi
    echo "ERROR: command or file not found: $x"
    exit 1
}

count_vcf_stats() {
    local vcf="$1"
    local label="$2"
    local nsamples
    local nvars

    nsamples=$("$BCFTOOLS_BIN" query -l "$vcf" | wc -l)
    nvars=$("$BCFTOOLS_BIN" view -H "$vcf" | wc -l)

    echo "$label: $nsamples samples, $nvars variants"
    printf "%s\t%s\t%s\n" "$label" "$nsamples" "$nvars" >> "$SUMMARY_FILE"
}

count_plink_stats() {
    local prefix="$1"
    local label="$2"
    local nsamples
    local nsnps

    nsamples=$(wc -l < "${prefix}.fam")
    nsnps=$(wc -l < "${prefix}.bim")

    echo "$label: $nsamples samples, $nsnps SNPs"
    printf "%s\t%s\t%s\n" "$label" "$nsamples" "$nsnps" >> "$SUMMARY_FILE"
}

detect_autosomes_from_header() {
    local vcf="$1"
    "$BCFTOOLS_BIN" view -h "$vcf" \
        | grep '^##contig=<ID=' \
        | sed -E 's/^##contig=<ID=([^,>]+).*/\1/' \
        | grep -E '^(chr)?([1-9]|1[0-9]|2[0-2])$' \
        | paste -sd, -
}

detect_chr_style() {
    local vcf="$1"
    if "$BCFTOOLS_BIN" view -h "$vcf" | grep -q '^##contig=<ID=chr1[,>]' ; then
        echo "chr"
    else
        echo "plain"
    fi
}

# ==================================================
# USAGE
# ==================================================
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input.vcf>"
    exit 1
fi

INPUT_VCF="$1"
require_file "$INPUT_VCF"

require_cmd_or_file "$PLINK_BIN"
require_cmd_or_file "$PLINK2_BIN"
require_cmd_or_file "$BCFTOOLS_BIN"
require_cmd_or_file "$VCFTOOLS_BIN"
require_cmd_or_file "$BGZIP_BIN"
require_cmd_or_file "$TABIX_BIN"
require_cmd_or_file "$JAVA_BIN"

require_file "$BEAGLE_JAR"
require_file "$HAPIBD_JAR"
require_dir "$GENETIC_MAP_DIR"

# ==================================================
# OUTPUT PATHS
# ==================================================
mkdir -p "$OUTDIR"

REPORT_DIR="$OUTDIR/reports"
TMP_DIR="$OUTDIR/tmp"
DIAG_DIR="$OUTDIR/diagnostics"
LOG_DIR="$OUTDIR/logs"
IBD_DIR="$OUTDIR/ibd"
CHR_VCF_DIR="$IBD_DIR/per_chr_vcf"
PHASED_DIR="$IBD_DIR/phased_beagle"
HAPIBD_DIR="$IBD_DIR/hapibd_results"

mkdir -p "$REPORT_DIR" "$TMP_DIR" "$DIAG_DIR" "$LOG_DIR" "$IBD_DIR" "$CHR_VCF_DIR" "$PHASED_DIR" "$HAPIBD_DIR"

INPUT_COPY="$OUTDIR/input.vcf"
INPUT_BGZ="$OUTDIR/input.vcf.gz"
NORM_VCF="$OUTDIR/input.norm.vcf.gz"
AUTO_VCF="$OUTDIR/input.autosomes.vcf.gz"
STRICT_VCF="$OUTDIR/input.autosomes.biallelic.vcf.gz"
FILTERED_VCF="$OUTDIR/input.filtered.vcf.gz"

PLINK_PREFIX="$OUTDIR/input_plink"
MIND_PREFIX="$OUTDIR/input_mind"
NOAMB_PREFIX="$OUTDIR/input_noamb"

PRUNED_PREFIX="$OUTDIR/pruned"
PRUNED_DATASET="$OUTDIR/pruned_dataset"
REL_PRUNED_PREFIX="$OUTDIR/final_unrelated_pruned"
FINAL_PREFIX="$OUTDIR/final_unrelated_full"

REMOVE_IDS_FILE="$REPORT_DIR/removed_related_samples.txt"
KEEP_IDS_FILE="$REPORT_DIR/keep_unrelated_samples.txt"

IBD_ALL="$IBD_DIR/ibd_all.ibd"
HBD_ALL="$IBD_DIR/hbd_all.hbd"

SUMMARY_FILE="$REPORT_DIR/summary.tsv"
printf "Stage\tSamples\tVariants\n" > "$SUMMARY_FILE"

# ==================================================
# STEP 0: BGZIP + INDEX INPUT
# ==================================================
log_step "STEP 0: BGZIP + INDEX INPUT VCF"

cp "$INPUT_VCF" "$INPUT_COPY"
"$BGZIP_BIN" -f "$INPUT_COPY"
"$TABIX_BIN" -f -p vcf "$INPUT_BGZ"

count_vcf_stats "$INPUT_BGZ" "Input VCF"

# ==================================================
# STEP 1: NORMALIZE MULTIALLELIC SITES
# ==================================================
log_step "STEP 1: NORMALIZE MULTIALLELIC SITES"

"$BCFTOOLS_BIN" norm -m -any "$INPUT_BGZ" -Oz -o "$NORM_VCF"
"$TABIX_BIN" -f -p vcf "$NORM_VCF"

count_vcf_stats "$NORM_VCF" "Normalized VCF"

# ==================================================
# STEP 2: DETECT AUTOSOMAL CONTIGS
# ==================================================
log_step "STEP 2: DETECT AUTOSOMAL CONTIGS"

AUTOSOME_LIST=$(detect_autosomes_from_header "$NORM_VCF")

if [ -z "$AUTOSOME_LIST" ]; then
    echo "ERROR: no autosomal contigs detected in $NORM_VCF"
    echo "Contigs found in header:"
    "$BCFTOOLS_BIN" view -h "$NORM_VCF" | grep '^##contig=<ID=' || true
    exit 1
fi

echo "Detected autosomal contigs:"
echo "$AUTOSOME_LIST" | tr ',' '\n' | tee "$DIAG_DIR/autosomes_detected.txt"

CHR_STYLE=$(detect_chr_style "$NORM_VCF")
echo "Chromosome naming style: $CHR_STYLE" | tee "$DIAG_DIR/chromosome_style.txt"

# ==================================================
# STEP 3: KEEP AUTOSOMES ONLY
# ==================================================
log_step "STEP 3: KEEP AUTOSOMES ONLY"

"$BCFTOOLS_BIN" view \
    -r "$AUTOSOME_LIST" \
    -Oz \
    -o "$AUTO_VCF" \
    "$NORM_VCF"

"$TABIX_BIN" -f -p vcf "$AUTO_VCF"

count_vcf_stats "$AUTO_VCF" "Autosomal VCF"

# ==================================================
# STEP 4: KEEP AUTOSOMAL BIALLELIC SNPS
# ==================================================
log_step "STEP 4: KEEP AUTOSOMAL BIALLELIC SNPS"

"$BCFTOOLS_BIN" view \
    -m2 -M2 -v snps \
    -Oz \
    -o "$STRICT_VCF" \
    "$AUTO_VCF"

"$TABIX_BIN" -f -p vcf "$STRICT_VCF"

count_vcf_stats "$STRICT_VCF" "Autosomal biallelic SNP VCF"

# ==================================================
# STEP 5: FILTER BY MAF + VARIANT MISSINGNESS
# ==================================================
log_step "STEP 5: FILTER BY MAF AND VARIANT MISSINGNESS"

MIN_NONMISSING=$(awk -v x="$MAX_MISSING" 'BEGIN{printf "%.6f", 1-x}')

"$VCFTOOLS_BIN" \
    --gzvcf "$STRICT_VCF" \
    --maf "$MAF" \
    --max-missing "$MIN_NONMISSING" \
    --recode --stdout \
    | "$BGZIP_BIN" -c > "$FILTERED_VCF"

"$TABIX_BIN" -f -p vcf "$FILTERED_VCF"

count_vcf_stats "$FILTERED_VCF" "Filtered VCF"

# ==================================================
# STEP 6: CONVERT TO PLINK
# ==================================================
log_step "STEP 6: CONVERT TO PLINK"

"$PLINK_BIN" \
    --vcf "$FILTERED_VCF" \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --make-bed \
    --out "$PLINK_PREFIX"

count_plink_stats "$PLINK_PREFIX" "PLINK after VCF conversion"

# ==================================================
# STEP 7: REMOVE SAMPLES WITH HIGH MISSINGNESS
# ==================================================
log_step "STEP 7: REMOVE SAMPLES WITH HIGH MISSINGNESS"

"$PLINK_BIN" \
    --bfile "$PLINK_PREFIX" \
    --mind "$MAX_MISSING" \
    --allow-extra-chr \
    --make-bed \
    --out "$MIND_PREFIX"

count_plink_stats "$MIND_PREFIX" "PLINK after sample missingness filter"

# ==================================================
# STEP 8: REMOVE AMBIGUOUS SNPs
# ==================================================
log_step "STEP 8: REMOVE AMBIGUOUS SNPs"

AMB_FILE="$REPORT_DIR/ambiguous_snps_removed.txt"
KEEP_NONAMB_FILE="$REPORT_DIR/non_ambiguous_snps_kept.txt"

: > "$AMB_FILE"
: > "$KEEP_NONAMB_FILE"

awk '
{
    a=toupper($5); b=toupper($6);
    if ((a=="A" && b=="T") || (a=="T" && b=="A") ||
        (a=="C" && b=="G") || (a=="G" && b=="C")) {
        print $2 >> amb
    } else {
        print $2 >> keep
    }
}
' amb="$AMB_FILE" keep="$KEEP_NONAMB_FILE" "${MIND_PREFIX}.bim"

AMB_N=$(wc -l < "$AMB_FILE" 2>/dev/null || echo 0)
KEEP_N=$(wc -l < "$KEEP_NONAMB_FILE" 2>/dev/null || echo 0)

echo "Ambiguous SNPs removed: $AMB_N"
echo "Non-ambiguous SNPs kept: $KEEP_N"

"$PLINK_BIN" \
    --bfile "$MIND_PREFIX" \
    --extract "$KEEP_NONAMB_FILE" \
    --snps-only just-acgt \
    --biallelic-only strict \
    --allow-extra-chr \
    --make-bed \
    --out "$NOAMB_PREFIX"

count_plink_stats "$NOAMB_PREFIX" "PLINK after ambiguous SNP removal"

# ==================================================
# STEP 9: LD-PRUNING FOR RELATEDNESS ONLY
# ==================================================
log_step "STEP 9: LD-PRUNING FOR RELATEDNESS ONLY"

"$PLINK_BIN" \
    --bfile "$NOAMB_PREFIX" \
    --allow-extra-chr \
    --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_R2" \
    --out "$PRUNED_PREFIX"

PRUNE_KEEP=$(wc -l < "${PRUNED_PREFIX}.prune.in")
PRUNE_DROP=$(wc -l < "${PRUNED_PREFIX}.prune.out")

echo "SNPs kept after pruning: $PRUNE_KEEP"
echo "SNPs removed by pruning: $PRUNE_DROP"

"$PLINK_BIN" \
    --bfile "$NOAMB_PREFIX" \
    --extract "${PRUNED_PREFIX}.prune.in" \
    --allow-extra-chr \
    --make-bed \
    --out "$PRUNED_DATASET"

count_plink_stats "$PRUNED_DATASET" "Pruned dataset for relatedness"

# ==================================================
# STEP 10: REMOVE RELATIVES ON PRUNED DATASET
# ==================================================
log_step "STEP 10: REMOVE RELATIVES ON PRUNED DATASET"

"$PLINK_BIN" \
    --bfile "$PRUNED_DATASET" \
    --allow-extra-chr \
    --rel-cutoff "$REL_CUTOFF" \
    --make-bed \
    --out "$REL_PRUNED_PREFIX"

count_plink_stats "$REL_PRUNED_PREFIX" "Final unrelated dataset (pruned)"

awk '{print $1"\t"$2}' "${PRUNED_DATASET}.fam" | sort -u > "$TMP_DIR/all_samples_before_rel.txt"
awk '{print $1"\t"$2}' "${REL_PRUNED_PREFIX}.fam" | sort -u > "$KEEP_IDS_FILE"

comm -23 "$TMP_DIR/all_samples_before_rel.txt" "$KEEP_IDS_FILE" > "$REMOVE_IDS_FILE"

REMOVED_REL=$(wc -l < "$REMOVE_IDS_FILE")
KEPT_REL=$(wc -l < "$KEEP_IDS_FILE")

echo "Related samples removed: $REMOVED_REL"
echo "Samples kept as unrelated: $KEPT_REL"

if [ "$REMOVED_REL" -gt 0 ]; then
    echo "Removed related sample IDs:"
    cat "$REMOVE_IDS_FILE"
fi

# ==================================================
# STEP 11: APPLY SAME SAMPLE REMOVAL TO FULL-DENSITY DATASET
# ==================================================
log_step "STEP 11: APPLY UNRELATED SAMPLE LIST TO FULL-DENSITY DATASET"

"$PLINK_BIN" \
    --bfile "$NOAMB_PREFIX" \
    --keep "$KEEP_IDS_FILE" \
    --allow-extra-chr \
    --make-bed \
    --out "$FINAL_PREFIX"

count_plink_stats "$FINAL_PREFIX" "Final unrelated full-density dataset"

# ==================================================
# STEP 12: EXPORT FINAL FULL-DENSITY DATASET TO PER-CHR VCF.GZ
# ==================================================
log_step "STEP 12: EXPORT FINAL FULL-DENSITY DATASET TO PER-CHROMOSOME VCF.GZ"

for c in $(seq 1 22); do
    echo "Exporting chromosome $c ..."

    if [ "$CHR_STYLE" = "chr" ]; then
        CHR_ARG="chr${c}"
    else
        CHR_ARG="${c}"
    fi

    "$PLINK2_BIN" \
        --bfile "$FINAL_PREFIX" \
        --allow-extra-chr \
        --chr "$CHR_ARG" \
        --export vcf bgz \
        --out "$CHR_VCF_DIR/pre_chr${c}"
done

# ==================================================
# STEP 13: PHASE EACH CHROMOSOME WITH BEAGLE
# ==================================================
log_step "STEP 13: PHASE EACH CHROMOSOME WITH BEAGLE"

for c in $(seq 1 22); do
    CHR_INPUT_VCF="$CHR_VCF_DIR/pre_chr${c}.vcf.gz"
    MAP_FILE="$GENETIC_MAP_DIR/plink.chr${c}.GRCh37.map"
    OUT_PREFIX="$PHASED_DIR/pre_chr${c}"

    require_file "$CHR_INPUT_VCF"
    require_file "$MAP_FILE"

    echo "Phasing chromosome $c ..."

    "$JAVA_BIN" -jar "$BEAGLE_JAR" \
        gt="$CHR_INPUT_VCF" \
        map="$MAP_FILE" \
        out="$OUT_PREFIX" \
        nthreads="$THREADS"
done

# ==================================================
# STEP 14: RUN HAP-IBD
# ==================================================
log_step "STEP 14: RUN HAP-IBD"

for c in $(seq 1 22); do
    PHASED_VCF="$PHASED_DIR/pre_chr${c}.vcf.gz"
    MAP_FILE="$GENETIC_MAP_DIR/plink.chr${c}.GRCh37.map"
    OUT_PREFIX="$HAPIBD_DIR/ibd_chr${c}"

    require_file "$PHASED_VCF"
    require_file "$MAP_FILE"

    echo "Running hap-ibd for chromosome $c ..."

    "$JAVA_BIN" -jar "$HAPIBD_JAR" \
        gt="$PHASED_VCF" \
        map="$MAP_FILE" \
        out="$OUT_PREFIX" \
        min-output="$HAPIBD_MIN_OUTPUT"
done

# ==================================================
# STEP 15: CONCATENATE IBD/HBD RESULTS
# ==================================================
log_step "STEP 15: CONCATENATE IBD/HBD RESULTS"

: > "$IBD_ALL"
: > "$HBD_ALL"

for c in $(seq 1 22); do
    IBD_GZ="$HAPIBD_DIR/ibd_chr${c}.ibd.gz"
    HBD_GZ="$HAPIBD_DIR/ibd_chr${c}.hbd.gz"

    if [ -f "$IBD_GZ" ]; then
        zcat "$IBD_GZ" >> "$IBD_ALL"
    fi

    if [ -f "$HBD_GZ" ]; then
        zcat "$HBD_GZ" >> "$HBD_ALL"
    fi
done

IBD_LINES=$(wc -l < "$IBD_ALL" 2>/dev/null || echo 0)
HBD_LINES=$(wc -l < "$HBD_ALL" 2>/dev/null || echo 0)

echo "Combined IBD segments: $IBD_LINES"
echo "Combined HBD segments: $HBD_LINES"

# ==================================================
# FINAL SUMMARY
# ==================================================
log_step "FINAL OUTPUT SUMMARY"

echo "Main outputs:"
echo "  Full-density unrelated PLINK: ${FINAL_PREFIX}.bed/.bim/.fam"
echo "  Removed relatives:            $REMOVE_IDS_FILE"
echo "  Unrelated keep list:          $KEEP_IDS_FILE"
echo "  Summary table:                $SUMMARY_FILE"
echo "  Per-chromosome VCFs:          $CHR_VCF_DIR"
echo "  Phased VCFs:                  $PHASED_DIR"
echo "  hap-ibd results:              $HAPIBD_DIR"
echo "  Combined IBD file:            $IBD_ALL"
echo "  Combined HBD file:            $HBD_ALL"

echo
echo "DONE"