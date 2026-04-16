#!/bin/bash
set -euo pipefail

# =========================
# CONFIG
# =========================
PLINK_BIN="${PLINK_BIN:-plink_install/plink}"
PLINK2_BIN="${PLINK2_BIN:-plink2}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
VCFTOOLS_BIN="${VCFTOOLS_BIN:-vcftools}"
ADMIXTURE_BIN="${ADMIXTURE_BIN:-admixture}"
BGZIP_BIN="${BGZIP_BIN:-bgzip}"
TABIX_BIN="${TABIX_BIN:-tabix}"
BEAGLE_JAR="${BEAGLE_JAR:-beagle.jar}"
HAPIBD_JAR="${HAPIBD_JAR:-hap-ibd.jar}"
GENETIC_MAP_DIR="${GENETIC_MAP_DIR:-genetic_map_GRCh37}"
JAVA_BIN="${JAVA_BIN:-java}"
THREADS="${THREADS:-4}"

RUN_ADMIXTURE="${RUN_ADMIXTURE:-1}"
RUN_ROH="${RUN_ROH:-1}"
RUN_IBD="${RUN_IBD:-1}"

K_MIN="${K_MIN:-3}"
K_MAX="${K_MAX:-15}"
N_REPS="${N_REPS:-10}"

MAF="${MAF:-0.01}"
MAX_MISSING="${MAX_MISSING:-0.05}"
REL_CUTOFF="${REL_CUTOFF:-0.25}"
PRUNE_WINDOW="${PRUNE_WINDOW:-200}"
PRUNE_STEP="${PRUNE_STEP:-25}"
PRUNE_R2="${PRUNE_R2:-0.2}"

ROH_SNP="${ROH_SNP:-50}"
ROH_KB="${ROH_KB:-500}"
ROH_DENSITY="${ROH_DENSITY:-50}"
ROH_GAP="${ROH_GAP:-1000}"
ROH_WINDOW_SNP="${ROH_WINDOW_SNP:-50}"
ROH_WINDOW_HET="${ROH_WINDOW_HET:-1}"
ROH_WINDOW_MISSING="${ROH_WINDOW_MISSING:-5}"
ROH_WINDOW_THRESHOLD="${ROH_WINDOW_THRESHOLD:-0.05}"

HAPIBD_MIN_OUTPUT="${HAPIBD_MIN_OUTPUT:-3}"
AUTOSOMES="${AUTOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
OUT_DIR_NAME="${OUT_DIR_NAME:-output_0204}"

# =========================
# FUNCTIONS
# =========================
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

extract_cv() {
    local log_file="$1"
    awk '/CV error/ {print $NF}' "$log_file" | tail -n1
}

extract_ll() {
    local log_file="$1"
    awk '
        /Loglikelihood/ {val=$NF}
        /Log-likelihood/ {val=$NF}
        END {
            if (val == "") print "NA";
            else print val;
        }
    ' "$log_file"
}

write_chr_map_strip_chr() {
    local out="$1"
    cat > "$out" <<'MAP'
chr1	1
chr2	2
chr3	3
chr4	4
chr5	5
chr6	6
chr7	7
chr8	8
chr9	9
chr10	10
chr11	11
chr12	12
chr13	13
chr14	14
chr15	15
chr16	16
chr17	17
chr18	18
chr19	19
chr20	20
chr21	21
chr22	22
chrX	X
chrY	Y
chrM	MT
chrMT	MT
MAP
}

activate_plink2_env() {
    if command -v "$PLINK2_BIN" >/dev/null 2>&1; then
        echo "plink2 available: $(command -v "$PLINK2_BIN")"
        return 0
    fi

    echo "ERROR: plink2 not found in current environment."
    echo "Activate your conda environment first, then rerun the script."
    echo "Example:"
    echo "  conda activate plink_env"
    echo "  bash $(basename "$0") /path/to/folder"
    exit 1
}

# =========================
# USAGE
# =========================
if [ $# -lt 1 ]; then
    echo "Usage: $0 <vcf_folder>"
    exit 1
fi

VCF_DIR="$1"
require_dir "$VCF_DIR"

ANNOT_FILE="$VCF_DIR/idpop.csv"
require_file "$ANNOT_FILE"

require_cmd_or_file "$PLINK_BIN"
require_cmd_or_file "$BCFTOOLS_BIN"
require_cmd_or_file "$VCFTOOLS_BIN"
require_cmd_or_file "$BGZIP_BIN"
require_cmd_or_file "$TABIX_BIN"

mapfile -t files < <(find "$VCF_DIR" -maxdepth 1 -type f -name '*.vcf.gz' | LC_ALL=C sort)

if [ "${#files[@]}" -ne 2 ]; then
    echo "ERROR: expected exactly 2 files matching *.vcf.gz in $VCF_DIR, found ${#files[@]}"
    printf 'Found files:\n'
    printf '  %s\n' "${files[@]:-<none>}"
    exit 1
fi

VCF_A="${files[0]}"
VCF_B="${files[1]}"
BASE_A=$(basename "$VCF_A" .vcf.gz)
BASE_B=$(basename "$VCF_B" .vcf.gz)

# =========================
# OUTPUT DIRS
# =========================
OUT_DIR="$VCF_DIR/$OUT_DIR_NAME"
REPORT_DIR="$OUT_DIR/reports"
LOG_DIR="$OUT_DIR/logs"
CV_DIR="$OUT_DIR/cv"
BEST_DIR="$OUT_DIR/best"
TMP_DIR="$OUT_DIR/tmp"
DIAG_DIR="$OUT_DIR/diagnostics"
ANALYSIS_DIR="$OUT_DIR/roh_ibd_analysis"
ROH_DIR="$ANALYSIS_DIR/roh"
IBD_DIR="$ANALYSIS_DIR/ibd"
CHR_VCF_DIR="$IBD_DIR/per_chr_vcf"
PHASED_DIR="$IBD_DIR/phased_beagle"
HAPIBD_DIR="$IBD_DIR/hapibd_results"

mkdir -p "$OUT_DIR" "$REPORT_DIR" "$LOG_DIR" "$CV_DIR" "$BEST_DIR" "$TMP_DIR" "$DIAG_DIR"
mkdir -p "$ROH_DIR" "$IBD_DIR" "$CHR_VCF_DIR" "$PHASED_DIR" "$HAPIBD_DIR"
cp "$ANNOT_FILE" "$OUT_DIR/"

SUMMARY_FILE="$REPORT_DIR/summary.tsv"
NORMALIZATION_STATS="$REPORT_DIR/normalization_stats.tsv"
SHARED_POSITIONS="$REPORT_DIR/shared_positions.tsv"
ALL_COMMON_IDS="$REPORT_DIR/shared_sample_annotation_order.tsv"
KEEP_IDS_FILE="$REPORT_DIR/keep_unrelated_samples.txt"
REMOVE_IDS_FILE="$REPORT_DIR/removed_related_samples.txt"

printf "Stage\tSamples\tVariants\n" > "$SUMMARY_FILE"
printf "Sample\tTotal\tSplit\tJoined\tRealigned\tSkipped\n" > "$NORMALIZATION_STATS"

# =========================
# STEP 0: FIND INPUT FILES
# =========================
log_step "STEP 0: FIND INPUT FILES"

echo "VCF A: $VCF_A"
echo "VCF B: $VCF_B"
echo "Annotation file: $ANNOT_FILE"
echo "Output directory: $OUT_DIR"

# =========================
# STEP 1: PREPARE ANNOTATED IDS
# =========================
log_step "STEP 1: PREPARE ANNOTATED SAMPLE LIST"

ANNOT_IDS="$OUT_DIR/annotated_ids.txt"

awk -F',' 'NR>1 {gsub(/\r/,"",$1); if($1!="") print $1}' "$ANNOT_FILE" \
    | sed 's/\r$//' \
    | awk '{gsub(/^[ \t]+|[ \t]+$/, "", $0); if($0!="") print $0}' \
    | LC_ALL=C sort -u \
    > "$ANNOT_IDS"

N_ANNOT=$(wc -l < "$ANNOT_IDS")
echo "Annotated IDs in idpop.csv: $N_ANNOT"

# =========================
# STEP 2: SUBSET + NORMALIZE BOTH VCFs
# =========================
log_step "STEP 2: SUBSET EACH VCF TO ANNOTATED SAMPLES AND NORMALIZE"

for INPUT_VCF in "$VCF_A" "$VCF_B"; do
    BASENAME=$(basename "$INPUT_VCF" .vcf.gz)
    SAMPLE_LIST="$OUT_DIR/${BASENAME}.samples.txt"
    KEEP_LIST="$OUT_DIR/${BASENAME}.keep.txt"
    SUBSET_VCF="$OUT_DIR/${BASENAME}.subset.vcf.gz"
    NORM_VCF="$OUT_DIR/${BASENAME}.norm.vcf.gz"

    "$BCFTOOLS_BIN" query -l "$INPUT_VCF" \
        | sed 's/\r$//' \
        | awk '{gsub(/^[ \t]+|[ \t]+$/, "", $0); if($0!="") print $0}' \
        | LC_ALL=C sort -u \
        > "$SAMPLE_LIST"

    grep -Fxf "$ANNOT_IDS" "$SAMPLE_LIST" > "$KEEP_LIST" || true

    N_TOTAL=$(wc -l < "$SAMPLE_LIST")
    N_KEEP=$(wc -l < "$KEEP_LIST")

    echo "$BASENAME: $N_KEEP / $N_TOTAL samples matched annotation"

    if [ "$N_KEEP" -eq 0 ]; then
        echo "ERROR: no annotated samples found in $INPUT_VCF"
        exit 1
    fi

    "$BCFTOOLS_BIN" view -S "$KEEP_LIST" -Oz -o "$SUBSET_VCF" "$INPUT_VCF"
    "$BCFTOOLS_BIN" index -f "$SUBSET_VCF"

    "$BCFTOOLS_BIN" norm -m -any "$SUBSET_VCF" -Ou \
        | "$BCFTOOLS_BIN" norm -d none -Oz -o "$NORM_VCF"

    "$BCFTOOLS_BIN" index -f "$NORM_VCF"

    TOTAL=$("$BCFTOOLS_BIN" view -H "$SUBSET_VCF" | wc -l)
    SPLIT=$("$BCFTOOLS_BIN" norm -m -any "$SUBSET_VCF" -Ou | "$BCFTOOLS_BIN" view -H | wc -l)
    JOINED="$SPLIT"
    REALIGNED=0
    SKIPPED=0

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$BASENAME" "$TOTAL" "$SPLIT" "$JOINED" "$REALIGNED" "$SKIPPED" \
        >> "$NORMALIZATION_STATS"

    count_vcf_stats "$NORM_VCF" "$BASENAME normalized"
done

NORM_A="$OUT_DIR/${BASE_A}.norm.vcf.gz"
NORM_B="$OUT_DIR/${BASE_B}.norm.vcf.gz"

# =========================
# STEP 3A: DIAGNOSTIC OF CHROMOSOME NAMES
# =========================
log_step "STEP 3A: DIAGNOSTIC OF CHROMOSOME NAMES"

CHROMS_A="$DIAG_DIR/${BASE_A}.chroms.txt"
CHROMS_B="$DIAG_DIR/${BASE_B}.chroms.txt"

"$BCFTOOLS_BIN" query -f '%CHROM\n' "$NORM_A" | LC_ALL=C sort -u > "$CHROMS_A"
"$BCFTOOLS_BIN" query -f '%CHROM\n' "$NORM_B" | LC_ALL=C sort -u > "$CHROMS_B"

echo "Unique chromosome names in A:"
cat "$CHROMS_A"
echo
echo "Unique chromosome names in B:"
cat "$CHROMS_B"

# =========================
# STEP 3B: HARMONIZE CHROMOSOME NAMES
# =========================
log_step "STEP 3B: HARMONIZE CHROMOSOME NAMES"

RENAMED_A="$NORM_A"
RENAMED_B="$NORM_B"
CHR_MODE="unchanged"

if grep -q '^chr' "$CHROMS_A" && ! grep -q '^chr' "$CHROMS_B"; then
    CHR_MAP="$TMP_DIR/chr_rename_A_strip_chr.txt"
    write_chr_map_strip_chr "$CHR_MAP"
    RENAMED_A="$TMP_DIR/${BASE_A}.norm.renamed.vcf.gz"
    "$BCFTOOLS_BIN" annotate --rename-chrs "$CHR_MAP" -Oz -o "$RENAMED_A" "$NORM_A"
    "$BCFTOOLS_BIN" index -f "$RENAMED_A"
    CHR_MODE="A: chr-prefixed -> plain"
elif ! grep -q '^chr' "$CHROMS_A" && grep -q '^chr' "$CHROMS_B"; then
    CHR_MAP="$TMP_DIR/chr_rename_B_strip_chr.txt"
    write_chr_map_strip_chr "$CHR_MAP"
    RENAMED_B="$TMP_DIR/${BASE_B}.norm.renamed.vcf.gz"
    "$BCFTOOLS_BIN" annotate --rename-chrs "$CHR_MAP" -Oz -o "$RENAMED_B" "$NORM_B"
    "$BCFTOOLS_BIN" index -f "$RENAMED_B"
    CHR_MODE="B: chr-prefixed -> plain"
fi

echo "Chromosome harmonization mode: $CHR_MODE" | tee "$DIAG_DIR/chrom_harmonization.txt"

NORM_A="$RENAMED_A"
NORM_B="$RENAMED_B"

"$BCFTOOLS_BIN" query -f '%CHROM\n' "$NORM_A" | LC_ALL=C sort -u > "$DIAG_DIR/${BASE_A}.chroms.after.txt"
"$BCFTOOLS_BIN" query -f '%CHROM\n' "$NORM_B" | LC_ALL=C sort -u > "$DIAG_DIR/${BASE_B}.chroms.after.txt"

# =========================
# STEP 3C: RESTRICT TO AUTOSOMES
# =========================
log_step "STEP 3C: RESTRICT BOTH DATASETS TO AUTOSOMES"

AUTO_A="$OUT_DIR/${BASE_A}.autosomes.vcf.gz"
AUTO_B="$OUT_DIR/${BASE_B}.autosomes.vcf.gz"

"$BCFTOOLS_BIN" view -r "$AUTOSOMES" -Oz -o "$AUTO_A" "$NORM_A"
"$BCFTOOLS_BIN" index -f "$AUTO_A"

"$BCFTOOLS_BIN" view -r "$AUTOSOMES" -Oz -o "$AUTO_B" "$NORM_B"
"$BCFTOOLS_BIN" index -f "$AUTO_B"

NORM_A="$AUTO_A"
NORM_B="$AUTO_B"

count_vcf_stats "$NORM_A" "$BASE_A autosomes"
count_vcf_stats "$NORM_B" "$BASE_B autosomes"

# =========================
# STEP 3D: DIAGNOSTIC OF SHARED POSITIONS
# =========================
log_step "STEP 3D: DIAGNOSTIC OF SHARED POSITIONS"

A_POS="$DIAG_DIR/${BASE_A}.positions.tsv"
B_POS="$DIAG_DIR/${BASE_B}.positions.tsv"
A_FULL="$DIAG_DIR/${BASE_A}.sites_full.tsv"
B_FULL="$DIAG_DIR/${BASE_B}.sites_full.tsv"

"$BCFTOOLS_BIN" query -f '%CHROM\t%POS\n' "$NORM_A" | LC_ALL=C sort -u > "$A_POS"
"$BCFTOOLS_BIN" query -f '%CHROM\t%POS\n' "$NORM_B" | LC_ALL=C sort -u > "$B_POS"
"$BCFTOOLS_BIN" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$NORM_A" | LC_ALL=C sort -u > "$A_FULL"
"$BCFTOOLS_BIN" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$NORM_B" | LC_ALL=C sort -u > "$B_FULL"

SHARED_POS_COUNT=$(LC_ALL=C comm -12 "$A_POS" "$B_POS" | tee "$SHARED_POSITIONS" | wc -l)
SHARED_FULL_COUNT=$(LC_ALL=C comm -12 "$A_FULL" "$B_FULL" | wc -l)

echo "Shared positions (CHROM POS): $SHARED_POS_COUNT" | tee "$DIAG_DIR/overlap_summary.txt"
echo "Shared exact sites (CHROM POS REF ALT): $SHARED_FULL_COUNT" | tee -a "$DIAG_DIR/overlap_summary.txt"

if [ "$SHARED_POS_COUNT" -eq 0 ]; then
    echo "ERROR: no shared autosomal positions found between the two VCFs after chromosome harmonization"
    exit 1
fi

# =========================
# STEP 4: SUBSET BOTH VCFs TO SHARED POSITIONS
# =========================
log_step "STEP 4: SUBSET BOTH VCFs TO SHARED POSITIONS"

A_COMMON="$OUT_DIR/${BASE_A}.sharedpos.vcf.gz"
B_COMMON="$OUT_DIR/${BASE_B}.sharedpos.vcf.gz"

"$BCFTOOLS_BIN" view -T "$SHARED_POSITIONS" -Oz -o "$A_COMMON" "$NORM_A"
"$BCFTOOLS_BIN" index -f "$A_COMMON"

"$BCFTOOLS_BIN" view -T "$SHARED_POSITIONS" -Oz -o "$B_COMMON" "$NORM_B"
"$BCFTOOLS_BIN" index -f "$B_COMMON"

count_vcf_stats "$A_COMMON" "$BASE_A shared-positions"
count_vcf_stats "$B_COMMON" "$BASE_B shared-positions"

# =========================
# STEP 5: MERGE SHARED-POSITION VCFs
# =========================
log_step "STEP 5: MERGE SHARED-POSITION VCFs"

MERGED="$OUT_DIR/merged_shared_positions.vcf.gz"
SORTED_MERGED="$OUT_DIR/merged_shared_positions.sorted.vcf.gz"

"$BCFTOOLS_BIN" merge \
    --force-samples \
    -m none \
    -Oz \
    -o "$MERGED" \
    "$A_COMMON" "$B_COMMON"

"$BCFTOOLS_BIN" index -f "$MERGED"
count_vcf_stats "$MERGED" "Merged shared-position VCF"

"$BCFTOOLS_BIN" sort "$MERGED" -Oz -o "$SORTED_MERGED"
"$BCFTOOLS_BIN" index -f "$SORTED_MERGED"
count_vcf_stats "$SORTED_MERGED" "Sorted merged shared-position VCF"

# =========================
# STEP 6: QC FILTERING
# =========================
log_step "STEP 6: FILTER TO AUTOSOMAL BIALLELIC SNPs WITH MAF AND MISSINGNESS"

STRICT_VCF="$OUT_DIR/merged_shared_positions.strict.vcf.gz"
FILTERED="$OUT_DIR/merged_shared_positions.strict.filtered.vcf.gz"
MIN_NONMISSING=$(awk -v x="$MAX_MISSING" 'BEGIN{printf "%.6f", 1-x}')

"$BCFTOOLS_BIN" view \
    -r "$AUTOSOMES" \
    -m2 -M2 -v snps \
    -Oz -o "$STRICT_VCF" \
    "$SORTED_MERGED"

"$BCFTOOLS_BIN" index -f "$STRICT_VCF"
count_vcf_stats "$STRICT_VCF" "Autosomal biallelic SNP VCF"

"$VCFTOOLS_BIN" \
    --gzvcf "$STRICT_VCF" \
    --maf "$MAF" \
    --max-missing "$MIN_NONMISSING" \
    --recode --stdout \
    | "$BGZIP_BIN" -c > "$FILTERED"

"$BCFTOOLS_BIN" index -f "$FILTERED"
count_vcf_stats "$FILTERED" "Filtered VCF"

# =========================
# STEP 7: CONVERT TO PLINK
# =========================
log_step "STEP 7: CONVERT FILTERED VCF TO PLINK"

FILTERED_PLINK="$OUT_DIR/merged_input"

"$PLINK_BIN" \
    --vcf "$FILTERED" \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --make-bed \
    --out "$FILTERED_PLINK"

count_plink_stats "$FILTERED_PLINK" "PLINK after VCF conversion"

# =========================
# STEP 8: REMOVE SAMPLES WITH HIGH MISSINGNESS
# =========================
log_step "STEP 8: REMOVE SAMPLES WITH HIGH MISSINGNESS"

MIND_PREFIX="$OUT_DIR/merged_input_mind"

"$PLINK_BIN" \
    --bfile "$FILTERED_PLINK" \
    --mind "$MAX_MISSING" \
    --allow-extra-chr \
    --make-bed \
    --out "$MIND_PREFIX"

count_plink_stats "$MIND_PREFIX" "PLINK after sample missingness filter"

# =========================
# STEP 9: REMOVE AMBIGUOUS SNPs
# =========================
log_step "STEP 9: REMOVE AMBIGUOUS SNPs"

NOAMB_PREFIX="$OUT_DIR/merged_input_noamb"
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

# =========================
# STEP 10: LD-PRUNING FOR RELATEDNESS ONLY
# =========================
log_step "STEP 10: LD-PRUNING FOR RELATEDNESS ONLY"

PRUNED_PREFIX="$OUT_DIR/pruned_for_relatedness"
PRUNED_REL_DATASET="$OUT_DIR/pruned_dataset_for_relatedness"

"$PLINK_BIN" \
    --bfile "$NOAMB_PREFIX" \
    --allow-extra-chr \
    --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_R2" \
    --out "$PRUNED_PREFIX"

PRUNE_KEEP=$(wc -l < "${PRUNED_PREFIX}.prune.in")
PRUNE_DROP=$(wc -l < "${PRUNED_PREFIX}.prune.out")

echo "SNPs kept after pruning for relatedness: $PRUNE_KEEP"
echo "SNPs removed by pruning for relatedness: $PRUNE_DROP"

"$PLINK_BIN" \
    --bfile "$NOAMB_PREFIX" \
    --extract "${PRUNED_PREFIX}.prune.in" \
    --allow-extra-chr \
    --make-bed \
    --out "$PRUNED_REL_DATASET"

count_plink_stats "$PRUNED_REL_DATASET" "Pruned dataset for relatedness"

# =========================
# STEP 11: REMOVE RELATIVES ON PRUNED DATASET
# =========================
log_step "STEP 11: REMOVE RELATIVES ON PRUNED DATASET"

REL_PRUNED_PREFIX="$OUT_DIR/final_unrelated_pruned_for_relatedness"

"$PLINK_BIN" \
    --bfile "$PRUNED_REL_DATASET" \
    --allow-extra-chr \
    --rel-cutoff "$REL_CUTOFF" \
    --make-bed \
    --out "$REL_PRUNED_PREFIX"

count_plink_stats "$REL_PRUNED_PREFIX" "Final unrelated dataset (pruned for relatedness)"

awk '{print $1"\t"$2}' "${PRUNED_REL_DATASET}.fam" | LC_ALL=C sort -u > "$TMP_DIR/all_samples_before_rel.txt"
awk '{print $1"\t"$2}' "${REL_PRUNED_PREFIX}.fam" | LC_ALL=C sort -u > "$KEEP_IDS_FILE"

LC_ALL=C comm -23 "$TMP_DIR/all_samples_before_rel.txt" "$KEEP_IDS_FILE" > "$REMOVE_IDS_FILE"

REMOVED_REL=$(wc -l < "$REMOVE_IDS_FILE")
KEPT_REL=$(wc -l < "$KEEP_IDS_FILE")

echo "Related samples removed: $REMOVED_REL"
echo "Samples kept as unrelated: $KEPT_REL"

if [ "$REMOVED_REL" -gt 0 ]; then
    echo "Removed related sample IDs:"
    cat "$REMOVE_IDS_FILE"
fi

# =========================
# STEP 12: APPLY SAME SAMPLE REMOVAL TO FULL-DENSITY DATASET
# =========================
log_step "STEP 12: CREATE FINAL FULL-DENSITY UNRELATED DATASET"

FINAL_FULL_PREFIX="$OUT_DIR/final_unrelated_full"

"$PLINK_BIN" \
    --bfile "$NOAMB_PREFIX" \
    --keep "$KEEP_IDS_FILE" \
    --allow-extra-chr \
    --make-bed \
    --out "$FINAL_FULL_PREFIX"

count_plink_stats "$FINAL_FULL_PREFIX" "Final unrelated full-density dataset"

# =========================
# STEP 13: SAMPLE ORDER TABLE
# =========================
log_step "STEP 13: PREPARE SAMPLE ORDER TABLE"

awk -F',' 'NR==1 {next} {gsub(/\r/,"",$1); print $0}' "$ANNOT_FILE" > "$TMP_DIR/annotation_rows_noheader.csv"
awk '{print $1}' "${FINAL_FULL_PREFIX}.fam" > "$TMP_DIR/final_ids.txt"

{
    echo -e "IID\tannotation_row"
    while read -r iid; do
        row=$(awk -F',' -v id="$iid" '$1==id {print $0; exit}' "$TMP_DIR/annotation_rows_noheader.csv")
        echo -e "${iid}\t${row:-NA}"
    done < "$TMP_DIR/final_ids.txt"
} > "$ALL_COMMON_IDS"

echo "Sample order table written to: $ALL_COMMON_IDS"

# =========================
# STEP 14: RUN ROH ON FULL-DENSITY UNRELATED DATASET
# =========================
if [ "$RUN_ROH" = "1" ]; then
    log_step "STEP 14: RUN ROH ANALYSIS"

    ROH_PREFIX="$ROH_DIR/final_unrelated_roh"
    ROH_FILE="${ROH_PREFIX}.hom"
    ROH_SUMMARY_FILE="$ROH_DIR/roh_summary.txt"

    "$PLINK_BIN" \
        --bfile "$FINAL_FULL_PREFIX" \
        --allow-extra-chr \
        --homozyg \
        --homozyg-snp "$ROH_SNP" \
        --homozyg-kb "$ROH_KB" \
        --homozyg-density "$ROH_DENSITY" \
        --homozyg-gap "$ROH_GAP" \
        --homozyg-window-snp "$ROH_WINDOW_SNP" \
        --homozyg-window-het "$ROH_WINDOW_HET" \
        --homozyg-window-missing "$ROH_WINDOW_MISSING" \
        --homozyg-window-threshold "$ROH_WINDOW_THRESHOLD" \
        --out "$ROH_PREFIX"

    if [ -f "$ROH_FILE" ]; then
        ROH_LINES=$(tail -n +2 "$ROH_FILE" | wc -l)
    else
        ROH_LINES=0
    fi

    NSAMPLES=$(wc -l < "${FINAL_FULL_PREFIX}.fam")

    {
        echo "Samples: $NSAMPLES"
        echo "ROH segments: $ROH_LINES"
        echo "ROH file: $ROH_FILE"
    } | tee "$ROH_SUMMARY_FILE"
fi

# =========================
# STEP 15: RUN IBD ON FULL-DENSITY UNRELATED DATASET
# =========================
if [ "$RUN_IBD" = "1" ]; then
    log_step "STEP 15: CHECK IBD REQUIREMENTS"

    require_cmd_or_file "$PLINK2_BIN"
    require_cmd_or_file "$JAVA_BIN"
    require_file "$BEAGLE_JAR"
    require_file "$HAPIBD_JAR"
    require_dir "$GENETIC_MAP_DIR"
    activate_plink2_env

    log_step "STEP 16: EXPORT FULL-DENSITY DATASET TO PER-CHROMOSOME VCF.GZ"

    for c in $(seq 1 22); do
        echo "Exporting chromosome $c ..."
        "$PLINK2_BIN" \
            --bfile "$FINAL_FULL_PREFIX" \
            --allow-extra-chr \
            --chr "$c" \
            --export vcf bgz \
            --out "$CHR_VCF_DIR/pre_chr${c}"
    done

    log_step "STEP 17: PHASE EACH CHROMOSOME WITH BEAGLE"

    for c in $(seq 1 22); do
        INPUT_VCF="$CHR_VCF_DIR/pre_chr${c}.vcf.gz"
        MAP_FILE="$GENETIC_MAP_DIR/plink.chr${c}.GRCh37.map"
        OUT_PREFIX="$PHASED_DIR/pre_chr${c}"

        require_file "$INPUT_VCF"
        require_file "$MAP_FILE"

        echo "Phasing chromosome $c ..."
        "$JAVA_BIN" -jar "$BEAGLE_JAR" \
            gt="$INPUT_VCF" \
            map="$MAP_FILE" \
            out="$OUT_PREFIX" \
            nthreads="$THREADS"
    done

    log_step "STEP 18: RUN HAP-IBD"

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

    log_step "STEP 19: CONCATENATE IBD/HBD"

    IBD_ALL="$IBD_DIR/ibd_all.ibd"
    HBD_ALL="$IBD_DIR/hbd_all.hbd"

    : > "$IBD_ALL"
    : > "$HBD_ALL"

    for c in $(seq 1 22); do
        IBD_FILE="$HAPIBD_DIR/ibd_chr${c}.ibd.gz"
        HBD_FILE="$HAPIBD_DIR/ibd_chr${c}.hbd.gz"

        if [ -f "$IBD_FILE" ]; then
            zcat "$IBD_FILE" >> "$IBD_ALL"
        fi

        if [ -f "$HBD_FILE" ]; then
            zcat "$HBD_FILE" >> "$HBD_ALL"
        fi
    done

    IBD_LINES=$(wc -l < "$IBD_ALL" 2>/dev/null || echo 0)
    HBD_LINES=$(wc -l < "$HBD_ALL" 2>/dev/null || echo 0)

    echo "Combined IBD segments: $IBD_LINES"
    echo "Combined HBD segments: $HBD_LINES"
fi

# =========================
# STEP 20: LD-PRUNE FINAL FULL-DENSITY DATASET FOR ADMIXTURE
# =========================
if [ "$RUN_ADMIXTURE" = "1" ]; then
    log_step "STEP 20: LD-PRUNE FINAL FULL-DENSITY DATASET FOR ADMIXTURE"

    require_cmd_or_file "$ADMIXTURE_BIN"

    ADMIX_PRUNE_PREFIX="$OUT_DIR/pruned_for_admixture"
    FINAL_ADMIX_PREFIX="$OUT_DIR/final_unrelated_pruned"

    "$PLINK_BIN" \
        --bfile "$FINAL_FULL_PREFIX" \
        --allow-extra-chr \
        --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_R2" \
        --out "$ADMIX_PRUNE_PREFIX"

    ADMIX_PRUNE_KEEP=$(wc -l < "${ADMIX_PRUNE_PREFIX}.prune.in")
    ADMIX_PRUNE_DROP=$(wc -l < "${ADMIX_PRUNE_PREFIX}.prune.out")

    echo "SNPs kept after pruning for ADMIXTURE: $ADMIX_PRUNE_KEEP"
    echo "SNPs removed by pruning for ADMIXTURE: $ADMIX_PRUNE_DROP"

    "$PLINK_BIN" \
        --bfile "$FINAL_FULL_PREFIX" \
        --extract "${ADMIX_PRUNE_PREFIX}.prune.in" \
        --allow-extra-chr \
        --make-bed \
        --out "$FINAL_ADMIX_PREFIX"

    count_plink_stats "$FINAL_ADMIX_PREFIX" "Final unrelated LD-pruned dataset for ADMIXTURE"

    log_step "STEP 21: RUN ADMIXTURE FOR MULTIPLE K AND MULTIPLE SEEDS"

    RESULTS_TSV="$CV_DIR/admixture_cv_ll_summary.tsv"
    printf "K\tseed\tCV\tLoglikelihood\tQ_file\tP_file\tlog_file\n" > "$RESULTS_TSV"

    ADMIX_WORKDIR="$OUT_DIR/admixture_runs"
    mkdir -p "$ADMIX_WORKDIR"

    FINAL_BASE=$(basename "$FINAL_ADMIX_PREFIX")

    for K in $(seq "$K_MIN" "$K_MAX"); do
        for SEED in $(seq 1 "$N_REPS"); do
            RUN_PREFIX="$LOG_DIR/K${K}_seed${SEED}"
            LOG_FILE="${RUN_PREFIX}.log"

            echo "Running ADMIXTURE: K=$K seed=$SEED"

            rm -f "./${FINAL_BASE}.${K}.Q" "./${FINAL_BASE}.${K}.P" 2>/dev/null || true

            "$ADMIXTURE_BIN" \
                --cv \
                -j"$THREADS" \
                -s "$SEED" \
                "${FINAL_ADMIX_PREFIX}.bed" "$K" \
                | tee "$LOG_FILE"

            Q_SRC="./${FINAL_BASE}.${K}.Q"
            P_SRC="./${FINAL_BASE}.${K}.P"
            Q_DST="${ADMIX_WORKDIR}/K${K}_seed${SEED}.Q"
            P_DST="${ADMIX_WORKDIR}/K${K}_seed${SEED}.P"

            require_file "$Q_SRC"
            require_file "$P_SRC"

            mv "$Q_SRC" "$Q_DST"
            mv "$P_SRC" "$P_DST"

            CV=$(extract_cv "$LOG_FILE")
            LL=$(extract_ll "$LOG_FILE")

            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$K" "$SEED" "${CV:-NA}" "${LL:-NA}" "$Q_DST" "$P_DST" "$LOG_FILE" \
                >> "$RESULTS_TSV"
        done
    done

    echo "ADMIXTURE summary written to: $RESULTS_TSV"

    log_step "STEP 22: SELECT BEST RUN PER K"

    BEST_TSV="$BEST_DIR/best_runs_per_k.tsv"
    printf "K\tseed\tCV\tLoglikelihood\tQ_file\tP_file\tlog_file\n" > "$BEST_TSV"

    awk 'NR==1 {next}
    {
        k=$1
        seed=$2
        cv=$3
        ll=$4
        line=$0

        if (cv == "NA") next
        if (ll == "NA") ll = -1e99

        if (!(k in best_cv) || cv < best_cv[k] || (cv == best_cv[k] && ll > best_ll[k])) {
            best_cv[k] = cv
            best_ll[k] = ll
            best_line[k] = line
        }
    }
    END {
        for (k in best_line) print best_line[k]
    }' "$RESULTS_TSV" | sort -n >> "$BEST_TSV"

    log_step "STEP 23: SUMMARIZE CV ERROR BY K"

    MEAN_CV_TSV="$BEST_DIR/mean_cv_by_k.tsv"
    printf "K\tmean_CV\tmin_CV\tmax_CV\tn_runs\n" > "$MEAN_CV_TSV"

    awk 'NR==1 {next}
    {
        k=$1
        cv=$3
        if (cv == "NA") next
        sum[k] += cv
        n[k] += 1
        if (!(k in min) || cv < min[k]) min[k] = cv
        if (!(k in max) || cv > max[k]) max[k] = cv
    }
    END {
        for (k in sum) {
            printf "%s\t%.8f\t%.8f\t%.8f\t%d\n", k, sum[k]/n[k], min[k], max[k], n[k]
        }
    }' "$RESULTS_TSV" | sort -n >> "$MEAN_CV_TSV"
fi

# =========================
# FINAL SUMMARY
# =========================
log_step "FINAL OUTPUT SUMMARY"

echo "Main outputs:"
echo "  Full-density unrelated PLINK:    ${FINAL_FULL_PREFIX}.bed/.bim/.fam"
echo "  Keep unrelated IDs:              $KEEP_IDS_FILE"
echo "  Removed relatives:               $REMOVE_IDS_FILE"
echo "  Pipeline summary:                $SUMMARY_FILE"
echo "  Normalization stats:             $NORMALIZATION_STATS"
echo "  Shared positions:                $SHARED_POSITIONS"
echo "  Chrom diagnostics:               $DIAG_DIR"
echo "  Sample order table:              $ALL_COMMON_IDS"

if [ "$RUN_ROH" = "1" ]; then
    echo
    echo "ROH outputs:"
    echo "  ${ROH_PREFIX}.hom"
    echo "  ${ROH_PREFIX}.hom.indiv"
    echo "  ${ROH_PREFIX}.hom.summary"
fi

if [ "$RUN_IBD" = "1" ]; then
    echo
    echo "IBD outputs:"
    echo "  Per-chromosome VCFs:             $CHR_VCF_DIR"
    echo "  Phased VCFs:                     $PHASED_DIR"
    echo "  hap-ibd results:                 $HAPIBD_DIR"
    echo "  Combined IBD file:               $IBD_ALL"
    echo "  Combined HBD file:               $HBD_ALL"
fi

if [ "$RUN_ADMIXTURE" = "1" ]; then
    echo
    echo "ADMIXTURE outputs:"
    echo "  Final pruned PLINK:              ${FINAL_ADMIX_PREFIX}.bed/.bim/.fam"
    echo "  Full run summary:                $RESULTS_TSV"
    echo "  Best run per K:                  $BEST_TSV"
    echo "  Mean CV by K:                    $MEAN_CV_TSV"
fi

echo
echo "DONE"
