set -euo pipefail
if [[ $# -lt 7 ]]; then
    echo "Usage: bash $0 input.vcf.gz 2303/idpop.csv GROUP1 GROUP2 PROJ1,PROJ2 gadma_template_or_dir outdir"
    exit 1
fi
VCF="$1"
IDPOP="$2"
GROUP1="$3"
GROUP2="$4"
PROJ="$5"
GADMA_SOURCE="$6"
OUTDIR="$7"
mkdir -p "$OUTDIR"
for x in python3 bcftools tabix bgzip; do
    if ! command -v "$x" >/dev/null 2>&1; then
        echo "ERROR: '$x' not found in PATH"
        exit 1
    fi
done
if ! command -v easySFS >/dev/null 2>&1; then
    echo "ERROR: easySFS not found in PATH"
    echo "Activate the environment where easySFS is installed."
    exit 1
fi
if ! command -v gadma >/dev/null 2>&1; then
    echo "ERROR: gadma not found in PATH"
    echo "Activate the environment where GADMA is installed."
    exit 1
fi
if [[ ! -f "$VCF" ]]; then
    echo "ERROR: VCF not found: $VCF"
    exit 1
fi
if [[ ! -f "$IDPOP" ]]; then
    echo "ERROR: idpop file not found: $IDPOP"
    exit 1
fi
if [[ ! -e "$GADMA_SOURCE" ]]; then
    echo "ERROR: GADMA template file or directory not found: $GADMA_SOURCE"
    exit 1
fi
if [[ ! -f "${VCF}.tbi" && ! -f "${VCF}.csi" ]]; then
    echo "Index not found for VCF. Creating tabix index..."
    tabix -p vcf "$VCF"
fi
echo "[1/8] Preparing sample lists and popfile"
python3 <<PY
import pandas as pd
from pathlib import Path
idpop = pd.read_csv("${IDPOP}")
required = {"IID", "groups"}
missing = required - set(idpop.columns)
if missing:
    raise ValueError(f"idpop.csv must contain columns {required}, missing: {missing}")
g1 = "${GROUP1}"
g2 = "${GROUP2}"
outdir = Path("${OUTDIR}")
sub = idpop[idpop["groups"].isin([g1, g2])].copy()
sub["IID"] = sub["IID"].astype(str)
if (sub["groups"] == g1).sum() == 0:
    raise ValueError(f"No samples found for group: {g1}")
if (sub["groups"] == g2).sum() == 0:
    raise ValueError(f"No samples found for group: {g2}")
sub1 = (
    sub[sub["groups"] == g1][["IID", "groups"]]
    .drop_duplicates()
    .sort_values("IID")
)
sub2 = (
    sub[sub["groups"] == g2][["IID", "groups"]]
    .drop_duplicates()
    .sort_values("IID")
)
sub1["IID"].to_csv(outdir / "group1.samples.txt", index=False, header=False)
sub2["IID"].to_csv(outdir / "group2.samples.txt", index=False, header=False)
merged = pd.concat([sub1, sub2], axis=0)
merged["IID"].to_csv(outdir / "all.samples.txt", index=False, header=False)
with open(outdir / "popfile.txt", "w") as f:
    for _, row in merged.iterrows():
        f.write(f"{row['IID']}\t{row['groups']}\n")
print(f"{g1}: {len(sub1)} samples")
print(f"{g2}: {len(sub2)} samples")
print(f"Popfile written: {outdir / 'popfile.txt'}")
PY
echo "[2/8] Checking sample intersection with VCF"
bcftools query -l "$VCF" | sort > "$OUTDIR/vcf.samples.txt"
sort "$OUTDIR/all.samples.txt" > "$OUTDIR/all.samples.sorted.txt"
comm -12 "$OUTDIR/all.samples.sorted.txt" "$OUTDIR/vcf.samples.txt" > "$OUTDIR/samples.in.vcf.txt"
comm -23 "$OUTDIR/all.samples.sorted.txt" "$OUTDIR/vcf.samples.txt" > "$OUTDIR/samples.missing.txt" || true
N_PRESENT=$(wc -l < "$OUTDIR/samples.in.vcf.txt")
if [[ "$N_PRESENT" -eq 0 ]]; then
    echo "ERROR: none of the selected samples are present in the VCF"
    exit 1
fi
echo "Samples present in VCF: $N_PRESENT"
if [[ -s "$OUTDIR/samples.missing.txt" ]]; then
    echo "WARNING: some samples are missing in VCF:"
    cat "$OUTDIR/samples.missing.txt"
fi
python3 <<PY
from pathlib import Path
outdir = Path("${OUTDIR}")
present = set(x.strip() for x in (outdir / "samples.in.vcf.txt").read_text().splitlines() if x.strip())
lines = []
for line in (outdir / "popfile.txt").read_text().splitlines():
    if not line.strip():
        continue
    s, g = line.split()[:2]
    if s in present:
        lines.append((s, g))
if not lines:
    raise ValueError("No samples remain in popfile after intersecting with VCF.")
with open(outdir / "popfile.present.txt", "w") as f:
    for s, g in lines:
        f.write(f"{s}\t{g}\n")
print(f"Samples retained for easySFS: {len(lines)}")
PY
echo "[3/8] Filtering VCF to selected samples, autosomes, biallelic SNPs"
FIRST_CHROM=$(bcftools query -f '%CHROM\n' "$VCF" | grep -m1 -E '^(chr)?[0-9]+$' || true)
if [[ -z "$FIRST_CHROM" ]]; then
    echo "ERROR: could not detect chromosome naming style in VCF"
    exit 1
fi
if [[ "$FIRST_CHROM" =~ ^chr ]]; then
    REGIONS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
else
    REGIONS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
fi
echo "Detected chromosome style: $FIRST_CHROM"
echo "Using autosomal regions: $REGIONS"
bcftools view \
    -S "$OUTDIR/samples.in.vcf.txt" \
    -r "$REGIONS" \
    -m2 -M2 -v snps \
    -Oz -o "$OUTDIR/selected.2pop.autosomes.snps.vcf.gz" \
    "$VCF"
tabix -f -p vcf "$OUTDIR/selected.2pop.autosomes.snps.vcf.gz"
NVAR=$(bcftools index -n "$OUTDIR/selected.2pop.autosomes.snps.vcf.gz")
NSAMPLES=$(bcftools query -l "$OUTDIR/selected.2pop.autosomes.snps.vcf.gz" | wc -l)
echo "Filtered VCF samples: $NSAMPLES"
echo "Filtered VCF variants: $NVAR"
if [[ "$NVAR" -eq 0 ]]; then
    echo "ERROR: filtered VCF contains 0 variants"
    exit 1
fi
echo "[3.5/8] Thinning VCF for easySFS to reduce memory usage"
if command -v vcftools >/dev/null 2>&1; then
    vcftools \
        --gzvcf "$OUTDIR/selected.2pop.autosomes.snps.vcf.gz" \
        --thin 10000 \
        --recode --stdout | bgzip -c > "$OUTDIR/selected.2pop.autosomes.snps.thin.vcf.gz"
    tabix -f -p vcf "$OUTDIR/selected.2pop.autosomes.snps.thin.vcf.gz"
    EASY_VCF="$OUTDIR/selected.2pop.autosomes.snps.thin.vcf.gz"
else
    echo "WARNING: vcftools not found, using unthinned VCF"
    EASY_VCF="$OUTDIR/selected.2pop.autosomes.snps.vcf.gz"
fi
NVAR_EASY=$(bcftools index -n "$EASY_VCF")
echo "easySFS input VCF: $EASY_VCF"
echo "Variants in easySFS VCF: $NVAR_EASY"
if [[ "$NVAR_EASY" -eq 0 ]]; then
    echo "ERROR: easySFS input VCF contains 0 variants"
    exit 1
fi
echo "[4/8] Running easySFS preview to inspect projection choices"
mkdir -p "$OUTDIR/easySFS_preview"
set +e
easySFS \
    -i "$EASY_VCF" \
    -p "$OUTDIR/popfile.present.txt" \
    -a \
    --preview \
    -o "$OUTDIR/easySFS_preview"
PREVIEW_STATUS=$?
set -e
echo
if [[ "$PREVIEW_STATUS" -ne 0 ]]; then
    echo "WARNING: easySFS preview failed or was killed."
    echo "Continuing with user-supplied projection: $PROJ"
else
    echo "Preview done."
    echo "Check terminal output to confirm projection is sensible."
fi
echo "Requested projection: $PROJ"
echo
echo "[5/8] Building folded SFS with projection $PROJ"
mkdir -p "$OUTDIR/easySFS_out"
easySFS \
    -i "$EASY_VCF" \
    -p "$OUTDIR/popfile.present.txt" \
    -a \
    -f \
    --proj "$PROJ" \
    -o "$OUTDIR/easySFS_out"
echo "[6/8] Locating folded SFS"
DADI_DIR="$OUTDIR/easySFS_out/dadi"
if [[ ! -d "$DADI_DIR" ]]; then
  echo "ERROR: dadi directory not found: $DADI_DIR" >&2
  exit 1
fi
echo "Available SFS files:"
find "$DADI_DIR" -maxdepth 1 -type f -name "*.sfs" | sort
CHOSEN_SFS=$(find "$DADI_DIR" -maxdepth 1 -type f -name "*.sfs" | \
  grep -E "/${GROUP1}-${GROUP2}([.-].*)?\.sfs$|/${GROUP2}-${GROUP1}([.-].*)?\.sfs$" | \
  head -n 1)
if [[ -z "$CHOSEN_SFS" ]]; then
  echo "ERROR: could not find 2-population SFS for ${GROUP1} and ${GROUP2}" >&2
  exit 1
fi
export CHOSEN_SFS
echo "Using SFS file: $CHOSEN_SFS"
echo "[7/8] Preparing one or multiple GADMA runs"
export GADMA_SOURCE OUTDIR GROUP1 GROUP2 PROJ CHOSEN_SFS
mapfile -t GADMA_JOBS < <(python3 <<'PY'
import os
from pathlib import Path
source = Path(os.environ["GADMA_SOURCE"])
outdir = Path(os.environ["OUTDIR"])
g1 = os.environ["GROUP1"]
g2 = os.environ["GROUP2"]
proj = [x.strip() for x in os.environ["PROJ"].split(",") if x.strip()]
fs_file = Path(os.environ["CHOSEN_SFS"])
if len(proj) != 2:
    raise ValueError(f"Projection must contain exactly two comma-separated values, got: {proj}")
if source.is_file():
    templates = [source]
elif source.is_dir():
    candidates = []
    for pattern in ("*.params", "*.param", "*.txt"):
        candidates.extend(source.glob(pattern))
    templates = sorted({p.resolve() for p in candidates})
    if not templates:
        raise FileNotFoundError(
            f"No template files found in {source}. Expected *.params, *.param, or *.txt"
        )
else:
    raise ValueError(f"Unsupported GADMA source: {source}")
generated_dir = outdir / "gadma_params_generated"
runs_dir = outdir / "gadma_runs"
generated_dir.mkdir(parents=True, exist_ok=True)
runs_dir.mkdir(parents=True, exist_ok=True)
updates_common = {
    "Input data": str(fs_file),
    "Population labels": f"[{g1}, {g2}]",
    "Projections": f"[{proj[0]}, {proj[1]}]",
    "Outgroup": "False",
}
manifest_lines = []
for template in templates:
    stem = template.stem
    safe_stem = stem.replace(" ", "_")
    run_outdir = runs_dir / safe_stem
    generated_params = generated_dir / f"{safe_stem}.generated.params"
    text = template.read_text()
    lines = text.splitlines()
    updates = dict(updates_common)
    updates["Output directory"] = str(run_outdir)
    seen = set()
    new_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or ":" not in line:
            new_lines.append(line)
            continue
        key = line.split(":", 1)[0].strip()
        if key in updates:
            new_lines.append(f"{key}: {updates[key]}")
            seen.add(key)
        else:
            new_lines.append(line)
    for key, value in updates.items():
        if key not in seen:
            new_lines.append(f"{key}: {value}")
    generated_params.write_text("\n".join(new_lines) + "\n")
    manifest_lines.append(f"{template}\t{generated_params}\t{run_outdir}")
    print(f"{generated_params}\t{run_outdir}")
(outdir / "gadma_jobs.tsv").write_text(
    "template\tgenerated_params\trun_output_dir\n" + "\n".join(manifest_lines) + "\n"
)
PY
)
if [[ ${#GADMA_JOBS[@]} -eq 0 ]]; then
    echo "ERROR: no GADMA jobs were prepared"
    exit 1
fi
echo "Prepared ${#GADMA_JOBS[@]} GADMA run(s)"
printf '%s\n' "${GADMA_JOBS[@]}"
echo "[8/8] Launching GADMA runs"
SUCCESS_LIST="$OUTDIR/gadma_successful_runs.txt"
FAILED_LIST="$OUTDIR/gadma_failed_runs.txt"
: > "$SUCCESS_LIST"
: > "$FAILED_LIST"
run_gadma() {
    local params_file="$1"
    local status
    set +e
    gadma -p "$params_file"
    status=$?
    set -e
    if [[ "$status" -ne 0 ]]; then
        echo "Primary launch mode failed for $params_file. Trying fallback: gadma $params_file"
        set +e
        gadma "$params_file"
        status=$?
        set -e
    fi
    return "$status"
}
JOB_IDX=0
for job in "${GADMA_JOBS[@]}"; do
    JOB_IDX=$((JOB_IDX + 1))
    IFS=$'\t' read -r GENERATED_PARAMS RUN_OUTDIR <<< "$job"
    echo
    echo "=== GADMA job ${JOB_IDX}/${#GADMA_JOBS[@]} ==="
    echo "Params:  $GENERATED_PARAMS"
    echo "Outdir:  $RUN_OUTDIR"
    if run_gadma "$GENERATED_PARAMS"; then
        echo "$GENERATED_PARAMS" >> "$SUCCESS_LIST"
        echo "Status: SUCCESS"
    else
        echo "$GENERATED_PARAMS" >> "$FAILED_LIST"
        echo "Status: FAILED"
    fi
done
N_SUCCESS=$(wc -l < "$SUCCESS_LIST")
N_FAILED=$(wc -l < "$FAILED_LIST")
echo
echo "Finished. Main outputs:"
echo "  Popfile:                 $OUTDIR/popfile.present.txt"
echo "  Filtered full VCF:       $OUTDIR/selected.2pop.autosomes.snps.vcf.gz"
echo "  Thinned VCF:             $EASY_VCF"
echo "  easySFS preview:         $OUTDIR/easySFS_preview"
echo "  easySFS result:          $OUTDIR/easySFS_out"
echo "  Chosen SFS:              $(cat "$OUTDIR/chosen_sfs.txt")"
echo "  GADMA jobs manifest:     $OUTDIR/gadma_jobs.tsv"
echo "  Generated params dir:    $OUTDIR/gadma_params_generated"
echo "  GADMA runs dir:          $OUTDIR/gadma_runs"
echo "  Successful runs list:    $SUCCESS_LIST"
echo "  Failed runs list:        $FAILED_LIST"
echo "  Successful runs:         $N_SUCCESS"
echo "  Failed runs:             $N_FAILED"
echo
echo "If needed, rerun any failed model manually using its generated params file."
if [[ "$N_SUCCESS" -eq 0 ]]; then
    echo "ERROR: all GADMA runs failed"
    exit 1
fi
