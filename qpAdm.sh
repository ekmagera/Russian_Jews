#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# qpAdm pipeline with target-specific automatic model search
# Direct qpAdm, no qpfstats
#
# Input:
#   - jews.bcf
#   - 2303/idpop.csv   (columns: IID, groups)
#   - HO reference in EIGENSTRAT format
#
# Main logic:
#   1. Filter BCF to autosomal biallelic SNPs
#   2. Convert to PLINK
#   3. Basic QC
#   4. Remove strand-ambiguous SNPs
#   5. Remove relatives
#   6. LD prune
#   7. Convert own data to EIGENSTRAT
#   8. Rewrite labels from idpop.csv
#   9. Remove Unknown / No_Tag individuals
#  10. Subset HO to only needed populations
#  11. Merge own data with HO subset
#  12. Generate target-specific qpAdm source models
#  13. Run direct qpAdm for all models
#  14. Parse summary and best models
# ============================================================

# -----------------------------
# paths
# -----------------------------
BCF="${BCF:-jews.bcf}"
IDPOP="${IDPOP:-2303/idpop.csv}"

ADMIXTOOLS_DIR="${ADMIXTOOLS_DIR:-AdmixTools-master/bin}"
REICHLAB_DIR="${REICHLAB_DIR:-${ADMIXTOOLS_DIR}/reichlab}"
REF_PREFIX="${REF_PREFIX:-${REICHLAB_DIR}/v54.1.p1_HO_public}"

PLINK_BIN="${PLINK_BIN:-plink_install/plink}"
PLINK2_BIN="${PLINK2_BIN:-plink2}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"

CONVERTF_BIN="${CONVERTF_BIN:-${ADMIXTOOLS_DIR}/convertf}"
MERGEIT_BIN="${MERGEIT_BIN:-${ADMIXTOOLS_DIR}/mergeit}"
QPADM_BIN="${QPADM_BIN:-${ADMIXTOOLS_DIR}/qpAdm}"

OUTDIR="${OUTDIR:-qpadm_jews_run}"

# перед mkdir -p
rm -rf "${OUTDIR}/qpadm" "${OUTDIR}/lists" "${OUTDIR}/par" "${OUTDIR}/stats"
mkdir -p "${OUTDIR}"/{tmp,logs,par,lists,qpadm,merged,eigen,plink,stats}

# -----------------------------
# targets
# -----------------------------
TARGETS=(
  "Ashkenazi_Jews"
  "Bukharan_Jews"
  "Georgian_Jews"
  "Mountain_Jews"
  "Kurdistan_Jews"
)

# -----------------------------
# qpAdm model sizes
# only 2-way and 3-way
# -----------------------------
MODEL_SIZES=(2 3)

# -----------------------------
# target-specific source candidates
# Lebanese.HO must be present in every model
# -----------------------------
# format:
#   TARGET<TAB>comma,separated,sources
#
# Ashkenazi: Lebanese + Italy + Russia
# Bukharan: Lebanese + Iranian + Uzbek + Russia
# Georgian: Lebanese + Georgian + Russia
# Mountain: Lebanese + Iranian + Azerbaijani + Russia
# Kurdistan: Lebanese + Kurd + Russia
# -----------------------------
TARGET_SOURCE_CONFIG_FILE="${OUTDIR}/lists/target_source_candidates.tsv"

# -----------------------------
# helpers
# -----------------------------
log() {
  echo "[$(date '+%F %T')] $*"
}

require_file() {
  local f="$1"
  [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
}

# -----------------------------
# sanity checks
# -----------------------------
require_file "${BCF}"
require_file "${IDPOP}"
require_file "${REF_PREFIX}.geno"
require_file "${REF_PREFIX}.snp"
require_file "${REF_PREFIX}.ind"

command -v "${PLINK_BIN}" >/dev/null 2>&1 || { echo "ERROR: ${PLINK_BIN} not found" >&2; exit 1; }
command -v "${PLINK2_BIN}" >/dev/null 2>&1 || { echo "ERROR: ${PLINK2_BIN} not found" >&2; exit 1; }
command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1 || { echo "ERROR: ${BCFTOOLS_BIN} not found" >&2; exit 1; }

[[ -x "${CONVERTF_BIN}" ]] || { echo "ERROR: ${CONVERTF_BIN} not executable" >&2; exit 1; }
[[ -x "${MERGEIT_BIN}" ]] || { echo "ERROR: ${MERGEIT_BIN} not executable" >&2; exit 1; }
[[ -x "${QPADM_BIN}" ]] || { echo "ERROR: ${QPADM_BIN} not executable" >&2; exit 1; }

# -----------------------------
# step 1. extract autosomes and keep only biallelic SNPs
# -----------------------------
log "STEP 1: Extract autosomes and keep only biallelic SNPs"

AUTOSOME_BCF="${OUTDIR}/tmp/jews.autosomes.biallelic.bcf"
CONTIGS_TXT="${OUTDIR}/tmp/contigs.txt"

"${BCFTOOLS_BIN}" view -h "${BCF}" | grep '^##contig' > "${CONTIGS_TXT}"

if grep -q 'ID=chr1' "${CONTIGS_TXT}"; then
  REGIONS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
else
  REGIONS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
fi

"${BCFTOOLS_BIN}" view \
  -r "${REGIONS}" \
  -m2 -M2 \
  -v snps \
  -Ob -o "${AUTOSOME_BCF}" \
  "${BCF}" \
  > "${OUTDIR}/logs/01a_bcftools_autosomes.log" 2>&1

"${BCFTOOLS_BIN}" index -f "${AUTOSOME_BCF}" \
  >> "${OUTDIR}/logs/01a_bcftools_autosomes.log" 2>&1

# -----------------------------
# step 1b. BCF -> PLINK
# -----------------------------
log "STEP 1b: Convert filtered BCF to PLINK"

"${PLINK2_BIN}" \
  --bcf "${AUTOSOME_BCF}" \
  --double-id \
  --set-missing-var-ids @:#:\$r:\$a \
  --make-bed \
  --out "${OUTDIR}/plink/raw" \
  > "${OUTDIR}/logs/01b_bcf_to_plink.log" 2>&1

# -----------------------------
# step 2. Basic QC
# -----------------------------
log "STEP 2: Basic QC"

"${PLINK_BIN}" \
  --bfile "${OUTDIR}/plink/raw" \
  --geno 0.05 \
  --mind 0.05 \
  --maf 0.01 \
  --snps-only just-acgt \
  --make-bed \
  --out "${OUTDIR}/plink/qc1" \
  > "${OUTDIR}/logs/02_qc.log" 2>&1

# -----------------------------
# step 2b. Remove strand-ambiguous SNPs
# -----------------------------
log "STEP 2b: Remove strand-ambiguous SNPs"

awk '{
  a=$5; b=$6;
  if ((a=="A" && b=="T") || (a=="T" && b=="A") || (a=="C" && b=="G") || (a=="G" && b=="C")) print $2
}' "${OUTDIR}/plink/qc1.bim" > "${OUTDIR}/tmp/ambiguous_snps.txt"

"${PLINK_BIN}" \
  --bfile "${OUTDIR}/plink/qc1" \
  --exclude "${OUTDIR}/tmp/ambiguous_snps.txt" \
  --make-bed \
  --out "${OUTDIR}/plink/qc2" \
  > "${OUTDIR}/logs/02b_remove_ambiguous.log" 2>&1

# -----------------------------
# step 3. Remove relatives
# -----------------------------
log "STEP 3: Remove relatives"

"${PLINK2_BIN}" \
  --bfile "${OUTDIR}/plink/qc2" \
  --king-cutoff 0.0884 \
  --make-bed \
  --out "${OUTDIR}/plink/unrelated" \
  > "${OUTDIR}/logs/03_remove_relatives.log" 2>&1

# -----------------------------
# step 4. LD pruning
# -----------------------------
log "STEP 4: LD pruning"

"${PLINK_BIN}" \
  --bfile "${OUTDIR}/plink/unrelated" \
  --indep-pairwise 200 25 0.4 \
  --out "${OUTDIR}/plink/prune" \
  > "${OUTDIR}/logs/04_prune.log" 2>&1

"${PLINK_BIN}" \
  --bfile "${OUTDIR}/plink/unrelated" \
  --extract "${OUTDIR}/plink/prune.prune.in" \
  --make-bed \
  --out "${OUTDIR}/plink/unrelated_pruned" \
  > "${OUTDIR}/logs/04b_pruned_dataset.log" 2>&1

# -----------------------------
# step 5. Convert own dataset to EIGENSTRAT
# -----------------------------
log "STEP 5: Prepare EIGENSTRAT conversion"

cat > "${OUTDIR}/par/convertf_own.par" <<EOF
genotypename:    ${OUTDIR}/plink/unrelated_pruned.bed
snpname:         ${OUTDIR}/plink/unrelated_pruned.bim
indivname:       ${OUTDIR}/plink/unrelated_pruned.fam
outputformat:    EIGENSTRAT
genotypeoutname: ${OUTDIR}/eigen/jews.geno
snpoutname:      ${OUTDIR}/eigen/jews.snp
indivoutname:    ${OUTDIR}/eigen/jews.ind
familynames:     NO
EOF

"${CONVERTF_BIN}" -p "${OUTDIR}/par/convertf_own.par" \
  > "${OUTDIR}/logs/05_convertf.log" 2>&1

# -----------------------------
# step 5b. Rewrite jews.ind with population labels
# -----------------------------
log "STEP 5b: Rewrite jews.ind with population labels from idpop.csv"

python3 <<'PY'
import csv
from pathlib import Path

outdir = Path("qpadm_jews_run")
idpop = Path("2303/idpop.csv")

mapping = {}
with open(idpop, newline="", encoding="utf-8-sig") as f:
    r = csv.DictReader(f)
    if "IID" not in r.fieldnames or "groups" not in r.fieldnames:
        raise SystemExit("idpop.csv must contain IID and groups columns")
    for row in r:
        iid = row["IID"].strip()
        grp = row["groups"].strip()
        if iid:
            mapping[iid] = grp

ind_in = outdir / "eigen" / "jews.ind"
ind_out = outdir / "eigen" / "jews_labeled.ind"

n_total = 0
n_mapped = 0

with open(ind_in, encoding="utf-8") as fin, open(ind_out, "w", encoding="utf-8") as fout:
    for line in fin:
        if not line.strip():
            continue
        parts = line.strip().split()
        iid = parts[0]
        sex = parts[1] if len(parts) > 1 else "U"
        grp = mapping.get(iid, "Unknown")
        if grp != "Unknown":
            n_mapped += 1
        n_total += 1
        fout.write(f"{iid} {sex} {grp}\n")

print(f"Total individuals in jews.ind: {n_total}")
print(f"Mapped via idpop.csv: {n_mapped}")
print(f"Wrote: {ind_out}")
PY

mv "${OUTDIR}/eigen/jews_labeled.ind" "${OUTDIR}/eigen/jews.ind"

# -----------------------------
# step 5c. Remove Unknown / No_Tag individuals
# -----------------------------
log "STEP 5c: Remove Unknown / No_Tag individuals"

KEEP_IDS="${OUTDIR}/tmp/keep_ids.txt"

awk '$3!="Unknown" && $3!="No_Tag" {print $1}' "${OUTDIR}/eigen/jews.ind" > "${KEEP_IDS}"

"${PLINK_BIN}" \
  --bfile "${OUTDIR}/plink/unrelated_pruned" \
  --keep <(awk '{print $1, $1}' "${KEEP_IDS}") \
  --make-bed \
  --out "${OUTDIR}/plink/unrelated_pruned_clean" \
  > "${OUTDIR}/logs/05c_filter_samples.log" 2>&1

cat > "${OUTDIR}/par/convertf_clean.par" <<EOF
genotypename:    ${OUTDIR}/plink/unrelated_pruned_clean.bed
snpname:         ${OUTDIR}/plink/unrelated_pruned_clean.bim
indivname:       ${OUTDIR}/plink/unrelated_pruned_clean.fam
outputformat:    EIGENSTRAT
genotypeoutname: ${OUTDIR}/eigen/jews.geno
snpoutname:      ${OUTDIR}/eigen/jews.snp
indivoutname:    ${OUTDIR}/eigen/jews.ind
familynames:     NO
EOF

"${CONVERTF_BIN}" -p "${OUTDIR}/par/convertf_clean.par" \
  > "${OUTDIR}/logs/05c_convertf_clean.log" 2>&1

python3 <<'PY'
import csv
from pathlib import Path

outdir = Path("qpadm_jews_run")
idpop = Path("2303/idpop.csv")

mapping = {}
with open(idpop, newline="", encoding="utf-8-sig") as f:
    r = csv.DictReader(f)
    for row in r:
        mapping[row["IID"].strip()] = row["groups"].strip()

ind_in = outdir / "eigen" / "jews.ind"
ind_out = outdir / "eigen" / "jews_labeled.ind"

with open(ind_in, encoding="utf-8") as fin, open(ind_out, "w", encoding="utf-8") as fout:
    for line in fin:
        parts = line.strip().split()
        iid = parts[0]
        sex = parts[1]
        grp = mapping.get(iid, "Unknown")
        fout.write(f"{iid} {sex} {grp}\n")
PY

mv "${OUTDIR}/eigen/jews_labeled.ind" "${OUTDIR}/eigen/jews.ind"

# -----------------------------
# step 6. Write target-specific source and right-pop candidate table
# -----------------------------
log "STEP 6: Write target-specific source and right-pop candidate table"

TARGET_MODEL_CONFIG_FILE="${OUTDIR}/lists/target_model_config.tsv"

cat > "${TARGET_MODEL_CONFIG_FILE}" <<EOF
target	sources	rights
Ashkenazi_Jews	Lebanese.HO,Italian_South.HO,Russian.HO	Mbuti.HO,Ami.HO,Basque.HO,Biaka.HO,Chukchi.HO,Eskimo_Naukan.HO,Han.HO,Ju_hoan_North.HO,Karitiana.HO,Papuan.HO,She.HO,Ulchi.HO,Yoruba.HO
Bukharan_Jews	Lebanese.HO,Iranian.HO,Uzbek.HO,Tajik.HO,Russian.HO	Mbuti.HO,Ami.HO,Basque.HO,Biaka.HO,Chukchi.HO,Eskimo_Naukan.HO,Han.HO,Ju_hoan_North.HO,Karitiana.HO,Papuan.HO,She.HO,Ulchi.HO,Yoruba.HO
Georgian_Jews	Lebanese.HO,Iranian.HO,Georgian.HO,Russian.HO	Mbuti.HO,Ami.HO,Basque.HO,Biaka.HO,Chukchi.HO,Eskimo_Naukan.HO,Han.HO,Ju_hoan_North.HO,Karitiana.HO,Papuan.HO,She.HO,Ulchi.HO,Yoruba.HO
Mountain_Jews	Lebanese.HO,Lezgin.HO,Azeri.HO,Russian.HO	Mbuti.HO,Ami.HO,Basque.HO,Biaka.HO,Chukchi.HO,Eskimo_Naukan.HO,Han.HO,Ju_hoan_North.HO,Karitiana.HO,Papuan.HO,She.HO,Ulchi.HO,Yoruba.HO
Kurdistan_Jews	Lebanese.HO,Kurd.HO,Georgian.HO,Russian.HO	Mbuti.HO,Ami.HO,Basque.HO,Biaka.HO,Chukchi.HO,Eskimo_Naukan.HO,Han.HO,Ju_hoan_North.HO,Karitiana.HO,Papuan.HO,She.HO,Ulchi.HO,Yoruba.HO
EOF
# -----------------------------
# step 7. Subset HO to all needed populations
# -----------------------------
log "STEP 7: Subset HO to needed populations"

HO_SUB_PREFIX="${OUTDIR}/merged/ho_subset"

python3 <<'PY'
from pathlib import Path
import csv

outdir = Path("qpadm_jews_run")
config_file = outdir / "lists" / "target_model_config.tsv"
out_file = outdir / "lists" / "ho_needed_pops.txt"

needed = set()

with open(config_file, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        for field in ("sources", "rights"):
            vals = [x.strip() for x in row[field].split(",") if x.strip()]
            needed.update(vals)

with open(out_file, "w", encoding="utf-8") as f:
    for pop in sorted(needed):
        f.write(pop + "\n")

print(f"Wrote {out_file}")
print(f"Total HO populations needed: {len(needed)}")
PY

cat > "${OUTDIR}/par/ho_subset_convertf.par" <<EOF
genotypename:    ${REF_PREFIX}.geno
snpname:         ${REF_PREFIX}.snp
indivname:       ${REF_PREFIX}.ind
outputformat:    EIGENSTRAT
genotypeoutname: ${HO_SUB_PREFIX}.geno
snpoutname:      ${HO_SUB_PREFIX}.snp
indivoutname:    ${HO_SUB_PREFIX}.ind
poplistname:     ${OUTDIR}/lists/ho_needed_pops.txt
familynames:     NO
EOF

"${CONVERTF_BIN}" -p "${OUTDIR}/par/ho_subset_convertf.par" \
  > "${OUTDIR}/logs/07a_ho_subset_convertf.log" 2>&1
# -----------------------------
# step 7b. Merge with HO subset
# -----------------------------
log "STEP 7b: Merge with HO subset"

cat > "${OUTDIR}/par/mergeit.par" <<EOF
geno1: ${HO_SUB_PREFIX}.geno
snp1:  ${HO_SUB_PREFIX}.snp
ind1:  ${HO_SUB_PREFIX}.ind
geno2: ${OUTDIR}/eigen/jews.geno
snp2:  ${OUTDIR}/eigen/jews.snp
ind2:  ${OUTDIR}/eigen/jews.ind
genooutfilename: ${OUTDIR}/merged/merged.geno
snpoutfilename:  ${OUTDIR}/merged/merged.snp
indoutfilename:  ${OUTDIR}/merged/merged.ind
outputformat: EIGENSTRAT
docheck: YES
hashcheck: YES
EOF

"${MERGEIT_BIN}" -p "${OUTDIR}/par/mergeit.par" \
  > "${OUTDIR}/logs/07_mergeit.log" 2>&1

# -----------------------------
# step 7c. Save merged population counts
# -----------------------------
log "STEP 7c: Save merged population counts"

awk '{print $3}' "${OUTDIR}/merged/merged.ind" | sort | uniq -c | sort -nr \
  > "${OUTDIR}/stats/merged_pop_counts.txt"

# -----------------------------
# step 8. Generate target-specific source combinations
# with target-specific right populations
# -----------------------------
log "STEP 8: Generate target-specific source combinations"

python3 <<'PY'
from itertools import combinations
from pathlib import Path
import csv

outdir = Path("qpadm_jews_run")
config_file = outdir / "lists" / "target_model_config.tsv"
out_file = outdir / "lists" / "source_combinations.tsv"

model_sizes = [2, 3]

rows = []
combo_counter = 1

with open(config_file, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        target = row["target"].strip()
        candidates = [x.strip() for x in row["sources"].split(",") if x.strip()]
        rights = [x.strip() for x in row["rights"].split(",") if x.strip()]

        if not candidates:
            raise RuntimeError(f"{target}: empty source list")
        if not rights:
            raise RuntimeError(f"{target}: empty right list")

        # target must not appear among sources or rights
        if target in candidates:
            raise RuntimeError(f"{target}: target appears among sources")
        if target in rights:
            raise RuntimeError(f"{target}: target appears among rights")

        # avoid left-right overlap
        overlap = sorted(set(candidates) & set(rights))
        if overlap:
            raise RuntimeError(f"{target}: overlap between sources and rights: {','.join(overlap)}")

        for k in model_sizes:
            if k > len(candidates):
                continue
            for combo in combinations(candidates, k):
                combo_id = f"C{combo_counter:03d}"
                rows.append({
                    "target": target,
                    "combo_id": combo_id,
                    "n_sources": str(k),
                    "sources": ",".join(combo),
                    "rights": ",".join(rights),
                    "n_rights": str(len(rights)),
                })
                combo_counter += 1

with open(out_file, "w", encoding="utf-8", newline="") as f:
    fieldnames = ["target", "combo_id", "n_sources", "sources", "n_rights", "rights"]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {out_file}")
print(f"Total target-specific models: {len(rows)}")
PY
# -----------------------------
# step 9. Run qpAdm over all target-specific combinations
# with target-specific right populations
# -----------------------------
log "STEP 9: Run qpAdm over all target-specific combinations"

while IFS=$'\t' read -r TARGET COMBO_ID N_SOURCES SOURCE_LIST N_RIGHTS RIGHT_LIST; do
  [[ "${TARGET}" == "target" ]] && continue

  IFS=',' read -r -a CUR_SOURCES <<< "${SOURCE_LIST}"
  IFS=',' read -r -a CUR_RIGHTS <<< "${RIGHT_LIST}"

  MODEL_TAG="${TARGET}.${COMBO_ID}"
  LEFT_FILE="${OUTDIR}/lists/${MODEL_TAG}.left.txt"
  RIGHT_FILE="${OUTDIR}/lists/${MODEL_TAG}.right.txt"
  PAR_FILE="${OUTDIR}/par/${MODEL_TAG}.qpadm.par"
  OUT_FILE="${OUTDIR}/qpadm/${MODEL_TAG}.qpadm.txt"

  log "Running qpAdm: ${MODEL_TAG} | left=${SOURCE_LIST} | right=${RIGHT_LIST}"

  {
    echo "${TARGET}"
    for s in "${CUR_SOURCES[@]}"; do
      echo "$s"
    done
  } > "${LEFT_FILE}"

  {
    for r in "${CUR_RIGHTS[@]}"; do
      echo "$r"
    done
  } > "${RIGHT_FILE}"

  cat > "${PAR_FILE}" <<EOF
genotypename: ${OUTDIR}/merged/merged.geno
snpname:      ${OUTDIR}/merged/merged.snp
indivname:    ${OUTDIR}/merged/merged.ind
popleft:      ${LEFT_FILE}
popright:     ${RIGHT_FILE}
details:      YES
allsnps:      YES
inbreed:      NO
EOF

  "${QPADM_BIN}" -p "${PAR_FILE}" > "${OUT_FILE}" 2> "${OUT_FILE}.stderr" || true
done < "${OUTDIR}/lists/source_combinations.tsv"

# -----------------------------
# step 10. Parse qpAdm summary
# -----------------------------

log "STEP 10: Parse qpAdm summary"

python3 <<'PY'
import re
from pathlib import Path
import csv

qpadm_dir = Path("qpadm_jews_run/qpadm")
combo_file = Path("qpadm_jews_run/lists/source_combinations.tsv")
summary_out = Path("qpadm_jews_run/stats/qpadm_summary.tsv")
best_out = Path("qpadm_jews_run/stats/qpadm_best_models.tsv")
best_fallback_out = Path("qpadm_jews_run/stats/qpadm_best_models_with_fallback.tsv")

SE_WARNING_THRESHOLD = 0.5
FLOAT_RE = r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?"

def extract_floats(text):
    if not text:
        return []
    return [float(x) for x in re.findall(FLOAT_RE, text)]

def fmt_float_list(vals):
    if not vals:
        return ""
    return ",".join(f"{x:.6g}" for x in vals)

def read_lines(text):
    return text.splitlines()

def parse_tail_prob(text):
    lines = read_lines(text)

    # 1) best source: summ line
    # summ: Kurdistan_Jews    3      0.506272     0.341 ...
    for line in lines:
        if line.lstrip().startswith("summ:"):
            nums = extract_floats(line)
            if nums:
                # first number is n_sources, second is tail_prob
                if len(nums) >= 2:
                    return nums[1]
                # fallback: if format strange, take first float after target
                if len(nums) >= 1:
                    return nums[0]

    # 2) codimension / main model line
    # f4rank: 2 dof: 11 chisq: 10.270 tail: 0.506271583 ...
    # avoid full rank block with dof 0 tail 1
    for line in lines:
        if "f4rank:" in line and "tail:" in line and "dof:" in line:
            m_dof = re.search(r"dof:\s*(\d+)", line)
            m_tail = re.search(r"tail:\s*(" + FLOAT_RE + r")", line)
            if m_dof and m_tail:
                dof = int(m_dof.group(1))
                if dof > 0:
                    return float(m_tail.group(1))

    # 3) best pat 000
    # best pat:          000         0.506272              -  -
    for line in lines:
        m = re.search(r"best pat:\s*0+\s+(" + FLOAT_RE + r")", line)
        if m:
            return float(m.group(1))

    # 4) fixed pat table, use 000 row if present
    # 000  0  11  10.270  0.506272  0.342  0.694 -0.036 infeasible
    in_fixed_table = False
    for line in lines:
        if re.search(r"fixed\s+pat.*tail\s+prob", line, flags=re.IGNORECASE):
            in_fixed_table = True
            continue
        if in_fixed_table:
            if not line.strip():
                break
            m = re.match(
                r"\s*([01]+)\s+\d+\s+\d+\s+" + FLOAT_RE + r"\s+(" + FLOAT_RE + r")\b",
                line
            )
            if m and set(m.group(1)) == {"0"}:
                return float(m.group(2))

    return None

def parse_coeffs(text, expected_n):
    lines = read_lines(text)
    candidates = []

    for line in lines:
        if re.search(r"best coefficients\s*:", line, flags=re.IGNORECASE):
            nums = extract_floats(line)
            if nums:
                candidates.append(nums)
        elif re.search(r"coeffs\s*:", line, flags=re.IGNORECASE):
            nums = extract_floats(line)
            if nums:
                candidates.append(nums)
        elif re.search(r"totmean\s*:", line, flags=re.IGNORECASE):
            nums = extract_floats(line)
            if nums:
                candidates.append(nums)
        elif re.search(r"boot mean\s*:", line, flags=re.IGNORECASE):
            nums = extract_floats(line)
            if nums:
                candidates.append(nums)
        elif re.search(r"zzjmean\b", line, flags=re.IGNORECASE):
            nums = extract_floats(line)
            if nums:
                candidates.append(nums)

    exact = [v for v in candidates if len(v) == expected_n]
    if exact:
        return exact[0]
    return candidates[0] if candidates else []

def parse_stderr(text, expected_n):
    lines = read_lines(text)
    candidates = []

    for line in lines:
        if re.search(r"std\.?\s*errors?\s*:", line, flags=re.IGNORECASE):
            nums = extract_floats(line)
            if nums:
                candidates.append(nums)
        elif re.search(r"stderr(?:ors?)?\s*:", line, flags=re.IGNORECASE):
            nums = extract_floats(line)
            if nums:
                candidates.append(nums)
        elif re.search(r"standard\s*errors?\s*:", line, flags=re.IGNORECASE):
            nums = extract_floats(line)
            if nums:
                candidates.append(nums)

    exact = [v for v in candidates if len(v) == expected_n]
    if exact:
        return exact[0]
    return candidates[0] if candidates else []

combo_map = {}
with open(combo_file, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        combo_map[(row["target"], row["combo_id"])] = {
            "n_sources": int(row["n_sources"]),
            "sources": row["sources"],
            "n_rights": int(row["n_rights"]),
            "rights": row["rights"],
        }
rows = []

for fp in sorted(qpadm_dir.glob("*.qpadm.txt")):
    name = fp.name
    m = re.match(r"(.+)\.(C\d+)\.qpadm\.txt$", name)
    if not m:
        continue

    target, combo_id = m.group(1), m.group(2)
    meta = combo_map.get((target, combo_id))
    if not meta:
        continue

    expected_n = meta["n_sources"]
    text = fp.read_text(encoding="utf-8", errors="ignore")

    tail_prob = parse_tail_prob(text)
    coeffs = parse_coeffs(text, expected_n)
    std_errors = parse_stderr(text, expected_n)

    has_negative = any(x < 0 for x in coeffs)
    has_gt_one = any(x > 1 for x in coeffs)

    coeff_len_ok = (len(coeffs) == expected_n)
    se_len_ok = (len(std_errors) == expected_n) if std_errors else False

    coeffs_valid = coeff_len_ok and (not has_negative) and (not has_gt_one)
    se_suspicious = (len(std_errors) > 0) and any(x > SE_WARNING_THRESHOLD for x in std_errors)

    status_parts = []
    if tail_prob is None:
        status_parts.append("missing_p")
    if not coeff_len_ok:
        status_parts.append("bad_coeff_len")
    if not coeffs_valid:
        status_parts.append("invalid_coefficients")
    if std_errors and not se_len_ok:
        status_parts.append("bad_se_len")
    if se_suspicious:
        status_parts.append("unstable_se")

    status = "ok" if not status_parts else ";".join(status_parts)

    rows.append({
        "target": target,
        "combo_id": combo_id,
        "n_sources": str(meta["n_sources"]),
        "sources": meta["sources"],
        "tail_prob": "" if tail_prob is None else f"{tail_prob:.6g}",
        "coefficients": fmt_float_list(coeffs),
        "std_errors": fmt_float_list(std_errors),
        "coeffs_valid": "true" if coeffs_valid else "false",
        "has_negative": "true" if has_negative else "false",
        "has_gt_one": "true" if has_gt_one else "false",
        "se_suspicious": "true" if se_suspicious else "false",
        "status": status,
        "outfile": fp.name,
        "n_rights": str(meta["n_rights"]),
        "rights": meta["rights"],        
    })

summary_fields = [
    "target",
    "combo_id",
    "n_sources",
    "sources",
    "n_rights",
    "rights",
    "tail_prob",
    "coefficients",
    "std_errors",
    "coeffs_valid",
    "has_negative",
    "has_gt_one",
    "se_suspicious",
    "status",
    "outfile",
]

def p_as_float(row):
    try:
        return float(row["tail_prob"]) if row["tail_prob"] != "" else float("-inf")
    except ValueError:
        return float("-inf")

with open(summary_out, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter="\t")
    writer.writeheader()
    for row in sorted(rows, key=lambda r: (r["target"], -p_as_float(r))):
        writer.writerow(row)

by_target = {}
for row in rows:
    by_target.setdefault(row["target"], []).append(row)

best_valid_rows = []
best_fallback_rows = []

for target, target_rows in sorted(by_target.items()):
    valid_rows = [r for r in target_rows if r["coeffs_valid"] == "true"]

    if valid_rows:
        best_valid = max(valid_rows, key=p_as_float)
        best_valid_rows.append(best_valid)
        best_fallback_rows.append(best_valid)
    elif target_rows:
        fallback = max(target_rows, key=p_as_float)
        best_fallback_rows.append(fallback)

with open(best_out, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter="\t")
    writer.writeheader()
    for row in best_valid_rows:
        writer.writerow(row)

with open(best_fallback_out, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter="\t")
    writer.writeheader()
    for row in best_fallback_rows:
        writer.writerow(row)

print(f"Wrote summary: {summary_out}")
print(f"Wrote best valid models: {best_out}")
print(f"Wrote best models with fallback: {best_fallback_out}")
print(f"Parsed qpAdm files: {len(rows)}")
PY
# -----------------------------
# final
# -----------------------------
log "DONE"
log "Main outputs:"
log "  merged dataset: ${OUTDIR}/merged/merged.{geno,snp,ind}"
log "  merged popcnts: ${OUTDIR}/stats/merged_pop_counts.txt"
log "  source config:  ${OUTDIR}/lists/target_source_candidates.tsv"
log "  source combos:  ${OUTDIR}/lists/source_combinations.tsv"
log "  qpAdm results:  ${OUTDIR}/qpadm/*.qpadm.txt"
log "  summary:        ${OUTDIR}/stats/qpadm_summary.tsv"
log "  best models:    ${OUTDIR}/stats/qpadm_best_models.tsv"
