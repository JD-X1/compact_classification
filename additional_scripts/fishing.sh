#!/usr/bin/env bash
set -euo pipefail

# get the number of threads
while getopts ':t:i:r:o:' option
do
    case "${option}" in
        t) TCores=${OPTARG} ;;
        i) Input=${OPTARG} ;;
        r) reSource=${OPTARG} ;;
        o) Outdir=${OPTARG} ;;
        \?) echo "Illegal option -${OPTARG}" >&2; exit 2 ;;
        :)  echo "Option -${OPTARG} requires an argument" >&2; exit 2 ;;
    esac
done

orig_pwd="$(pwd)"

outdir="${Outdir:-.}"; outdir="${outdir%/}/"
[[ "${outdir}" = /* ]] || outdir="${orig_pwd}/${outdir}"
mkdir -p "${outdir}"

if [[ "${Input}" = /* ]]; then
  input_abs="${Input}"
else
  input_abs="${orig_pwd}/${Input}"
fi

resource_abs=""
if [[ -n "${reSource:-}" ]]; then
  if [[ "${reSource}" = /* ]]; then
    resource_abs="${reSource%/}"
  else
    resource_abs="${orig_pwd}/${reSource}"
  fi
fi

cd "${outdir}"
mag="$(basename "${input_abs}" _input_metadata.tsv)"
log_dir="logs/FishingLogs"
phyloscratch_dir="${mag}_PhyloFishScratch"
fish_out="${mag}_fish_out"
working_dataset="${mag}_working_dataset"

mkdir -p "${log_dir}"

echo "CWD: $(pwd)"
echo "outdir: ${outdir}"
echo "input abs: ${input_abs}"
echo "Fishing for ${mag}"
echo "log outdir: ${log_dir}"
echo "Fisher Out: ${fish_out}"
echo "phyloscratch dir: ${phyloscratch_dir}"

# A predictor that returned no residues (e.g. metaeuk found no genes) makes
# fisher.py divide by zero in is_aa_seq. Detect that here and fall through to an
# empty working_dataset so the workflow routes the MAG to UNCLASSIFIABLE.
residues=0
while IFS=$'\t' read -r pdir pfile _; do
  [[ -z "${pfile:-}" ]] && continue
  pf="${pdir}/${pfile}"
  [[ -f "${pf}" ]] || continue
  r=$(grep -v '^>' "${pf}" 2>/dev/null | tr -d '[:space:]' | wc -c)
  residues=$((residues + r))
done < "${input_abs}"
if [[ "${residues}" -eq 0 ]]; then
  echo "No predicted residues for ${mag}; skipping fisher, emitting empty working_dataset" >&2
  mkdir -p "${working_dataset}"
  exit 0
fi

phyloDB=""

for cand in \
  "${resource_abs:+${resource_abs}/}PhyloFisherDatabase_v1.0/database" \
  "/compact_classification/resources/PhyloFisherDatabase_v1.0/database" \
  "/PhyloFisherDatabase_v1.0/database" \
  "PhyloFisherDatabase_v1.0/database"
do
  [[ -n "${cand}" && -d "${cand}" ]] && { phyloDB="${cand}"; break; }
done

if [[ -z "${phyloDB}" ]]; then
  echo "ERROR: Could not find PhyloFisher DB (tried -r, /compact_classification/resources, /)." >&2
  exit 2
fi

echo "Gathering Bait"
if [[ ! -d "${phyloscratch_dir}" ]]; then
  mkdir -p "${phyloscratch_dir}"
  echo "Creating the PhyloFishScratch database"
  cp -r "${phyloDB}" "${phyloscratch_dir}/."
  echo "PhyloFishScratch database created"
fi

if [[ ! -d "${outdir}logs" ]]; then
  mkdir -p "${outdir}logs"
fi

echo "Casting Lines"

cp -f "${input_abs}" "${phyloscratch_dir}/${mag}_metadata.tsv"
config.py -d "${phyloscratch_dir}" -i "${phyloscratch_dir}/${mag}_metadata.tsv"
echo "Configuration of PhyloFisher Modules Complete"

echo "Waiting for the Fish to Bite"
fisher.py --threads "${TCores}" -o "${fish_out}" --keep_tmp \
    1> "${log_dir}/${mag}_fisher.log" 2>&1
echo "Fish Caught"


# informant.py is diagnostic-only and raises on a zero-ortholog matrix
# ("cannot set a frame with no defined columns"). A marker-empty MAG must not
# hard-fail here -- let it fall through to an empty working_dataset so the
# workflow can route it to a graceful UNCLASSIFIABLE verdict.
informant.py -i "${fish_out}" --orthologs_only \
    1> "${log_dir}/${mag}_informant.log" 2>&1 \
    || echo "informant.py exited non-zero (likely zero ortholog hits); continuing" >&2
echo "Informant Complete"

echo "Choosing the best fish"
working_dataset_constructor.py -i "${fish_out}" -o "${working_dataset}" \
    1> "${log_dir}/${mag}_working_dataset_constructor.log" 2>&1 \
    || echo "working_dataset_constructor.py exited non-zero; emitting empty working_dataset" >&2
# Guarantee the checkpoint output dir exists even when no markers were recovered,
# so goneFishing succeeds and the empty gene list drives the UNCLASSIFIABLE ramp.
mkdir -p "${working_dataset}"
echo "Fish on the grill"

[[ -f config.ini ]] && mv -f config.ini "${phyloscratch_dir}/config.ini"