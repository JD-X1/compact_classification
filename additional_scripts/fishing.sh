#!usr/bin/env bash
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

# get the file basename
mag="$(basename "${Input}" _input_metadata.tsv)"
outdir="${Outdir:-.}"; outdir="${outdir%/}/"

Input="$(basename "${Input}")"
log_dir="logs/FishingLogs"
phyloscratch_dir="${mag}_PhyloFishScratch"
fish_out="${mag}_fish_out"
working_dataset="${mag}_working_dataset"

mkdir -p "${log_dir}" "${phyloscratch_dir}"

echo "Fishing for ${mag}"
echo "log outdir: ${log_dir}"
echo "Fisher Out: ${fish_out}"
echo "phyloscratch dir: ${phyloscratch_dir}"

resource_root="${reSource:-}"; resources_root="${resource_root%/}/"
phyloDB=""

for cand in \
  "${resource_root}/PhyloFisherDatabase_v1.0/database" \
  "/compact_classification/resources/PhyloFisherDatabase_v1.0/database" \
  "/PhyloFisherDatabase_v1.0/database" \
  "PhyloFisherDatabase_v1.0/database"
do
  if [[ -d "${cand}" ]]; then phyloDB="${cand}"; break; fi
done

if [[ -z "${phyloDB}" ]]; then
  echo "ERROR: Could not find PhyloFisher DB (tried -r, /compact_classification/resources, /)." >&2
  exit 2
fi

cd "${outdir}"

echo "Gathering Bait"
if [[ ! -d "${phyloscratch_dir}/database" ]]; then
  echo "Creating the PhyloFishScratch database"
  cp -r "../${phyloDB}" "${phyloscratch_dir}"
  echo "PhyloFishScratch database created"
fi

echo "Casting Lines"

cp -f "${Input}" "${phyloscratch_dir}/metadata.tsv"
config.py -d "${phyloscratch_dir}/database" -i "${Input}"
echo "Configuration of PhyloFisher Modules Complete"

echo "Waiting for the Fish to Bite"
fisher.py --threads "${TCores}" -o "${fish_out}" --keep_tmp \
    1> "${log_dir}/${mag}_fisher.log" 2>&1
echo "Fish Caught"

informant.py -i "${fish_out}" --orthologs_only \
    1> "${log_dir}/${mag}_informant.log" 2>&1
echo "Informant Complete"

echo "Choosing the best fish"
working_dataset_constructor.py -i "${fish_out}" -o "${working_dataset}" \
    1> "${log_dir}/${mag}_working_dataset_constructor.log" 2>&1
echo "Fish on the grill"

mv config.ini "${phyloscratch_dir}/config.ini"

cd ..