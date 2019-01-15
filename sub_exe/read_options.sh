#! /bin/bash

: <<"#__LICENSE__"
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
#__LICENSE__

set -eu

check_existence() {
  file=${1}
  if [ ! -e ${file} ]; then
    echo "ERROR: no such file or directory: ${file}" >&2
    exit 1
  fi
}

check_absence() {
  file=${1}
  if [ -e ${file} ]; then
    echo "ERROR: the file or directory already exists: ${file}" >&2
    exit 1
  fi
}

# set default parameters
fa_lst=""
fq_lst=""
ref=""
ref_strain=""
out_dir=""
thr=1
jobs=1
allele_freq=0.9
depth=10
input_region=""
dist_from_indel=5
no_clean=false

if [ $# -eq 0 ]; then
  cat ${ROOT_PATH}/README
  exit 0
fi

while [ $# -gt 0 ]; do
  opt=${1}
  case ${opt} in
    -h | --help )
      cat ${ROOT_PATH}/README
      exit 0
      ;;
    -v | --version )
      cat ${ROOT_PATH}/VERSION
      exit 0
      ;;
    -q | --fastq_list )
      fq_lst=$(readlink -m ${2})
      shift
      ;;
    -r | --reference )
      ref=$(readlink -m ${2})
      shift
      ;;
    -o | --out_dir )
      out_dir=$(readlink -m ${2})
      shift
      ;;
    -j | --jobs )
      jobs=${2}
      shift
      ;;
    -t | --thread )
      thr=${2}
      shift
      ;;
    -a | --fasta_list )
      fa_lst=$(readlink -m ${2})
      shift
      ;;
    --ref_strain )
      ref_strain=${2}
      shift
      ;;
    --mask_region )
      input_region=$(readlink -m ${2})
      shift
      ;;
    --dist_from_indel )
      dist_from_indel=${2}
      shift
      ;;
    --allele_freq )
      allele_freq=${2}
      shift
      ;;
    --depth )
      depth=${2}
      shift
      ;;
    --no_clean )
      no_clean=true
      ;;
    * )
      echo "ERROR" >&2
      echo "  ${1}: unrecognized option or file linked to none of the options." >&2
      exit 1
      ;;
  esac
  shift
done

if [ -z "${fq_lst}" ] && [ -z "${fa_lst}" ]; then
  echo "ERROR: either -q | --fastq_list or -a | --fasta_list is required" >&2
  exit 1
fi
if [ -z "${ref}" ] && [ -z "${ref_strain}" ]; then
  ref_strain=$(head -n 1 ${fq_lst} | cut -f 1)
fi
if [ -z "${out_dir}" ]; then
  echo "ERROR: -o | --out_dir is required" >&2
  exit 1
fi
if [ $(echo "${allele_freq} < 0" | bc) -eq 1 ] || [ $(echo "${allele_freq} > 1" | bc) -eq 1 ]; then
  echo "ERROR: minimum allele frequency must be specified between 0 and 1" >&2
  exit 1
fi
if [ $(echo "${depth} < 0" | bc) -eq 1 ] || [ $(echo ${depth} | grep '\.' | wc -l) -eq 1 ]; then
  echo "ERROR: minimum coverage depth must be specified as a positive integer" >&2
  exit 1
fi

# Do not set FQ_LST=${fq_lst}, because ${FQ_LST} is edited if the user input assembly data
# Do not set REF=${ref}. index files for the reference sequence should not be created in the directory where ${ref} exists
  readonly OUT_DIR=${out_dir};                                             export OUT_DIR
  readonly TMP_DIR=${OUT_DIR}/tmp;                                         export TMP_DIR
  readonly LOG_DIR=${TMP_DIR}/log;                                         export LOG_DIR
  readonly FA_LST=${fa_lst};                                               export FA_LST
  readonly FQ_LST=${TMP_DIR}/fq_list;                                      export FQ_LST
if [ -n "${ref}" ]; then
  readonly REF=${TMP_DIR}/ref.fa;                                          export REF
else
  readonly REF=${TMP_DIR}/assemble_reads/${ref_strain}/platanus_contig.fa; export REF
  readonly REF_STRAIN=${ref_strain};                                       export REF_STRAIN
fi
  readonly THR=${thr};                                                     export THR
  readonly JOBS=${jobs};                                                   export JOBS
  readonly ALLELE_FREQ=${allele_freq};                                     export ALLELE_FREQ
  readonly INPUT_REGION=${input_region};                                   export INPUT_REGION
  readonly DEPTH=${depth};                                                 export DEPTH
  readonly DIST_FROM_INDEL=${dist_from_indel};                             export DIST_FROM_INDEL
  readonly NO_CLEAN=${no_clean};                                           export NO_CLEAN
  readonly BACTSNP_OUT=${OUT_DIR}/bactsnp_out;                             export BACTSNP_OUT
  readonly BACTSNP_ERR=${OUT_DIR}/bactsnp_err;                             export BACTSNP_ERR
  readonly ALL_ISOLATE_LST=${TMP_DIR}/all_isolate_list;                    export ALL_ISOLATE_LST
  readonly FQ_ISOLATE_LST=${TMP_DIR}/fq_isolate_list;                      export FQ_ISOLATE_LST
  readonly FA_ISOLATE_LST=${TMP_DIR}/fa_isolate_list;                      export FA_ISOLATE_LST

which platanus_trim >/dev/null
which bwa >/dev/null
which samtools >/dev/null
which platanus >/dev/null
which nucmer >/dev/null
if [ -n "${FA_LST}" ]; then
  which art_illumina >/dev/null
fi

if [ -n "${ref}" ]; then
  check_existence ${ref}
fi
check_absence ${OUT_DIR}
if [ -n "${fq_lst}" ]; then
  check_existence ${fq_lst}
  # check if ${fq_lst} has 3 columns per line
  if [ -n "$(awk -F '\t' 'NF != 3' ${fq_lst})" ]; then
    echo "ERROR: the format of fastq list file should be as follows" >&2
    cat ${OTHERS_PATH}/fq_list_format                                >&2
    exit 1
  fi
  # check if fastq files in ${fq_lst} exist
  cat ${fq_lst} | while read isolate fwd_fq rev_fq; do
    check_existence $(readlink -m ${fwd_fq})
    check_existence $(readlink -m ${rev_fq})
  done
fi
if [ -n "${FA_LST}" ]; then
  check_existence ${FA_LST}
  # check if ${FA_LST} has 2 columns per line
  if [ -n "$(awk -F '\t' 'NF != 2' ${FA_LST})" ]; then
    echo "ERROR: the format of fasta list file should be as follows" >&2
    cat ${OTHERS_PATH}/fa_list_format                                >&2
    exit 1
  fi
  # check if fasta files in ${FA_LST} exist
  cat ${FA_LST} | while read isolate fa; do
    check_existence $(readlink -m ${fa})
  done
fi
if [ -n "${INPUT_REGION}" ]; then
  check_existence ${INPUT_REGION}
  # check if ${INPUT_REGION} has 3 columns per line
  if [ -n "$(awk -F '\t' 'NF != 3' ${INPUT_REGION})" ]; then
    echo "ERROR: the format of mask-region file should be as follows" >&2
    cat ${OTHERS_PATH}/masked_region_format                           >&2
    exit 1
  fi
fi
# check if ${fq_lst} and ${FA_LST} contain the same isolates
  if [ $(cat ${fq_lst} ${FA_LST} | cut -f 1 | sort -u | wc -l) -lt $(cat ${fq_lst} ${FA_LST} | cut -f 1 | wc -l) ]; then
    echo "ERROR: ${fq_lst} and ${FA_LST} contain the same isolates" >&2
    exit 1
  fi
# check if ${fq_lst} and ${FA_LST} contain an isolate named "ref"
  if [ $(cat ${fq_lst} ${FA_LST} | cut -f 1 | grep "^ref$" | wc -l) -ne 0 ]; then
    echo "ERROR: \"ref\" can not be used as an isolate name" >&2
    exit 1
  fi

mkdir -p ${LOG_DIR}
mkdir -p ${TMP_DIR}
# Do not use 'if [ -n "${REF}" ]; then'
if [ -n "${ref}" ]; then
  cp ${ref} ${REF}
fi
{
    echo   "fastq list                                         : ${fq_lst}"
    echo   "fasta list                                         : ${FA_LST}"
  # ${REF_STRAIN} is not set when uesr input a ref
  if [ -z "${ref_strain}" ]; then
    echo   "reference                                          : ${ref}"
  else
    echo   "reference                                          : ${OUT_DIR}/assembly_results/${REF_STRAIN}.fa"
  fi
    echo   "output directory                                   : ${OUT_DIR}"
    echo   "number of threads                                  : ${THR}"
    echo   "number of concurrent jobs                          : ${JOBS}"
    echo   "no_clean option                                    : ${NO_CLEAN}"
  if [ -n "${INPUT_REGION}" ]; then
    echo   "masked region                                      : ${INPUT_REGION}"
  fi
    echo   "thresholds:"
    echo   "    minimum distance from the nearest indels       : ${DIST_FROM_INDEL}"
    echo   "    minimum coverage depth                         : ${DEPTH}"
    echo   "    minimum allele frequency                       : ${ALLELE_FREQ}"
    echo
} >${BACTSNP_OUT}

if [ -n "${fq_lst}" ]; then
  cp ${fq_lst} ${FQ_LST}
  cut -f 1 ${FQ_LST} >${FQ_ISOLATE_LST}
else
  touch ${FQ_LST}
  touch ${FQ_ISOLATE_LST}
fi
if [ -n "${FA_LST}" ]; then
  cut -f 1 ${FA_LST} >${FA_ISOLATE_LST}
else
  touch ${FA_ISOLATE_LST}
fi
cat ${FQ_ISOLATE_LST} ${FA_ISOLATE_LST} >${ALL_ISOLATE_LST}
if [ -n "${FA_LST}" ]; then
  for isolate in $(cat ${FA_ISOLATE_LST}); do
    echo "${isolate}	${TMP_DIR}/simulate_reads/${isolate}/R1.fq	${TMP_DIR}/simulate_reads/${isolate}/R2.fq" >>${FQ_LST}
  done
fi
