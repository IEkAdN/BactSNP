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

module=$(basename $0 .sh)
submodule=$1
if [ $# -eq 2 ]; then
  isolate=$2
fi

if [ ${submodule} == index_ref ]; then
  out_dir=${TMP_DIR}/${module}
  log_dir=${LOG_DIR}/${module}
  index_prefix_=${out_dir}/$(basename ${REF})
  bwa_index_log_prefix=${log_dir}/bwa_index
  mkdir -p ${out_dir}
  mkdir -p ${log_dir}
  echo "bwa index -p ${index_prefix_} ${REF} >${bwa_index_log_prefix}_stdout 2>${bwa_index_log_prefix}_stderr" >>${BACTSNP_OUT}
        bwa index -p ${index_prefix_} ${REF} >${bwa_index_log_prefix}_stdout 2>${bwa_index_log_prefix}_stderr
  touch ${bwa_index_log_prefix}.done

elif [ ${submodule} == map ]; then
  fwd_trim_fq_=${TMP_DIR}/trim_reads/${isolate}/R1.fq.trimmed
  rev_trim_fq_=${TMP_DIR}/trim_reads/${isolate}/R2.fq.trimmed
  index_prefix_=${TMP_DIR}/${module}/$(basename ${REF})
  out_dir=${TMP_DIR}/${module}/${isolate}
  log_dir=${LOG_DIR}/${module}/${isolate}
  raw_bam_=${out_dir}/raw.bam
  mem_log_prefix=${log_dir}/bwa_mem
  sort_log_prefix=${log_dir}/samtools_sort
  mkdir -p ${out_dir}
  mkdir -p ${log_dir}
  {
    echo -n "bwa mem -t ${THR} -R "@RG	ID:${isolate}	SM:${isolate}" -M ${index_prefix_} ${fwd_trim_fq_} ${rev_trim_fq_} 2>${mem_log_prefix}_stderr | "
    echo    "samtools sort -@ ${THR} -o ${raw_bam_} - 2>${sort_log_prefix}_stderr"
  } >>${BACTSNP_OUT}
             bwa mem -t ${THR} -R "@RG	ID:${isolate}	SM:${isolate}" -M ${index_prefix_} ${fwd_trim_fq_} ${rev_trim_fq_} 2>${mem_log_prefix}_stderr | \
             samtools sort -@ ${THR} -o ${raw_bam_} - 2>${sort_log_prefix}_stderr
  if [ "${PIPESTATUS[*]}" == "0 0" ]; then
    touch ${mem_log_prefix}.done
  else
    exit 1
  fi
  echo "samtools index ${raw_bam_}" >>${BACTSNP_OUT}
        samtools index ${raw_bam_}

elif [ ${submodule} == rm_duplicate ]; then
  input_bam_=${TMP_DIR}/${module}/${isolate}/raw.bam
  out_dir=${TMP_DIR}/${module}/${isolate}
  log_dir=${LOG_DIR}/${module}/${isolate}
  markduplicates_bam_=${out_dir}/markduplicates.bam
  markduplicates_log_prefix=${log_dir}/markduplicates
  markduplicates_metrics=${out_dir}/markduplicates_metrics
  echo "java -Xmx2g -jar ${OTHERS_PATH}/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=${input_bam_} OUTPUT=${markduplicates_bam_} METRICS_FILE=${markduplicates_metrics} VALIDATION_STRINGENCY=LENIENT 2>${markduplicates_log_prefix}_stderr" >>${BACTSNP_OUT}
        java -Xmx2g -jar ${OTHERS_PATH}/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=${input_bam_} OUTPUT=${markduplicates_bam_} METRICS_FILE=${markduplicates_metrics} VALIDATION_STRINGENCY=LENIENT 2>${markduplicates_log_prefix}_stderr
  touch ${markduplicates_log_prefix}.done
  echo "samtools index ${markduplicates_bam_}" >>${BACTSNP_OUT}
        samtools index ${markduplicates_bam_}
  if [ ${NO_CLEAN} == false ]; then
    echo "rm ${input_bam_}" >>${BACTSNP_OUT}
          rm ${input_bam_}
  fi
fi
