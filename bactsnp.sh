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
trap 'kill $(jobs -p)' 2

export ROOT_PATH="ROOT_PATH_DUMMY"
export SUB_EXE_PATH=${ROOT_PATH}/sub_exe
export OTHERS_PATH=${ROOT_PATH}/others
export PATH=${SUB_EXE_PATH}:$PATH

source ${SUB_EXE_PATH}/read_options.sh
echo "detailed log is printed in <output directory>/bactsnp_out"

module=simulate_reads
  if [ -n "${FA_LST}" ]; then
    echo "started '${module}' module"
    echo "started '${module}' module" >>${BACTSNP_OUT}
    set +e; xargs -a ${FA_LST} -r -i -P ${JOBS} ${module}.sh {}
    set -e; check_job_status.sh ${module} art ${FA_ISOLATE_LST}
    set -e; check_job_status.sh ${module} sam2fq ${FA_ISOLATE_LST}
    echo "finished '${module}' module"
    echo "finished '${module}' module" >>${BACTSNP_OUT}
  fi

module=trim_reads
  echo "started '${module}' module"
  echo "started '${module}' module" >>${BACTSNP_OUT}
  set +e; xargs -a ${FQ_LST} -r -i -P ${JOBS} ${module}.sh {}
  set -e; check_job_status.sh ${module} platanus_trim ${ALL_ISOLATE_LST}
  echo "finished '${module}' module"
  echo "finished '${module}' module" >>${BACTSNP_OUT}

module=assemble_reads
  echo "started '${module}' module"
  echo "started '${module}' module" >>${BACTSNP_OUT}
  set +e; xargs -a ${ALL_ISOLATE_LST} -r -i -P ${JOBS} ${module}.sh {}
  set -e; check_job_status.sh ${module} platanus ${ALL_ISOLATE_LST}
  echo "finished '${module}' module"
  echo "finished '${module}' module" >>${BACTSNP_OUT}

module=map_reads
  echo "started '${module}' module"
  echo "started '${module}' module" >>${BACTSNP_OUT}
  set +e; ${module}.sh index_ref
  set -e; check_job_status.sh ${module} bwa_index
  set +e; xargs -a ${ALL_ISOLATE_LST} -r -i -P ${JOBS} ${module}.sh map {}
  set -e; check_job_status.sh ${module} bwa_mem ${ALL_ISOLATE_LST}
  set +e; xargs -a ${ALL_ISOLATE_LST} -r -i -P ${JOBS} ${module}.sh rm_duplicate {}
  set -e; check_job_status.sh ${module} markduplicates ${ALL_ISOLATE_LST}
  echo "finished '${module}' module"
  echo "finished '${module}' module" >>${BACTSNP_OUT}

module=get_pseudo_genomes
  echo "started '${module}' module"
  echo "started '${module}' module" >>${BACTSNP_OUT}
  set +e; xargs -a ${ALL_ISOLATE_LST} -r -i -P ${JOBS} ${module}.sh run_mapping_mode {}
  set -e; check_job_status.sh ${module} samtools_view ${ALL_ISOLATE_LST}
  set -e; check_job_status.sh ${module} sam2pseudo ${ALL_ISOLATE_LST}
  set +e; xargs -a ${ALL_ISOLATE_LST} -r -i -P ${JOBS} ${module}.sh run_assembly_mode {}
  set -e; check_job_status.sh ${module} nucmer ${ALL_ISOLATE_LST}
  set -e; check_job_status.sh ${module} delta2pseudo ${ALL_ISOLATE_LST}
  set +e; xargs -a ${ALL_ISOLATE_LST} -r -i -P ${JOBS} ${module}.sh get_consensus {}
  set -e; check_job_status.sh ${module} merge_fa ${ALL_ISOLATE_LST}
  if [ -n "${INPUT_REGION}" ]; then
    set +e; ${module}.sh merge_input_region
    set -e; check_job_status.sh ${module} merge_region_input
    set +e; xargs -a ${ALL_ISOLATE_LST} -r -i -P ${JOBS} ${module}.sh mask_input_region {}
    set -e; check_job_status.sh ${module} mask_region_input ${ALL_ISOLATE_LST}
  fi
  set -e; ${module}.sh get_multi_fa
  echo "finished '${module}' module"
  echo "finished '${module}' module" >>${BACTSNP_OUT}

module=get_snps
  echo "started '${module}' module"
  echo "started '${module}' module" >>${BACTSNP_OUT}
  set +e; ${module}.sh
  set -e; check_job_status.sh ${module} fa2snp_wo_ref
  set -e; check_job_status.sh ${module} fa2snp_w_ref
  echo "finished '${module}' module"
  echo "finished '${module}' module" >>${BACTSNP_OUT}

module=get_replaced_pseudo_genomes
  echo "started '${module}' module"
  echo "started '${module}' module" >>${BACTSNP_OUT}
  set +e; ${module}.sh
  set -e; check_job_status.sh ${module} snp2fa
  echo "finished '${module}' module"
  echo "finished '${module}' module" >>${BACTSNP_OUT}

mkdir -p ${OUT_DIR}/mapping_results
mkdir -p ${OUT_DIR}/assembly_results
for isolate in $(cat ${ALL_ISOLATE_LST}); do
  echo "mv ${TMP_DIR}/assemble_reads/${isolate}/platanus_contig.fa ${OUT_DIR}/assembly_results/${isolate}.fa"    >>${BACTSNP_OUT}
        mv ${TMP_DIR}/assemble_reads/${isolate}/platanus_contig.fa ${OUT_DIR}/assembly_results/${isolate}.fa
  echo "mv ${TMP_DIR}/map_reads/${isolate}/markduplicates.bam     ${OUT_DIR}/mapping_results/${isolate}.bam"     >>${BACTSNP_OUT}
        mv ${TMP_DIR}/map_reads/${isolate}/markduplicates.bam     ${OUT_DIR}/mapping_results/${isolate}.bam
  echo "mv ${TMP_DIR}/map_reads/${isolate}/markduplicates.bam.bai ${OUT_DIR}/mapping_results/${isolate}.bam.bai" >>${BACTSNP_OUT}
        mv ${TMP_DIR}/map_reads/${isolate}/markduplicates.bam.bai ${OUT_DIR}/mapping_results/${isolate}.bam.bai
done

if [ "${NO_CLEAN}" == false ]; then
  if [ $REF == ${TMP_DIR}/ref.fa ]; then
    # for copy of the write-protected reference
    rm -rf $REF
  fi
  rm -r ${TMP_DIR}
fi

echo >>${BACTSNP_OUT}
echo "bactsnp completed!"
echo "bactsnp completed!" >>${BACTSNP_OUT}
