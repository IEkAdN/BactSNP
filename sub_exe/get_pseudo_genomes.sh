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

if [ ${submodule} == run_mapping_mode ]; then
  isolate=$2
  markduplicates_bam_=${TMP_DIR}/map_reads/${isolate}/markduplicates.bam
  out_dir=${TMP_DIR}/${module}/${isolate}
  log_dir=${LOG_DIR}/${module}/${isolate}
  markduplicates_sam=${out_dir}/markduplicates.sam
  view_log_prefix=${log_dir}/samtools_view
  sam2pseudo_dir=${out_dir}/sam2pseudo
  sam2pseudo_log_prefix=${log_dir}/sam2pseudo
  pseudo_genome_=${out_dir}/map_pseudo_genome.fa
  mkdir -p ${out_dir}
  mkdir -p ${log_dir}
  echo "samtools view ${markduplicates_bam_} >${markduplicates_sam} 2>${view_log_prefix}_stderr" >>${BACTSNP_OUT}
        samtools view ${markduplicates_bam_} >${markduplicates_sam} 2>${view_log_prefix}_stderr
  touch ${view_log_prefix}.done
  mkdir -p ${sam2pseudo_dir}
  echo "sam2pseudo ${markduplicates_sam} ${REF} ${sam2pseudo_dir} ${ALLELE_FREQ} ${DEPTH} >${pseudo_genome_} 2>${sam2pseudo_log_prefix}_stderr" >>${BACTSNP_OUT}
        sam2pseudo ${markduplicates_sam} ${REF} ${sam2pseudo_dir} ${ALLELE_FREQ} ${DEPTH} >${pseudo_genome_} 2>${sam2pseudo_log_prefix}_stderr
  touch ${sam2pseudo_log_prefix}.done
  echo "rm ${markduplicates_sam}" >>${BACTSNP_OUT}
        rm ${markduplicates_sam}
elif [ ${submodule} == run_assembly_mode ]; then
  isolate=$2
  out_dir=${TMP_DIR}/${module}/${isolate}
  log_dir=${LOG_DIR}/${module}/${isolate}
  mkdir -p ${out_dir}
  mkdir -p ${log_dir}
  platanus_assembly_=${TMP_DIR}/assemble_reads/${isolate}/platanus_contig.fa
  dist_from_edge=0
  nucmer_prefix=${out_dir}/nucmer
  nucmer_log_prefix=${log_dir}/nucmer
  delta=${nucmer_prefix}.delta
  pseudo_genome_=${out_dir}/ass_pseudo_genome.fa
  delta2pseudo_log_prefix=${log_dir}/delta2pseudo
  echo "nucmer -p ${nucmer_prefix} ${REF} ${platanus_assembly_} >${nucmer_log_prefix}_stdout 2>${nucmer_log_prefix}_stderr" >>${BACTSNP_OUT}
        nucmer -p ${nucmer_prefix} ${REF} ${platanus_assembly_} >${nucmer_log_prefix}_stdout 2>${nucmer_log_prefix}_stderr
  touch ${nucmer_log_prefix}.done
  echo "delta2pseudo ${REF} ${platanus_assembly_} ${delta} ${DIST_FROM_INDEL} >${pseudo_genome_} 2>${delta2pseudo_log_prefix}_stderr" >>${BACTSNP_OUT}
        delta2pseudo ${REF} ${platanus_assembly_} ${delta} ${DIST_FROM_INDEL} >${pseudo_genome_} 2>${delta2pseudo_log_prefix}_stderr
  touch ${delta2pseudo_log_prefix}.done
elif [ ${submodule} == get_consensus ]; then
  isolate=$2
  out_dir=${TMP_DIR}/${module}/${isolate}
  log_dir=${LOG_DIR}/${module}/${isolate}
  ass_pseudo_genome_=${out_dir}/ass_pseudo_genome.fa
  map_pseudo_genome_=${out_dir}/map_pseudo_genome.fa
  consensus_pseudo_genome_=${out_dir}/pseudo_genome.fa
  merge_fa_log_prefix=${log_dir}/merge_fa
  echo "merge_fa ${ass_pseudo_genome_} ${map_pseudo_genome_} >${consensus_pseudo_genome_} 2>${merge_fa_log_prefix}_stderr" >>${BACTSNP_OUT}
        merge_fa ${ass_pseudo_genome_} ${map_pseudo_genome_} >${consensus_pseudo_genome_} 2>${merge_fa_log_prefix}_stderr
  touch ${merge_fa_log_prefix}.done
elif [ ${submodule} == merge_input_region ]; then
  out_dir=${TMP_DIR}/${module}
  log_dir=${LOG_DIR}/${module}
  merge_region_log_prefix=${log_dir}/merge_region_input
  merged_input_region_=${out_dir}/merged_input.region
  echo "merge_region ${REF} ${INPUT_REGION} >${merged_input_region_} 2>${merge_region_log_prefix}" >>${BACTSNP_OUT}
        merge_region ${REF} ${INPUT_REGION} >${merged_input_region_} 2>${merge_region_log_prefix}
  touch ${merge_region_log_prefix}.done
elif [ ${submodule} == mask_input_region ]; then
  isolate=$2
  out_dir=${TMP_DIR}/${module}/${isolate}
  log_dir=${LOG_DIR}/${module}/${isolate}
  pseudo_genome_=${TMP_DIR}/${module}/${isolate}/pseudo_genome.fa
  masked_pseudo_genome_=${out_dir}/masked_pseudo_genome.fa
  merged_input_region_=${TMP_DIR}/${module}/merged_input.region
  mask_region_log_prefix=${log_dir}/mask_region_input
  echo "mask_region ${pseudo_genome_} ${merged_input_region_} >${masked_pseudo_genome_} 2>${mask_region_log_prefix}_stderr" >>${BACTSNP_OUT}
        mask_region ${pseudo_genome_} ${merged_input_region_} >${masked_pseudo_genome_} 2>${mask_region_log_prefix}_stderr
  touch ${mask_region_log_prefix}.done
elif [ ${submodule} == get_multi_fa ]; then
  out_dir=${TMP_DIR}/${module}
  mkdir -p ${OUT_DIR}/pseudo_genome
  { for isolate in $(cat ${ALL_ISOLATE_LST}); do
      if [ -n "${INPUT_REGION}" ]; then
        pseudo_genome_=${TMP_DIR}/${module}/${isolate}/masked_pseudo_genome.fa
      else
        pseudo_genome_=${TMP_DIR}/${module}/${isolate}/pseudo_genome.fa
      fi
      echo "mv ${pseudo_genome_} ${OUT_DIR}/pseudo_genome/${isolate}.fa" >>${BACTSNP_OUT}
            mv ${pseudo_genome_} ${OUT_DIR}/pseudo_genome/${isolate}.fa
      echo ">${isolate}"
      grep -v '>' ${OUT_DIR}/pseudo_genome/${isolate}.fa | perl -pe 's/\n//' | perl -pe 's/(.{80})/\1\n/g'
      echo
    done
  } >${OUT_DIR}/pseudo_genomes_wo_ref.fa
  { echo ">ref"
    grep -v '>' ${REF} | perl -pe 's/\n//' | perl -pe 's/(.{80})/\1\n/g'
    echo
    cat ${OUT_DIR}/pseudo_genomes_wo_ref.fa
  } >${OUT_DIR}/pseudo_genomes_w_ref.fa
fi
