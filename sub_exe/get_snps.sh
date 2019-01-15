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

out_dir=${TMP_DIR}/${module}
log_dir=${LOG_DIR}/${module}
fa2snp_wo_ref_log_prefix=${log_dir}/fa2snp_wo_ref
fa2snp_w_ref_log_prefix=${log_dir}/fa2snp_w_ref
mkdir -p ${out_dir}
mkdir -p ${log_dir}
for isolate in $(cat ${ALL_ISOLATE_LST}); do
  echo "${isolate}	${OUT_DIR}/pseudo_genome/${isolate}.fa"
done >${out_dir}/wo_ref_pseudo_genome_list
{ echo "ref	${REF}"
  cat ${out_dir}/wo_ref_pseudo_genome_list
} >${out_dir}/w_ref_pseudo_genome_list
echo "fa2snp ${out_dir}/wo_ref_pseudo_genome_list >${OUT_DIR}/snps_wo_ref.tsv 2>${fa2snp_wo_ref_log_prefix}_stderr" >>${BACTSNP_OUT}
      fa2snp ${out_dir}/wo_ref_pseudo_genome_list >${OUT_DIR}/snps_wo_ref.tsv 2>${fa2snp_wo_ref_log_prefix}_stderr
touch ${fa2snp_wo_ref_log_prefix}.done
echo "fa2snp ${out_dir}/w_ref_pseudo_genome_list  >${OUT_DIR}/snps_w_ref.tsv  2>${fa2snp_w_ref_log_prefix}_stderr"  >>${BACTSNP_OUT}
      fa2snp ${out_dir}/w_ref_pseudo_genome_list  >${OUT_DIR}/snps_w_ref.tsv  2>${fa2snp_w_ref_log_prefix}_stderr
touch ${fa2snp_w_ref_log_prefix}.done
