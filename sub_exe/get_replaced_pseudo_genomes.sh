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

snp_=${OUT_DIR}/snps_wo_ref.tsv
#out_dir=${TMP_DIR}/${module}
log_dir=${LOG_DIR}/${module}
snp2fa_log_prefix=${log_dir}/snp2fa
#mkdir -p ${out_dir}
mkdir -p ${log_dir}
mkdir -p ${OUT_DIR}/replaced_pseudo_genome
echo "snp2fa ${snp_} ${REF} ${OUT_DIR}/replaced_pseudo_genome 2>${snp2fa_log_prefix}_stderr" >>${BACTSNP_OUT}
      snp2fa ${snp_} ${REF} ${OUT_DIR}/replaced_pseudo_genome 2>${snp2fa_log_prefix}_stderr
for isolate in $(cat ${ALL_ISOLATE_LST}); do
  echo ">${isolate}"
  grep -v '>' ${OUT_DIR}/replaced_pseudo_genome/${isolate}.fa | perl -pe 's/\n//' | perl -pe 's/(.{80})/\1\n/g'
  echo
done >${OUT_DIR}/replaced_pseudo_genomes_wo_ref.fa
{ echo ">ref"
  grep -v '>' ${REF} | perl -pe 's/\n//' | perl -pe 's/(.{80})/\1\n/g'
  echo
  cat ${OUT_DIR}/replaced_pseudo_genomes_wo_ref.fa
} >${OUT_DIR}/replaced_pseudo_genomes_w_ref.fa
touch ${snp2fa_log_prefix}.done
