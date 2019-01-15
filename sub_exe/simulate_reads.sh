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
fa_lst_line=$1

isolate=$(echo "${fa_lst_line}" | cut -f 1)
fa=$(echo "${fa_lst_line}" | cut -f 2)
eval fa=${fa}  # expand '~'
fa=$(readlink -f ${fa})
out_dir=${TMP_DIR}/${module}/${isolate}
log_dir=${LOG_DIR}/${module}/${isolate}
art_log_prefix=${log_dir}/art
sam2fq_log_prefix=${log_dir}/sam2fq
mkdir -p ${out_dir}
mkdir -p ${log_dir}
echo "art_illumina -i ${fa} -l 250 -f 40 -ss MSv3 -m 500 -s 50 -na -rs 1 -ef -o ${out_dir}/R >${art_log_prefix}_stdout 2>${art_log_prefix}_stderr" >>${BACTSNP_OUT}
      art_illumina -i ${fa} -l 250 -f 40 -ss MSv3 -m 500 -s 50 -na -rs 1 -ef -o ${out_dir}/R >${art_log_prefix}_stdout 2>${art_log_prefix}_stderr
touch ${art_log_prefix}.done
echo "rm ${out_dir}/R1.fq" >>${BACTSNP_OUT}
      rm ${out_dir}/R1.fq
echo "rm ${out_dir}/R2.fq" >>${BACTSNP_OUT}
      rm ${out_dir}/R2.fq
echo "rm ${out_dir}/R.sam" >>${BACTSNP_OUT}
      rm ${out_dir}/R.sam
echo "sam2fq ${out_dir}/R_errFree.sam ${out_dir}/R1.fq ${out_dir}/R2.fq 2>${sam2fq_log_prefix}_stderr" >>${BACTSNP_OUT}
      sam2fq ${out_dir}/R_errFree.sam ${out_dir}/R1.fq ${out_dir}/R2.fq 2>${sam2fq_log_prefix}_stderr
touch ${sam2fq_log_prefix}.done
echo "rm ${out_dir}/R_errFree.sam" >>${BACTSNP_OUT}
      rm ${out_dir}/R_errFree.sam
