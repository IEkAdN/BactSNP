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
fq_lst_line=$1

isolate=$(echo "${fq_lst_line}" | cut -f 1)
fwd_fq=$(echo  "${fq_lst_line}" | cut -f 2)
rev_fq=$(echo  "${fq_lst_line}" | cut -f 3)
eval fwd_fq=${fwd_fq}  # expand '~'
eval rev_fq=${rev_fq}  # expand '~'
fwd_fq=$(readlink -f ${fwd_fq})
rev_fq=$(readlink -f ${rev_fq})
out_dir=${TMP_DIR}/${module}/${isolate}
log_dir=${LOG_DIR}/${module}/${isolate}
fwd_cp_fq=${out_dir}/R1.fq
rev_cp_fq=${out_dir}/R2.fq
platanus_trim_log_prefix=${log_dir}/platanus_trim
mkdir -p ${out_dir}
mkdir -p ${log_dir}
if [ ${fwd_fq: -3} == .gz ]; then
  echo "gunzip -c ${fwd_fq} >${fwd_cp_fq}" >>${BACTSNP_OUT}
        gunzip -c ${fwd_fq} >${fwd_cp_fq}
else
  echo "cp ${fwd_fq} ${fwd_cp_fq}" >>${BACTSNP_OUT}
        cp ${fwd_fq} ${fwd_cp_fq}
fi
if [ ${rev_fq: -3} == .gz ]; then
  echo "gunzip -c ${rev_fq} >${rev_cp_fq}" >>${BACTSNP_OUT}
        gunzip -c ${rev_fq} >${rev_cp_fq}
else
  echo "cp ${rev_fq} ${rev_cp_fq}" >>${BACTSNP_OUT}
        cp ${rev_fq} ${rev_cp_fq}
fi
echo "platanus_trim ${fwd_cp_fq} ${rev_cp_fq} -t ${THR} >${platanus_trim_log_prefix}_stdout 2>${platanus_trim_log_prefix}_stderr" >>${BACTSNP_OUT}
      platanus_trim ${fwd_cp_fq} ${rev_cp_fq} -t ${THR} >${platanus_trim_log_prefix}_stdout 2>${platanus_trim_log_prefix}_stderr
touch ${platanus_trim_log_prefix}.done
rm -f ${fwd_cp_fq}
rm -f ${rev_cp_fq}
