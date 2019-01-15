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
isolate=$1

fwd_trim_fq_=${TMP_DIR}/trim_reads/${isolate}/R1.fq.trimmed
rev_trim_fq_=${TMP_DIR}/trim_reads/${isolate}/R2.fq.trimmed
out_dir=${TMP_DIR}/${module}/${isolate}
log_dir=${LOG_DIR}/${module}/${isolate}
platanus_log_prefix=${log_dir}/platanus
mkdir -p ${out_dir}
mkdir -p ${log_dir}
echo "platanus assemble -o ${out_dir}/platanus -f ${fwd_trim_fq_} ${rev_trim_fq_} -u 0 -m 2 -t ${THR} >${platanus_log_prefix}_stdout 2>${platanus_log_prefix}_stderr" >>${BACTSNP_OUT}
      platanus assemble -o ${out_dir}/platanus -f ${fwd_trim_fq_} ${rev_trim_fq_} -u 0 -m 2 -t ${THR} >${platanus_log_prefix}_stdout 2>${platanus_log_prefix}_stderr
touch ${platanus_log_prefix}.done
