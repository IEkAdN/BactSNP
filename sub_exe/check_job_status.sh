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

set -e

module=$1
job=$2
if [ $# -eq 3 ]; then
  isolate_lst=$3
else
  isolate_lst=""
fi

PrintErrorMessage() {
  job_printed="$1"
  echo "ERROR: ${job_printed} in '${module}' module for isolate '${isolate}' failed"
  echo "       check log files below"
}

PrintLogName() {
  prefix=$1
  for log in ${prefix}*; do
    if [ -s ${log} ]; then
      echo "       ${log}"
    fi
  done
}


job_status=0
if [ -n "${isolate_lst}" ]; then
  for isolate in $(cat ${isolate_lst}); do
    if [ ! -e ${LOG_DIR}/${module}/${isolate}/${job}.done ]; then
      if [ ${job} == bwa_mem ]; then
        PrintErrorMessage "bwa_mem or samtools_sort"
        for job_ in bwa_mem samtools_sort; do
          PrintLogName ${LOG_DIR}/${module}/${isolate}/${job_}
        done
      else
        if [ ${job} == platanus ]; then
          echo "cp ${OTHERS_PATH}/platanus_warning ${LOG_DIR}/${module}/${isolate}/" >>${BACTSNP_OUT}
                cp ${OTHERS_PATH}/platanus_warning ${LOG_DIR}/${module}/${isolate}/
        fi
        PrintErrorMessage ${job}
        PrintLogName ${LOG_DIR}/${module}/${isolate}/${job}
      fi
      job_status=1
    fi
  done
else
  if [ ! -e ${LOG_DIR}/${module}/${job}.done ]; then
    echo "ERROR: ${job} in '${module}' step failed"
    echo "       check log files below"
    PrintLogName ${LOG_DIR}/${module}/${job}
    job_status=1
  fi
fi >>${BACTSNP_ERR}

if [ ${job_status} -eq 1 ]; then
  echo "BactSNP has failed. See ${BACTSNP_ERR}" >&2
  exit 1
fi
