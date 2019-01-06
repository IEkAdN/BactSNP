#! /bin/bash

if [ $# != 3 ]; then
  echo "usage: `basename $0` [stdout of EvolveAGene] [root_genome.fasta] <branch length between reference and root>"
  exit 1
fi

evolveagene_stdout=$1
root_genome_fasta=$2
branch_length_bw_ref_root=$3

root_genome_size=$(grep -v '>' $root_genome_fasta | perl -pe 's/\n//g' | wc -c)

sed -n '21,21'p $evolveagene_stdout | \
# convert the number of substitutions into the ratio of substitution to the root genome size
perl -pe 's/:(\d+)/:\n\1\n/g' | awk -v root_genome_size=$root_genome_size 'BEGIN {OFMT="%.13f"} {if ($0 ~ /^[0-9]+$/) {print $0 / root_genome_size} else {print $0}}' | perl -pe 's/\n//g' | perl -pe 's/(.+)/\1\n/' | \
# added a branch from the root to the reference isolate
perl -pe "s/(.+);/(\1, Ref:$branch_length_bw_ref_root);/"
