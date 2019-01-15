#! /bin/bash

if [ $# != 5 ]; then
  echo "usage: `basename $0` [in.newick] [in.Ref.fa] [in.SubstFreq] <# variable sites> [out.dir]"
  exit 1
fi

Tree=$1
RefFa=$(readlink -f $2)
SubstFreq=$(readlink -f $3)
VarSiteNum=$4
OutD=$(readlink -f $5)

RefId=Ref
Cov=40
ReadLen=250
FragmentSize=500
FragmentSizeStdev=50
PercentClustered=0.25
ExponentailMean=125
IndelModel="LAV 1.7 541"
IndelRate=0.1

cat <<EOS
treefile_path = $Tree
base_genome_name = $RefId
base_genome_path = $RefFa
number_of_variable_sites = $VarSiteNum
output_dir = $OutD
rate_matrix = $(cat $SubstFreq)
coverage = $Cov
read_length = $ReadLen
fragment_size = $FragmentSize
stdev_frag_size = $FragmentSizeStdev
mutation_clustering = ON
percent_clustered = $PercentClustered
exponential_mean = $ExponentailMean
indel_model = $IndelModel
indel_rate = $IndelRate
EOS
