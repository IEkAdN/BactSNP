#! /bin/bash

if [ $# != 3 ]; then
  echo "usage: `basename $0` [full path of paup executable] [stdout of show-aligns] <output prefix>"
  exit 1
fi

Paup=$1
Aligns=$2
Prfx=$3
Alignment=$Prfx.alignment
AlignmentFa=$Prfx.fa
Nex=$Prfx.nex
SeqNex=$Prfx\_seq.nex
PaupLog=$Prfx.paup_log
PaupOut=$Prfx.paup_out
PaupErr=$Prfx.paup_err
EvolutionaryModelParam=$Prfx.evolutionary_model_param
RefRootBranchLen=$Prfx.ref_root_branch_len

tail -n +3 $Aligns | grep -v "^ *$" | grep -v "=" | grep -v "^\-\-" | grep -v "\^" | awk '{print $2}' | perl -pe 's/\./-/g' | perl -pe 's/[^-acgtn\n]/N/g' | perl -pe 's/a/A/g' | perl -pe 's/c/C/g' | perl -pe 's/g/G/g' | perl -pe 's/t/T/g' >$Alignment
{ echo ">Root"
  awk 'NR % 2 == 1' $Alignment | perl -pe 's/\n//g' | perl -pe 's/(.+)/\1\n/'
  echo ">Ref"
  awk 'NR % 2 == 0' $Alignment | perl -pe 's/\n//g' | perl -pe 's/(.+)/\1\n/'
  echo ">Dummy"
  awk 'NR % 2 == 0' $Alignment | perl -pe 's/\n//g' | perl -pe 's/(.+)/\1\n/'
} >$AlignmentFa

cat <<EOS >$Nex
begin paup;
  set autoclose=yes warntree=no warnreset=no;
  tonexus replace=yes format=FASTA fromfile=$AlignmentFa tofile=$SeqNex;
  execute $SeqNex;
  nj;
  lset nst=6 rmatrix=estimate;
  log file=$PaupLog;
  lscores 1;
  showdist;
  quit warnTSave=no;
end;
EOS

$Paup -f $Nex >$PaupOut 2>$PaupErr

for Nuc in AC AG AT CG CT GT; do
  awk -v Nuc=$Nuc '$1 == Nuc {print $2}' $PaupLog | perl -pe 's/\n/,/'
done | perl -pe 's/(.+),/\1\n/' >$EvolutionaryModelParam

awk '$2 == "Dummy" {print $3}' $PaupLog >$RefRootBranchLen
