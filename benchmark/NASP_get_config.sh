#! /bin/bash

if [ $# != 5 ]; then
  echo "usage: `basename $0` [reference.fa] [read_directory] [fastq_list] [output_directory] <run name>"
  exit 1
fi

RefFa=$(echo $1 | perl -pe 's/\//\\\//g')
ReadDir=$(echo $2 | perl -pe 's/\//\\\//g')
Iso_R1Fq_R2Fq=$3
OutDir=$(echo $4 | perl -pe 's/\//\\\//g')
RunName=$5

if [ `echo $RefFa | perl -pe 's/.+(...)/\1/'` != ".fa" ]; then
  echo "ERROR: invalid reference fasta: $RefFa" >&2
  exit 1
fi
RefPrfx=`basename $RefFa .fa`

cat <<'EOS' | perl -pe "s/__RunName__/<RunName>$RunName<\/RunName>/" | \
              perl -pe "s/__OutDir__/<OutputFolder>$OutDir<\/OutputFolder>/" | \
              perl -pe "s/__Ref__/<Reference name=\"$RefPrfx\" path=\"$RefFa\">/" | \
              perl -pe "s/__ReadDir__/<ReadFolder path=\"$ReadDir\">/"
<?xml version="1.0" ?>
<NaspInputData>
    <Options>
        __RunName__
        __OutDir__
        __Ref__
            <FindDups>True</FindDups>
        </Reference>
        <Filters>
            <ProportionFilter>0.9</ProportionFilter>
            <CoverageFilter>3</CoverageFilter>
        </Filters>
        <JobSubmitter>SGE</JobSubmitter>
    </Options>
    <Files>
        __ReadDir__
EOS

cat $Iso_R1Fq_R2Fq | while read Iso R1Fq R2Fq; do
  echo "            <ReadPair sample=\"$Iso\">"
  echo "                <Read1Filename>$R1Fq</Read1Filename>"
  echo "                <Read2Filename>$R2Fq</Read2Filename>"
  echo "            </ReadPair>"
done
  
cat <<'EOS'
        </ReadFolder>
    </Files>
    <ExternalApplications>
        <Index name="Index" path="/data/dai/tools/nasp_1.1.2/bin">
            <AdditionalArguments/>
            <JobParameters name="nasp_index">
                <MemRequested>2</MemRequested>
                <NumCPUs>1</NumCPUs>
                <Walltime>4</Walltime>
                <Queue>sge100.q</Queue>
                <JobSubmitterArgs/>
            </JobParameters>
        </Index>
        <MatrixGenerator name="MatrixGenerator" path="/data/dai/tools/nasp_1.1.2/lib/python3.4/site-packages/nasp/nasptool_linux_64">
            <AdditionalArguments/>
            <JobParameters name="nasp_matrix">
                <MemRequested>8</MemRequested>
                <NumCPUs>24</NumCPUs>
                <Walltime>48</Walltime>
                <Queue>sge100.q</Queue>
                <JobSubmitterArgs/>
            </JobParameters>
        </MatrixGenerator>
        <Picard name="Picard" path="/data/dai/tools/picard-2.18.17/picard.jar">
            <AdditionalArguments/>
        </Picard>
        <Samtools name="Samtools" path="/data/dai/tools/samtools-1.2/samtools">
            <AdditionalArguments/>
        </Samtools>
        <DupFinder name="DupFinder" path="/data/dai/tools/MUMmer3.23/nucmer">
            <AdditionalArguments/>
            <JobParameters>
                <MemRequested>4</MemRequested>
                <NumCPUs>1</NumCPUs>
                <Walltime>4</Walltime>
                <Queue>sge100.q</Queue>
                <JobSubmitterArgs/>
            </JobParameters>
        </DupFinder>
        <Aligner name="BWA-mem" path="/data/dai/tools/bwa-0.7.15/bwa">
            <AdditionalArguments/>
            <JobParameters>
                <MemRequested>3</MemRequested>
                <NumCPUs>1</NumCPUs>
                <Walltime>36</Walltime>
                <Queue>sge100.q</Queue>
                <JobSubmitterArgs/>
            </JobParameters>
        </Aligner>
        <SNPCaller name="GATK" path="/data/dai/tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar">
            <AdditionalArguments>-stand_call_conf 100 -ploidy 1</AdditionalArguments>
            <JobParameters>
                <MemRequested>3</MemRequested>
                <NumCPUs>1</NumCPUs>
                <Walltime>36</Walltime>
                <Queue>sge100.q</Queue>
                <JobSubmitterArgs/>
            </JobParameters>
        </SNPCaller>
    </ExternalApplications>
</NaspInputData>
EOS
