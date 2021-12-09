#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=01:30:00
#SBATCH --constraint=haswell
#SBATCH --output=BATCH_OUTPUT
#SBATCH --error=BATCH_OUTPUT


EXE=/global/homes/c/cjain7/shared/tools/seqrequester/seqrequester
GENOME=/global/homes/c/cjain7/shared/projects/assembly/chirag/data/genomes/chm13.draft_v1.1.fasta
GENOMESIZE=/global/homes/c/cjain7/shared/projects/assembly/chirag/data/genomes/chm13.draft_v1.1.fasta.size
COVERAGE=20
DIST=/global/homes/c/cjain7/shared/seq/human/ONT/HG02080_HPRC/HG02080_2.fastq.hist.txt
SIZE=`grep -v ">" $GENOME | wc | awk '{print $3-$1}'`
OUTPUTSEQ=reads.fasta
SEQKIT=/global/homes/c/cjain7/shared/tools/seqKit-v2.0.0/seqkit

$EXE simulate -truncate -genome $GENOME -genomesize $SIZE -coverage $COVERAGE -distribution $DIST > $OUTPUTSEQ

grep ">" $OUTPUTSEQ > ${OUTPUTSEQ}.headers
cat ${OUTPUTSEQ}.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > ${OUTPUTSEQ}.headers.bed
python countcontainedreads.py ${OUTPUTSEQ}.headers.bed

bedtools genomecov -i ${OUTPUTSEQ}.headers.bed -g $GENOMESIZE -bga | awk '{if ($4==0) print $0}' > ${OUTPUTSEQ}.headers.nocov.bed
bedtools genomecov -i ${OUTPUTSEQ}.headers.bed -g $GENOMESIZE > ${OUTPUTSEQ}.headers.cov
cat $OUTPUTSEQ | head -n 1000000 | $SEQKIT seq -g -m 1000 | $SEQKIT watch --fields ReadLen -Q -O readlen_min1k.hist.pdf
cat $OUTPUTSEQ | $SEQKIT stats -a > ${OUTPUTSEQ}.stats
