#Purpose: Simulate reads from a diploid genome using seqrequester 

HAP1=data/genomes/HG002.hifiasm.trio.0.16.1.hap1.fa
HAP2=data/genomes/HG002.hifiasm.trio.0.16.1.hap2.fa
COVERAGE_HAP=15
DIST=ONT/HG02080_HPRC/HG02080_2.fastq.hist.txt

#https://www.biostars.org/p/173963/#174150
GENOMESIZE=data/genomes/HG002.hifiasm.trio.0.16.1.size

SIZE_SUM=`grep -v ">" $HAP1 | wc | awk '{print $3-$1}'`
seqrequester simulate -truncate -genome $HAP1 -genomesize $SIZE_SUM -coverage $COVERAGE_HAP -distribution $DIST > reads.fasta

SIZE_SUM=`grep -v ">" $HAP2 | wc | awk '{print $3-$1}'`
seqrequester simulate -truncate -genome $HAP2 -genomesize $SIZE_SUM -coverage $COVERAGE_HAP -distribution $DIST >> reads.fasta

grep ">" reads.fasta > reads.fasta.headers
cat reads.fasta.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > reads.fasta.headers.bed
bedtools genomecov -i ${OUTPUTSEQ}.headers.bed -g $GENOMESIZE -bga | awk '{if ($4==0) print $0}' > reads.fasta.headers.nocov.bed
