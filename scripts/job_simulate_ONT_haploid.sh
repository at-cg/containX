#Purpose: Simulate reads from a haploid genome using seqrequester 

GENOME=data/genomes/chm13.draft_v1.1.fasta
COVERAGE=20
DIST=ONT/HG02080_HPRC/HG02080_2.fastq.hist.txt
SIZE_SUM=`grep -v ">" $GENOME | wc | awk '{print $3-$1}'`

#https://www.biostars.org/p/173963/#174150
GENOMESIZE=data/genomes/chm13.draft_v1.1.fasta.size

seqrequester simulate -truncate -genome $GENOME -genomesize $SIZE_SUM -coverage $COVERAGE -distribution $DIST > reads.fasta
grep ">" reads.fasta > reads.fasta.headers
cat reads.fasta.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > reads.fasta.headers.bed
python countcontainedreads.py reads.fasta.headers.bed
bedtools genomecov -i reads.fasta.headers.bed -g $GENOMESIZE -bga | awk '{if ($4==0) print $0}' > reads.fasta.headers.nocov.bed
