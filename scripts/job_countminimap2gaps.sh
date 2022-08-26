HAP1=data/genomes/HG002.hifiasm.trio.0.16.1.hap1.fa
HAP2=data/genomes/HG002.hifiasm.trio.0.16.1.hap2.fa

#https://www.biostars.org/p/173963/#174150
GENOMESIZE=data/genomes/HG002.hifiasm.trio.0.16.1.size

cat overlaps.paf | awk '{if ($3 == 0 && $2 == $4 && $2 < $7 && $10*100.0/$11 == 100) print $0}' | cut -f1 > tmp1
cat overlaps.paf | awk '{if ($8 == 0 && $7 == $9 && $7 < $2 && $10*100.0/$11 == 100) print $0}' | cut -f6 > tmp2
cat tmp1 tmp2 | sort | uniq > contained.txt
cat reads.fasta.headers | sort | sed 's/>//g' > reads.fasta.sorted.headers
comm -23 reads.fasta.sorted.headers contained.txt > non-contained.txt
seqtk subseq reads.fasta non-contained.txt > non-contained.fasta
minimap2 -t 32 -N 50 -cx map-ont $HAP1 non-contained.fasta > mm2.paf
minimap2 -t 32 -N 50 -cx map-ont $HAP2 non-contained.fasta >> mm2.paf
cat mm2.paf | awk '{if ($3 == 0 && $2 == $4 && $2 == $10) print $6"\t"$8"\t"$9}' > tmp3
cat mm2.paf | awk '{if ($8 == 0 && $7 == $9 && $7 == $10) print $6"\t"$8"\t"$9}' >> tmp3
cat tmp3 | sort -k 1,1 -k2,2n -k3,3nr >  mm2.exactmapped.bed
bedtools genomecov -i mm2.exactmapped.bed -g $GENOMESIZE -bga | awk '{if ($4==0) print $0}' > mm2.exactmapped.nocov.bed
bedtools subtract -A -a mm2.exactmapped.nocov.bed -b reads.fasta.headers.nocov.bed > mm2.exactmapped.nocov.subtracted.bed
python3 useful_scripts/printEdgeBedIntervals.py $GENOMESIZE 25000 > genome_25kbp_ends.bed
bedtools subtract -A -a mm2.exactmapped.nocov.subtracted.bed -b genome_25kbp_ends.bed > mm2.exactmapped.nocov.subtracted.noends.bed
#combine adjacent intervals that are likely segregrated due to FN calls for contained reads
bedtools merge -d 1000 -i mm2.exactmapped.nocov.subtracted.noends.bed > mm2.exactmapped.nocov.subtracted.noends.merged.bed
rm tmp1 tmp2 tmp3 non-contained.fasta  
