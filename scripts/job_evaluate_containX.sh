#Purpose: Evaluate containX and other methods for selecting useful contained reads

HAP1=data/genomes/HG002.hifiasm.trio.0.16.1.hap1.fa
HAP2=data/genomes/HG002.hifiasm.trio.0.16.1.hap2.fa
READS=human_diploid_30x_HG02080_HPRC_hifiasm.trio.0.16.1/reads.fasta
READIDS=nonRedundantContainedReads.txt
GAPS=map_mm2_noncontained/mm2.exactmapped.nocov.subtracted.noends.100.bed

seqtk subseq $READS $READIDS > non-redundant.fasta
minimap2 -t 32 -N 50 -cx map-ont $HAP1 non-redundant.fasta > mm2.paf
minimap2 -t 32 -N 50 -cx map-ont $HAP2 non-redundant.fasta >> mm2.paf
cat mm2.paf | awk '{if ($3 == 0 && $2 == $4 && $2 == $10) print $6"\t"$8"\t"$9}' > tmp
cat mm2.paf | awk '{if ($8 == 0 && $7 == $9 && $7 == $10) print $6"\t"$8"\t"$9}' >> tmp
cat tmp | sort -k 1,1 -k2,2n -k3,3nr >  mm2.exactmapped.bed
bedtools subtract -f 1 -a $GAPS -b mm2.exactmapped.bed > unresolved_gaps.bed
