#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --clusters=escori
#SBATCH --qos=bigmem
#SBATCH --output=BATCH_OUTPUT
#SBATCH --error=BATCH_OUTPUT

HAP1=/global/homes/c/cjain7/shared/projects/assembly/chirag/data/genomes/HG002.hifiasm.trio.0.16.1.hap1.fa
HAP2=/global/homes/c/cjain7/shared/projects/assembly/chirag/data/genomes/HG002.hifiasm.trio.0.16.1.hap2.fa
SEQTK=/global/homes/c/cjain7/shared/tools/seqtk-1.3/seqtk
READS=/project/projectdirs/m3788/projects/assembly/chirag/data/simulated_reads/ONT/human_diploid_30x_HG02080_HPRC_hifiasm.trio.0.16.1/reads.fasta
EXE=/project/projectdirs/m3788/tools/minimap2-2.23_x64-linux/minimap2
READIDS=../nonRedundantContainedReads.3.txt
GENOMESIZE=/global/homes/c/cjain7/shared/projects/assembly/chirag/data/genomes/HG002.hifiasm.trio.0.16.1.size
GAPS=/project/projectdirs/m3788/projects/assembly/chirag/data/simulated_reads/ONT/human_diploid_30x_HG02080_HPRC_hifiasm.trio.0.16.1/map_mm2_noncontained/mm2.exactmapped.nocov.subtracted.noends.100.bed

/usr/bin/time $SEQTK subseq $READS $READIDS > non-redundant.fasta
/usr/bin/time $EXE -t 32 -N 50 -cx map-ont $HAP1 non-redundant.fasta > mm2.paf
/usr/bin/time $EXE -t 32 -N 50 -cx map-ont $HAP2 non-redundant.fasta >> mm2.paf
cat mm2.paf | awk '{if ($3 == 0 && $2 == $4 && $2 == $10) print $6"\t"$8"\t"$9}' > tmp
cat mm2.paf | awk '{if ($8 == 0 && $7 == $9 && $7 == $10) print $6"\t"$8"\t"$9}' >> tmp
cat tmp | sort -k 1,1 -k2,2n -k3,3nr >  mm2.exactmapped.bed
bedtools subtract -f 1 -a $GAPS -b mm2.exactmapped.bed > unresolved_gaps.bed
