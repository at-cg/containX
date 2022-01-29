#!/bin/bash
#SBATCH --clusters=escori
#SBATCH --qos=bigmem
#SBATCH --time=24:00:00
#SBATCH --constraint=haswell
#SBATCH --output=BATCH_OUTPUT
#SBATCH --error=BATCH_OUTPUT

READS=../reads.fasta
OVERLAPS=../minimap2_trial/overlaps.paf
SEQTK=/global/homes/c/cjain7/shared/tools/seqtk-1.3/seqtk
EXE=/project/projectdirs/m3788/tools/minimap2-2.23_x64-linux/minimap2
HAP1=/global/homes/c/cjain7/shared/projects/assembly/chirag/data/genomes/HG002.hifiasm.trio.0.16.1.hap1.fa
HAP2=/global/homes/c/cjain7/shared/projects/assembly/chirag/data/genomes/HG002.hifiasm.trio.0.16.1.hap2.fa
GENOMESIZE=/global/homes/c/cjain7/shared/projects/assembly/chirag/data/genomes/HG002.hifiasm.trio.0.16.1.size

for MINIDENTITY in 100
do
  echo "MINIDENTITY=" ${MINIDENTITY}
  cat ${OVERLAPS} | awk -v minidnty="$MINIDENTITY" '{if ($3 == 0 && $2 == $4 && $2 < $7 && $10*100.0/$11 >= minidnty) print $0}' | cut -f1 > tmp1
  cat ${OVERLAPS} | awk -v minidnty="$MINIDENTITY" '{if ($8 == 0 && $7 == $9 && $7 < $2 && $10*100.0/$11 >= minidnty) print $0}' | cut -f6 > tmp2
  cat tmp1 tmp2 | sort | uniq > contained.${MINIDENTITY}.txt
  cat ../reads.fasta.headers | sort | sed 's/>//g' > reads.fasta.sorted.headers
  comm -23 reads.fasta.sorted.headers contained.${MINIDENTITY}.txt > non-contained.${MINIDENTITY}.txt
  wc -l contained.${MINIDENTITY}.txt
  wc -l non-contained.${MINIDENTITY}.txt
  $SEQTK subseq $READS non-contained.${MINIDENTITY}.txt > non-contained.fasta
  /usr/bin/time $EXE -t 32 -N 50 -cx map-ont $HAP1 non-contained.fasta > mm2.${MINIDENTITY}.paf
  /usr/bin/time $EXE -t 32 -N 50 -cx map-ont $HAP2 non-contained.fasta >> mm2.${MINIDENTITY}.paf
  cat mm2.${MINIDENTITY}.paf | awk '{if ($3 == 0 && $2 == $4 && $2 == $10) print $6"\t"$8"\t"$9}' > tmp3
  cat mm2.${MINIDENTITY}.paf | awk '{if ($8 == 0 && $7 == $9 && $7 == $10) print $6"\t"$8"\t"$9}' >> tmp3
  cat tmp3 | sort -k 1,1 -k2,2n -k3,3nr >  mm2.exactmapped.${MINIDENTITY}.bed
  bedtools genomecov -i mm2.exactmapped.${MINIDENTITY}.bed -g $GENOMESIZE -bga | awk '{if ($4==0) print $0}' > mm2.exactmapped.nocov.${MINIDENTITY}.bed
  bedtools subtract -A -a mm2.exactmapped.nocov.${MINIDENTITY}.bed -b ../reads.fasta.headers.nocov.bed > mm2.exactmapped.nocov.subtracted.${MINIDENTITY}.bed
  bedtools genomecov -i mm2.exactmapped.${MINIDENTITY}.bed -g $GENOMESIZE > mm2.exactmapped.cov.${MINIDENTITY}.hist
  python3 /global/homes/c/cjain7/shared/projects/assembly/chirag/code/useful_scripts/printEdgeBedIntervals.py $GENOMESIZE 25000 > genome_25kbp_ends.bed
  bedtools subtract -A -a mm2.exactmapped.nocov.subtracted.${MINIDENTITY}.bed -b genome_25kbp_ends.bed > mm2.exactmapped.nocov.subtracted.noends.${MINIDENTITY}.bed
done

rm tmp1 tmp2 tmp3 non-contained.fasta  
