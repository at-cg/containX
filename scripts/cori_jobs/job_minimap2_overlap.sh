#!/bin/bash
#SBATCH --clusters=escori
#SBATCH --qos=bigmem
#SBATCH --time=36:00:00
#SBATCH --constraint=haswell
#SBATCH --output=BATCH_OUTPUT
#SBATCH --error=BATCH_OUTPUT

EXE=/project/projectdirs/m3788/tools/minimap2-2.23_x64-linux/minimap2
READS=../reads.fasta
/usr/bin/time $EXE -t 32 -w 101 -k 27 -g 500 -B 8 -O 8,48 -E 4,2 -cx ava-ont $READS $READS > overlaps.paf
