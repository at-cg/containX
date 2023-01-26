READS=sample.fa
# either have minimap2 executable in your PATH or provide absolute path to minimap2 executable
minimap2 -t 32 -w 101 -k 27 -g 500 -B 8 -O 8,48 -E 4,2 -cx ava-ont $READS $READS > overlaps.paf
../../containX -l 10 -n nonRedundantContainedReads.txt -N nonRedundantReads.txt -d graph.gfa $READS overlaps.paf #remove selected contained reads
