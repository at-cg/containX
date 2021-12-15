READS=sample.fa

/usr/bin/time $MINIMAP2 -t 32 -w 101 -k 27 -g 500 -B 8 -O 8,48 -E 4,2 -cx ava-ont $READS $READS > overlaps.paf
/usr/bin/time $CONTAINX -l 10 -i 100 -d dump.gfa $READS overlaps.paf
/usr/bin/time $CONTAINX -l 10 -i 100 -t 100 -d dump2.gfa $READS overlaps.paf
/usr/bin/time $CONTAINX -l 10 -i 100 -t 100 -c -d dump3.gfa $READS overlaps.paf
