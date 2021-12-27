# containX
```
Usage: containX [options] <input-reads.fq> <in.paf>
Options:
  -l NUM      min overlap length, default 5000
  -i NUM      min overlap percentage identity [0.0-100.0], default 100
  -t NUM      thread count, default 1
  -I NUM      count of iterations, default 2
  -s NUM      sample k-mer with NUM probability, default 0.25
  -m NUM      min fraction of minimizer matches for redundant contained reads, default 1
  -w NUM      walk length cutoff as a factor of read length, default 2
  -W NUM      walk length cutoff in terms of absoute base count, default 4294967295
  -H          use homopolymer-compressed k-mer
  -c          simply mark all contained reads as redundant
  -C NUM      mark reads contained in >NUM reads as redundant, default 4294967295
  -f NUM      fuzz value during transitive reduction, default 0
  -T NUM      threshold for tip length removal, default 3, set 0 to disable
  -n FILE     dump read ids of non-redundant contained reads
  -L FILE     dump algorithm log
  -d FILE     dump graph in gfa format without sequences
  -D FILE     dump graph in gfa format with sequences
```
