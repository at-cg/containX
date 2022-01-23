# containX
```
Usage: containX [options] <input-reads.fq> <in.paf>
Options:
  -l NUM      min overlap length, default 5000
  -i NUM      min overlap percentage identity [0.0-100.0], default 100
  -t NUM      thread count, default 1
  -s NUM      sample k-mer with NUM probability, default 0.25
  -m NUM      min fraction of minimizer matches for redundant contained reads, default 1
  -w NUM      walk length cutoff as a factor of read length, default 2
  -H          use homopolymer-compressed k-mer
  -c          simply mark all contained reads as redundant
  -f NUM      fuzz value during transitive reduction, default 0 (-1 disables reduction)
  -T NUM      threshold for tip length removal, default 0
  -p FILE     list of read ids to ignore
  -n FILE     dump read ids of non-redundant contained reads
  -N FILE     dump read ids of non-redundant reads
  -L FILE     dump algorithm log
  -d FILE     dump graph in gfa format without sequences
  -D FILE     dump graph in gfa format with sequences
```
