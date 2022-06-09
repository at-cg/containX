#include <iostream>
#include <cassert>
#include <chrono>
#include <zlib.h>

#include "graph.hpp"
#include "algo.hpp"
#include "param.hpp"
#include "ketopt.h"

int main(int argc, char *argv[])
{
  ketopt_t o = KETOPT_INIT;
  int c;
  algoParams param;

  //initialize default values of various parameters
  param.initParams();

  while ((c = ketopt(&o, argc, argv, 1, "cd:D:f:g:hHi:I:l:L:m:M:n:N:p:s:t:T:w:", 0)) >= 0)
  {
    if (c == 'c') param.removeAllContainedReads = true;
    else if (c == 'd') param.gfadumpfilename = o.arg;
    else if (c == 'D') param.gfadumpfilename = o.arg, param.printReadStrings = true;
    else if (c == 'f') param.fuzz = atoi(o.arg);
    else if (c == 'g') param.pacedumpfilename = o.arg;
    else if (c == 'h') param.runHui2016 = true;
    else if (c == 'H') param.hpc = true;
    else if (c == 'i') param.min_ovlp_identity = atof(o.arg);
    else if (c == 'I') param.iter = atoi(o.arg);
    else if (c == 'l') param.min_ovlp_len = atoi(o.arg);
    else if (c == 'L') param.logFileName = o.arg;
    else if (c == 'm') param.cutoff = atof(o.arg);
    else if (c == 'M') param.min_read_len = atoi(o.arg);
    else if (c == 'n') param.dumpNonRedudantContainedReads = o.arg;
    else if (c == 'N') param.dumpNonRedudantReads = o.arg;
    else if (c == 'p') param.hetReads = o.arg;
    else if (c == 's') param.d = atof(o.arg);
    else if (c == 't') param.threads = atoi(o.arg);
    else if (c == 'T') param.maxTipLen = atoi(o.arg);
    else if (c == 'w') param.depthReadLen = atoi(o.arg);
  } 

  //print usage
  if (argc <= o.ind + 1) {
    std::cerr << "Usage: containX [options] <input-reads.fq> <in.paf>\n";
    std::cerr << "Options:\n";
    std::cerr << "  -M NUM      min read length, default " << param.min_read_len << "\n";
    std::cerr << "  -l NUM      min overlap length, default " << param.min_ovlp_len << "\n";
    std::cerr << "  -i NUM      min overlap percentage identity [0.0-100.0], default " << param.min_ovlp_identity << "\n";
    std::cerr << "  -t NUM      thread count, default " << param.threads << "\n";
    std::cerr << "  -I NUM      count of iterations, default " << param.iter << "\n";
    std::cerr << "  -s NUM      sample k-mer with NUM probability, default " << param.d << "\n";
    std::cerr << "  -m NUM      min fraction of minimizer matches for redundant contained reads, default " << param.cutoff << "\n";
    std::cerr << "  -w NUM      walk length cutoff as a factor of read length, default " << param.depthReadLen << "\n";
    std::cerr << "  -H          use homopolymer-compressed k-mer\n";
    std::cerr << "  -c          simply mark all contained reads as redundant\n";
    std::cerr << "  -f NUM      fuzz value during transitive reduction, default " << param.fuzz << " (-1 disables reduction)\n";
    std::cerr << "  -T NUM      threshold for tip length removal, default " << param.maxTipLen << "\n";
    std::cerr << "  -p FILE     list of non-repetitive heterozygous read ids (from hifiasm)\n";
    std::cerr << "  -n FILE     dump read ids of non-redundant contained reads\n";
    std::cerr << "  -N FILE     dump read ids of non-redundant reads\n";
    std::cerr << "  -L FILE     dump algorithm log\n";
    std::cerr << "  -d FILE     dump graph in gfa format without sequences\n";
    std::cerr << "  -D FILE     dump graph in gfa format with sequences\n";
    std::cerr << "  -g FILE     dump undirected graph in PACE format\n";
    return 1;
  }

  assert (param.min_ovlp_identity >= 0.0);
  assert (param.min_ovlp_identity <= 100.0);
  assert (param.cutoff >= 0.0 && param.cutoff <= 1.0);
  assert (param.d > 0.0 && param.d <= 1.0);
  assert (param.depthReadLen > 0);

  param.printParams();
  omp_set_num_threads (param.threads);

  auto tStart = std::chrono::system_clock::now();
  std::cerr << "INFO, main(), started timer\n";

  graphcontainer g;
  ovlgraph_gen (argv[o.ind], argv[o.ind+1], param, g);


  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
  std::cerr << "INFO, main(), graph generation completed after " << wctduration.count() << " seconds\n";

  ovlgraph_simplify (g, param);

  wctduration = (std::chrono::system_clock::now() - tStart);
  std::cerr << "INFO, main(), graph simplification completed after " << wctduration.count() << " seconds\n";

  g.outputGFA (param.gfadumpfilename, param.printReadStrings);
  g.exportPACE (param.pacedumpfilename);
  g.outputNonRedudantContainedReads (param.dumpNonRedudantContainedReads);
  g.outputNonRedudantReads (param.dumpNonRedudantReads);

  //log complete command given by user
  fprintf(stderr, "INFO, %s(), CMD:", __func__);
  for (int i = 0; i < argc; ++i)
    fprintf(stderr, " %s", argv[i]);
  std::cerr << "\n";

  return 0;
}
