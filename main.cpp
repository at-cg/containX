#include <iostream>
#include <cassert>
#include <zlib.h>

#include "graph.hpp"
#include "algo.hpp"
#include "param.hpp"
#include "ketopt.h"

int main(int argc, char *argv[])
{
  ketopt_t o = KETOPT_INIT;
  float min_ovlp_identity = 100.0; //[0-100]
  int min_ovlp_len = 5000;
  int c;
  bool printReadStrings = false;
  bool removeAllContainedReads = false;
  algoParams param;
  param.hpc = false;
  param.fuzz = 100; //disable transitive reduction of edges by default
  param.max_iter = 5; //upper bound for graph simplification iterations
  param.cutoff = 1.0; //[0-1]
  param.maxTipLen = 3;
  param.depthReadLen = 5, param.depthBaseCount = 100000;
  param.threads = 1;

  while ((c = ketopt(&o, argc, argv, 1, "cd:D:Hi:I:l:L:m:n:t:T:w:W:", 0)) >= 0)
  {
    if (c == 'c') removeAllContainedReads = true;
    else if (c == 'd') param.gfadumpfilename = o.arg;
    else if (c == 'D') param.gfadumpfilename = o.arg, printReadStrings = true;
    else if (c == 'f') param.fuzz = atoi(o.arg);
    else if (c == 'H') param.hpc = true;
    else if (c == 'i') min_ovlp_identity = atof(o.arg);
    else if (c == 'I') param.max_iter = atoi(o.arg);
    else if (c == 'l') min_ovlp_len = atoi(o.arg);
    else if (c == 'm') param.cutoff = atof(o.arg);
    else if (c == 'n') param.dumpNonRedudantContainedReads = o.arg;
    else if (c == 'L') param.logFileName = o.arg;
    else if (c == 't') param.threads = atoi(o.arg);
    else if (c == 'T') param.maxTipLen = atoi(o.arg);
    else if (c == 'w') param.depthReadLen = atoi(o.arg);
    else if (c == 'W') param.depthBaseCount = atoi(o.arg);
  }

  //print usage
  if (argc <= o.ind + 1) {
    std::cerr << "Usage: containX [options] <input-reads.fq> <in.paf>\n";
    std::cerr << "Options:\n";
    std::cerr << "  -l NUM      min overlap length, default " << min_ovlp_len << "\n";
    std::cerr << "  -i NUM      min overlap percentage identity [0.0-100.0], default " << min_ovlp_identity << "\n";
    std::cerr << "  -t NUM      thread count, default " << param.threads << "\n";
    std::cerr << "  -I NUM      max count of iterations, default " << param.max_iter << "\n";
    std::cerr << "  -m NUM      min fraction of minimizer matches for redundant contained reads, default " << param.cutoff << "\n";
    std::cerr << "  -w NUM      walk length cutoff as a factor of read length, default " << param.depthReadLen << "\n";
    std::cerr << "  -W NUM      walk length cutoff in terms of absoute base count, default " << param.depthBaseCount << "\n";
    std::cerr << "  -H          use homopolymer-compressed k-mer\n";
    std::cerr << "  -c          simply mark all contained reads as redundant and remove\n";
    std::cerr << "  -f NUM      fuzz value during transitive reduction, default " << param.fuzz << "\n";
    std::cerr << "  -T NUM      threshold for tip length removal, default " << param.maxTipLen << ", set 0 to disable\n";
    std::cerr << "  -n FILE     dump read ids of non-redundant contained reads\n";
    std::cerr << "  -L FILE     dump algorithm log\n";
    std::cerr << "  -d FILE     dump graph in gfa format without sequences\n";
    std::cerr << "  -D FILE     dump graph in gfa format with sequences\n";
    return 1;
  }

  assert (min_ovlp_identity >= 0.0);
  assert (min_ovlp_identity <= 100.0);
  assert (param.cutoff >= 0.0 && param.cutoff <= 1.0);
  assert (param.depthReadLen > 0);
  assert (param.depthBaseCount > 0);

  //set parameters
  param.maxContainmentDegree = 5; //warning: this would need adjustment with ploidy
  param.k = 16;
  param.d = 1.0/param.k;
  param.printParams();
  omp_set_num_threads (param.threads);

  graphcontainer g;
  ovlgraph_gen (argv[o.ind], argv[o.ind+1], min_ovlp_identity, min_ovlp_len, g);

#ifdef VERBOSE
  printContainmentDegreeDistribution (g, "ContainmentDegree.beforeSimplify.txt");
  printDegreeDistribution (g, "Degree.beforeSimplify.txt");
  printEdgesDOTFormat (g, "edges.beforeSimplify.DOT");
#endif

  ovlgraph_simplify (removeAllContainedReads, g, param);


#ifdef VERBOSE
  printContainmentDegreeDistribution (g, "ContainmentDegree.afterSimplify.txt");
  printDegreeDistribution (g, "Degree.afterSimplify.txt");
  printEdgesDOTFormat (g, "edges.afterSimplify.DOT");
#endif

  g.outputGFA (param.gfadumpfilename, printReadStrings);
  g.outputNonRedudantContainedReads (param.dumpNonRedudantContainedReads);


  //log complete command given by user
  fprintf(stderr, "INFO, %s(), CMD:", __func__);
  for (int i = 0; i < argc; ++i)
    fprintf(stderr, " %s", argv[i]);
  std::cerr << "\n";

  return 0;
}
