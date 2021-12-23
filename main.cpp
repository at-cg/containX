#include <iostream>
#include <cassert>
#include <zlib.h>

#include "graph.hpp"
#include "algo.hpp"
#include "ketopt.h"

int main(int argc, char *argv[])
{
  ketopt_t o = KETOPT_INIT;
  float min_ovlp_identity = 0.0; //[0-100]
  int min_ovlp_len = 0;
  int c;
  std::string gfadumpfilename;
  std::string dumpNonRedudantContainedReads;
  bool printReadStrings = true;
  bool removeAllContainedReads = false;
  algoParams param;
  param.hpc = false;
  param.fuzz = 100; //disable transitive reduction of edges by default
  param.max_iter = 1; //upper bound for graph simplification iterations
  param.cutoff = 1.0; //[0-1]
  param.maxTipLen = 3;

  while ((c = ketopt(&o, argc, argv, 1, "cd:D:Hi:I:l:m:n:t:T:", 0)) >= 0)
  {
    if (c == 'c') removeAllContainedReads = true;
    else if (c == 'd') gfadumpfilename = o.arg, printReadStrings = false;
    else if (c == 'D') gfadumpfilename = o.arg;
    else if (c == 'H') param.hpc = true;
    else if (c == 'i') min_ovlp_identity = atof(o.arg);
    else if (c == 'I') param.max_iter = atoi(o.arg);
    else if (c == 'l') min_ovlp_len = atoi(o.arg);
    else if (c == 'm') param.cutoff = atof(o.arg);
    else if (c == 'n') dumpNonRedudantContainedReads = o.arg;
    else if (c == 't') param.fuzz = atoi(o.arg);
    else if (c == 'T') param.maxTipLen = atoi(o.arg);
  }

  //print usage
  if (argc <= o.ind + 1) {
    std::cerr << "Usage: containX [options] <input-reads.fq> <in.paf>\n";
    std::cerr << "Options:\n";
    std::cerr << "  -l NUM      min overlap length, default " << min_ovlp_len << "\n";
    std::cerr << "  -i NUM      min overlap percentage identity [0.0-100.0], default " << min_ovlp_identity << "\n";
    std::cerr << "  -I NUM      max count of iterations, default " << param.max_iter << "\n";
    std::cerr << "  -m NUM      min fraction of minimizer matches for redundant contained reads, default " << param.cutoff << "\n";
    std::cerr << "  -H          use homopolymer-compressed k-mer\n";
    std::cerr << "  -p          mark contained reads at ends as redundant\n";
    std::cerr << "  -c          simply mark all contained reads as redundant and remove\n";
    std::cerr << "  -t NUM      fuzz value during transitive reduction, default " << param.fuzz << "\n";
    std::cerr << "  -T NUM      threshold for tip length removal, default " << param.maxTipLen << ", set 0 to disable\n";
    std::cerr << "  -n FILE     dump read ids of non-redundant contained reads\n";
    std::cerr << "  -d FILE     dump graph in gfa format without sequences\n";
    std::cerr << "  -D FILE     dump graph in gfa format with sequences\n";
    return 1;
  }

  assert (min_ovlp_identity >= 0.0);
  assert (min_ovlp_identity <= 100.0);
  assert (param.cutoff >= 0.0 && param.cutoff <= 1.0);

  //set parameters
  param.maxContainmentDegree = 5; //warning: this would need adjustment with ploidy
  param.depth = 5;
  param.k = 16;
  param.d = 1.0/param.k;
  param.printParams();

  graphcontainer g;
  ovlgraph_gen (argv[o.ind], argv[o.ind+1], min_ovlp_identity, min_ovlp_len, g);

#ifdef VERBOSE
  printContainmentDegreeDistribution (g, "ContainmentDegree.beforeSimplify.txt");
  printDegreeDistribution (g, "Degree.beforeSimplify.txt");
  printEdgesDOTFormat (g, "edges.beforeSimplify.DOT");
#endif

  ovlgraph_simplify (removeAllContainedReads, g, param, "log_simplify.txt");


#ifdef VERBOSE
  printContainmentDegreeDistribution (g, "ContainmentDegree.afterSimplify.txt");
  printDegreeDistribution (g, "Degree.afterSimplify.txt");
  printEdgesDOTFormat (g, "edges.afterSimplify.DOT");
#endif

  if (!gfadumpfilename.empty())
    g.outputGFA (gfadumpfilename, printReadStrings);

  if (!dumpNonRedudantContainedReads.empty())
    g.outputNonRedudantContainedReads (dumpNonRedudantContainedReads);


  //log complete command given by user
  fprintf(stderr, "INFO, %s(), CMD:", __func__);
  for (int i = 0; i < argc; ++i)
    fprintf(stderr, " %s", argv[i]);
  std::cerr << "\n";

  return 0;
}
