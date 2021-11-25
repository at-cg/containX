#include <iostream>
#include <cassert>
#include <zlib.h>
#include "graph.hpp"
#include "ketopt.h"

int main(int argc, char *argv[])
{
  ketopt_t o = KETOPT_INIT;
  float min_ovlp_identity = 0.0; //[0-100]
  int min_ovlp_len = 0; 
  int c;

  while ((c = ketopt(&o, argc, argv, 1, "l:i:", 0)) >= 0)
  {
    if (c == 'i') min_ovlp_identity = atof(o.arg);
    else if (c == 'l') min_ovlp_len = atoi(o.arg);
  }

  //print usage
  if (argc <= o.ind + 1) {
    std::cerr << "Usage: containX [options] <input-reads.fq> <in.paf>\n";
    std::cerr << "Options:\n";
    std::cerr << "  -l NUM      min overlap length, default " << min_ovlp_len << "\n";
    std::cerr << "  -o NUM      min overlap percentage identity [0.0-100.0], default " << min_ovlp_identity << "\n";
    return 1;
  }

  assert (min_ovlp_identity >= 0.0);
  assert (min_ovlp_identity <= 100.0);

  graphcontainer g; 
  ovlgraph_gen (argv[o.ind], argv[o.ind+1], min_ovlp_identity, min_ovlp_len, g);

  return 0;
}
