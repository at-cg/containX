#include <iostream>
#include <cassert>
#include <zlib.h>
#include "paf.hpp"
#include "ketopt.h"

int main(int argc, char *argv[])
{
  ketopt_t o = KETOPT_INIT;
  float min_ovlp_identity = 0.0; //[0-100]
  int c;

  while ((c = ketopt(&o, argc, argv, 1, "o:", 0)) >= 0)
  {
    if (c == 'o') min_ovlp_identity = atof(o.arg);
  }

  //print usage
  if (argc == o.ind) {
    std::cerr << "Usage: containX [options] <in.paf>\n";
    std::cerr << "Options:\n";
    std::cerr << "  -o NUM      min overlap percentage identity [0.0-100.0], default " << min_ovlp_identity << "\n";
    return 1;
  }

  assert (min_ovlp_identity >= 0.0);
  assert (min_ovlp_identity <= 100.0);

  //read input paf file
  paf_rec_t r;
  paf_file_t *fp = paf_open(argv[o.ind]);
  if (!fp) {
    std::cerr << "could not open PAF file\n";
    exit(1);
  }

  uint64_t tot = 0;
  while (paf_read(fp, &r) >= 0) {
    tot++;
  }

  std::cerr << "read " << tot << " paf records\n";

  paf_close(fp);

  return 0;
}
