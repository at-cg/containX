#ifndef PARAM_CONTAINX_H
#define PARAM_CONTAINX_H

struct algoParams
{
  float min_ovlp_identity;        //minimum overlap identity
  uint32_t min_ovlp_len;          //minimum overlap length (absolute)
  float min_ovlp_ratio;           //minimum overlap ratio relative to longest overlap for a read
  uint32_t minContainedReadLength;//minimum length of contained read
  uint32_t maxContainmentDegree;  //above this, all contained reads will be considered redundant
  uint32_t depthReadLen;          //how far to traverse while collecting minimizers, factor of read length
  uint32_t depthBaseCount;        //how far to traverse while collecting minimizers, absolute base count
                                  //minimum of the above two is used
  uint32_t k;                     //k-mer length (<=16)
  int fuzz;                       //fuzz parameter for transitive reduction (Myers 2005)
  uint32_t iter;                  //user-specified threshold on count of iterations
  uint32_t maxTipLen;             //threshold for tip length (0 means disabled)
  uint32_t threads;               //thread count
  float d;                        //[0,1] sampling density
  float cutoff;                   //what fraction of minimizers must match for redundancy
  bool hpc;                       //parse k-mers in homopolymer compressed space
  bool runHui2016;                //for benchmarking
  bool removeAllContainedReads;
  bool printReadStrings;

  //files for logging output or progress
  std::string gfadumpfilename;
  std::string dumpNonRedudantContainedReads;
  std::string logFileName;
  std::string hetContainedReads;

  void initParams ()
  {
    min_ovlp_identity = 100.0;
    min_ovlp_len = 5000;
    min_ovlp_ratio = 0.0; //disabled
    minContainedReadLength = 10000;
    maxContainmentDegree = UINT32_MAX; //disabled
    hpc = false;
    fuzz = 0;     //for transitive reduction
    iter = 1;     //iterations
    cutoff = 1.0; //[0-1]
    maxTipLen = 0;
    depthReadLen = 2;
    depthBaseCount = UINT32_MAX; //unused by default
    threads = 1;
    k = 16;
    d = 0.25;
    runHui2016 = false;
    removeAllContainedReads = false;
    printReadStrings = false;
  }

  void printParams ()
  {
    std::cerr << "INFO, printParams(), min_ovlp_len = " << min_ovlp_len << "\n";
    std::cerr << "INFO, printParams(), min_ovlp_identity = " << min_ovlp_identity << "\n";
    std::cerr << "INFO, printParams(), min_ovlp_ratio = " << min_ovlp_ratio << "\n";
    std::cerr << "INFO, printParams(), minContainedReadLength = " << minContainedReadLength << "\n";
    std::cerr << "INFO, printParams(), maxContainmentDegree = " << maxContainmentDegree << "\n";
    std::cerr << "INFO, printParams(), depthReadLen = " << depthReadLen << "\n";
    std::cerr << "INFO, printParams(), depthBaseCount = " << depthBaseCount << "\n";
    std::cerr << "INFO, printParams(), k = " << k << "\n";
    std::cerr << "INFO, printParams(), threads = " << threads << "\n";
    std::cerr << "INFO, printParams(), fuzz = " << fuzz << "\n";
    std::cerr << "INFO, printParams(), tip length cutoff = " << maxTipLen << "\n";
    std::cerr << "INFO, printParams(), d (density) = " << d << "\n";
    std::cerr << "INFO, printParams(), minimum fraction of minimizer match [0-1] = " << cutoff << "\n";
    std::cerr << "INFO, printParams(), iter = " << iter << "\n";
    std::cerr << "INFO, printParams(), homopolymer compression = " << std::boolalpha << hpc << "\n";
    if (runHui2016) std::cerr << "INFO, printParams(), runHui2016 = " << std::boolalpha << runHui2016 << "\n";
  }
};

#endif

