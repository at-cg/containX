#ifndef PARAM_CONTAINX_H
#define PARAM_CONTAINX_H

struct algoParams
{
  float min_ovlp_identity;        //minimum overlap identity
  uint32_t min_ovlp_len;          //minimum overlap length (absolute)
  uint32_t min_read_len;          //minimum read length
  uint32_t depthReadLen;          //how far to traverse while collecting minimizers, factor of read length
  uint32_t k;                     //k-mer length (<=16)
  int fuzz;                       //fuzz parameter for transitive reduction (Myers 2005)
  uint32_t iter;                  //count of algo iterations for benchmarking
  uint32_t maxTipLen;             //threshold for tip length (set 0 to remove unconnected vertices only)
  uint32_t threads;               //thread count
  float d;                        //[0,1] sampling density
  float cutoff;                   //what fraction of minimizers must match for redundancy
  bool hpc;                       //parse k-mers in homopolymer compressed space
  bool runHui2016;                //for benchmarking
  bool removeAllContainedReads;   //for benchmarking
  bool printReadStrings;          //print read strings in gfa file

  //files for logging output or progress
  std::string gfadumpfilename;
  std::string pacedumpfilename;
  std::string dumpNonRedudantContainedReads;
  std::string dumpNonRedudantReads;
  std::string logFileName;
  std::string hetReads;           //non-repetitive heterozygous reads predicted using hifiasm

  void initParams ()
  {
    min_ovlp_identity = 100.0;
    min_ovlp_len = 5000;
    min_read_len = 5000;
    hpc = false;
    fuzz = 0;     //for transitive reduction
    iter = 1;     //iterations
    cutoff = 1.0; //[0-1]
    maxTipLen = 0; //0 removes unconnected vertices only
    depthReadLen = 2;
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
    std::cerr << "INFO, printParams(), min_read_len = " << min_read_len << "\n";
    std::cerr << "INFO, printParams(), min_ovlp_identity = " << min_ovlp_identity << "\n";
    std::cerr << "INFO, printParams(), depthReadLen = " << depthReadLen << "\n";
    std::cerr << "INFO, printParams(), k = " << k << "\n";
    std::cerr << "INFO, printParams(), threads = " << threads << "\n";
    std::cerr << "INFO, printParams(), fuzz = " << fuzz << "\n";
    std::cerr << "INFO, printParams(), tip length cutoff = " << maxTipLen << "\n";
    std::cerr << "INFO, printParams(), d (density) = " << d << "\n";
    std::cerr << "INFO, printParams(), minimum fraction of minimizer match [0-1] = " << cutoff << "\n";
    std::cerr << "INFO, printParams(), iter = " << iter << "\n";
    std::cerr << "INFO, printParams(), homopolymer compression = " << std::boolalpha << hpc << "\n";
    if (runHui2016) std::cerr << "INFO, printParams(), runHui2016 = " << std::boolalpha << runHui2016 << "\n";
    if (removeAllContainedReads) std::cerr << "INFO, printParams(), removeAllContainedReads = " << std::boolalpha << removeAllContainedReads << "\n";
  }
};

#endif

