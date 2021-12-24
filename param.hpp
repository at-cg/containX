#ifndef PARAM_CONTAINX_H
#define PARAM_CONTAINX_H

struct algoParams
{
  uint32_t maxContainmentDegree;  //above this, all contained reads will be considered redundant
  uint32_t depthReadLen;          //how far to traverse while collecting minimizers, factor of read length
  uint32_t depthBaseCount;        //how far to traverse while collecting minimizers, absolute base count
                                  //minimum of the above two is used
  uint32_t k;                     //k-mer length (<=16)
  uint32_t fuzz;                  //fuzz parameter for transitive reduction (Myers 2005)
  uint32_t max_iter;              //user-specified threshold on count of iterations
  uint32_t maxTipLen;             //threshold for tip length (0 means disabled)
  float d;                        //[0,1] sampling density
  float cutoff;                   //what fraction of minimizers must match for redundancy
  bool hpc;                       //parse k-mers in homopolymer compressed space

  //files for logging output or progress
  std::string gfadumpfilename; 
  std::string dumpNonRedudantContainedReads;
  std::string logFileName;

  void printParams ()
  {
    std::cerr << "INFO, printParams(), maxContainmentDegree = " << maxContainmentDegree << "\n";
    std::cerr << "INFO, printParams(), depthReadLen = " << depthReadLen << "\n";
    std::cerr << "INFO, printParams(), depthBaseCount = " << depthBaseCount << "\n";
    std::cerr << "INFO, printParams(), k = " << k << "\n";
    std::cerr << "INFO, printParams(), fuzz = " << fuzz << "\n";
    std::cerr << "INFO, printParams(), tip length cutoff = " << maxTipLen << "\n";
    std::cerr << "INFO, printParams(), d (density) = " << d << "\n";
    std::cerr << "INFO, printParams(), cutoff = " << cutoff << "\n";
    std::cerr << "INFO, printParams(), max_iter = " << max_iter << "\n";
    std::cerr << "INFO, printParams(), hpc = " << std::boolalpha << hpc << "\n";
  }
};

#endif

