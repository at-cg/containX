#ifndef ASM_ALGO_H
#define ASM_ALGO_H

#include "common.hpp"


struct algoParams
{
  uint32_t maxContainmentDegree;  //above this, all contained reads will be considered redundant
  uint32_t depth;                 //how far to traverse while collecting minimizers, factor of read length
  uint32_t k;                     //k-mer length (<=16)
  float d;                        //[0,1] sampling density
  float cutoff;                   //what fraction of minimizers must match for redundancy
};

/**
 * @param[out]  container     push new minimizers here
 * @param[in]   str           original read string
 * @param[in]   beg           begin offset (inclusive), leftmost offset in string
 * @param[in]   end           end offset (exclusive), rightmost offset in string
 * @param[in]   rev           consider reverse complement of string if true
 */
void computeMinimizersFromString(std::vector<uint32_t> &container, const std::string &str, uint32_t beg, uint32_t end, bool rev, const algoParams &param)
{
  assert (param.k >= 10 && param.k <= 16);
  assert (param.d > 0 && param.d <= 1.0);

  //TODO: complete this function
}

void identifyRedundantReads(graphcontainer &g, const algoParams &param)
{
  g.removeContainedReadsAboveDegree (param.maxContainmentDegree);
  g.index();

  std::vector<uint32_t> mmWalkRead; //walk along same orientation as the read, and collect minimizers
  std::vector<uint32_t> mmWalkParentReads; //collect minimizers by walking from parent reads containing the current read
  //TODO: should we also walk in opposite direction?

  //TODO: reconsider if we should process contained reads in a particular order
  for (uint32_t i = 0; i < g.readCount; i++)
  {
    if (g.contained[i] == true && g.redundant[i] == false)
    {
    }
  }
}

#endif
