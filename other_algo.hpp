#ifndef ASM_OTHER_ALGO_H
#define ASM_OTHER_ALGO_H

/**
 * procedure to decide whether to retain or ignore contained read from Hui et al. 2016
 * "Overlap-Based Genome Assembly from Variable-Length Reads"
 */

#include <set>
#include "common.hpp"
#include "graph.hpp"
#include "param.hpp"
#include <omp.h>

//check consistency (matches) of characters followed by contained read
bool checkConsistencyForward (const graphcontainer &g, const containmentTuple &c1, const containmentTuple &c2)
{
  uint32_t l1, r1, readLen1 = g.readseq[c1.dst].length();
  uint32_t l2, r2, readLen2 = g.readseq[c2.dst].length();
  uint32_t beg1, beg2;

  if (c1.rev) { //contained read had matched parent read in rev orientation, see prefix of parent read
    l1 = 0, r1 = c1.dst_start_offset;
    if (r1 == 0) return true;
    beg1 = r1 - 1;
  }
  else { //contained read had matched parent read in fwd orientation, see suffix of parent read
    l1 = c1.dst_end_offset, r1 = readLen1;
    if (l1 == readLen1) return true;
    beg1 = l1;
  }

  if (c2.rev) {
    l2 = 0, r2 = c2.dst_start_offset;
    if (r2 == 0) return true;
    beg2 = r2 - 1;
  }
  else {
    l2 = c2.dst_end_offset, r2 = readLen2;
    if (l2 == readLen2) return true;
    beg2 = l2;
  }

  assert (r1 > l1 && r2 > l2);
  uint32_t charactersToMatch = std::min (r1-l1, r2-l2);

  for (uint32_t i = 0; i < charactersToMatch; i++)
  {
    int s1 = seq_nt4_table[(uint8_t) g.readseq[c1.dst][beg1]];
    if (c1.rev)  s1 = (3U^s1);

    int s2 = seq_nt4_table[(uint8_t) g.readseq[c2.dst][beg2]];
    if (c2.rev)  s2 = (3U^s2);

    if (s1 != s2) return false;

    //shift offsets
    if (c1.rev) beg1--; else beg1++;
    if (c2.rev) beg2--; else beg2++;
  }

  return true;
}


//check consistency (matches) of characters prior to contained read
//symmetric to checkConsistencyForward()
bool checkConsistencyBackward (const graphcontainer &g, const containmentTuple &c1, const containmentTuple &c2)
{
  uint32_t l1, r1, readLen1 = g.readseq[c1.dst].length();
  uint32_t l2, r2, readLen2 = g.readseq[c2.dst].length();
  uint32_t beg1, beg2;

  if (c1.rev) {
    l1 = c1.dst_end_offset, r1 = readLen1;
    if (l1 == readLen1) return true;
    beg1 = l1;
  }
  else {
    l1 = 0, r1 = c1.dst_start_offset;
    if (r1 == 0) return true;
    beg1 = r1 - 1;
  }

  if (c2.rev) {
    l2 = c2.dst_end_offset, r2 = readLen2;
    if (l2 == readLen2) return true;
    beg2 = l2;
  }
  else {
    l2 = 0, r2 = c2.dst_start_offset;
    if (r2 == 0) return true;
    beg2 = r2 - 1;
  }

  assert (r1 > l1 && r2 > l2);
  uint32_t charactersToMatch = std::min (r1-l1, r2-l2);

  for (uint32_t i = 0; i < charactersToMatch; i++)
  {
    int s1 = seq_nt4_table[(uint8_t) g.readseq[c1.dst][beg1] ];
    if (c1.rev)  s1 = (3U^s1);

    int s2 = seq_nt4_table[(uint8_t) g.readseq[c2.dst][beg2] ];
    if (c2.rev)  s2 = (3U^s2);

    if (s1 != s2) return false;

    //shift offsets
    if (c1.rev) beg1++; else beg1--;
    if (c2.rev) beg2++; else beg2--;
  }

  return true;
}

//check characters in parent reads on both sides of contained read
bool checkConsistency (const graphcontainer &g, const containmentTuple &c1, const containmentTuple &c2)
{
  return checkConsistencyForward (g, c1, c2) && checkConsistencyBackward (g, c1, c2);
}

void identifyRedundantReadsHuiEtAl(graphcontainer &g, const algoParams &param, std::ofstream& log)
{
#pragma omp for schedule(static, 1)
  for (uint32_t i = 0; i < g.readCount; i++)
  {
    if (g.contained[i] == true && g.deletedReads[i] == false)
    {
      bool consistent = true;
      //do pairwise comparisons among parent reads
      for (uint32_t j = g.containment_offsets[i]; j < g.containment_offsets[i+1]; j++)
      {
        for (uint32_t k = j+1; k < g.containment_offsets[i+1]; k++)
        {
          consistent = checkConsistency (g, g.containments[j], g.containments[k]);
          if (!consistent)
            j = k = g.containment_offsets[i+1]; //exit the two inner loops
        }
      }

      if (!consistent)
      {
        g.deletedReads[i] = true;
        if (!param.logFileName.empty()) {
#pragma omp critical
          log << g.umap_inverse[i] << "\tidentifyRedundantReadsHuiEtAl()\tREDUNDANT=T\n";
        }
      }
      else
      {
        if (!param.logFileName.empty()) {
#pragma omp critical
          log << g.umap_inverse[i] << "\tidentifyRedundantReadsHuiEtAl()\tREDUNDANT=F\n";
        }

      }
    }
  }
  std::cerr << "INFO, identifyRedundantReadsHuiEtAl() finished\n";
}

#endif
