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
  bool hpc;                       //parse k-mers in homopolymer compressed space
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
  assert (beg >= 0 && beg < str.length());
  assert (end >= beg && end <= str.length());

  if (end == beg)
    return;

  //'U' is for integer literal of type unsigned int
  uint32_t shift1 = 2*(param.k - 1), mask = (1U<<2*param.k) - 1, kmer = 0;
  float hashkmer;

  for (uint32_t i = beg, l = 0; i < end; ++i) {
    int c = seq_nt4_table[(uint8_t)str[i]];
    if (c < 4) { // not an ambiguous base
      if (param.hpc) {
        int skip_len = 1;
        if (i + 1 < end && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
          for (skip_len = 2; i + skip_len < end; ++skip_len)
            if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
              break;
          i += skip_len - 1;
        }
      }
      if (rev == false)
        kmer = (kmer << 2 | c) & mask;           // forward k-mer
      else
        kmer = (kmer >> 2) | (3U^c) << shift1;   // reverse k-mer
      l++;

      if (l >= k)
      {
        hashkmer = hash32(kmer, mask) * 1.0 / UINT32_MAX;
        if (hashkmer < d) container.emplace_back(kmer);
      }
    }
    else {
      l = 0;
    }
  }
}

/**
 * @param[in] src_vertex  starting vertex
 * @param[in] start_pos   0-based, leftmost position in read string
 */
void dfs_procedure (graphcontainer &g, uint32_t src_vertex, uint32_t start_pos, uint32_t remaining_depth_bases, std::set<uint32_t> &visited_vertices, std::vector<uint32_t> &minimizers, const algoParams &param)
{
  if (remaining_depth_bases == 0) return;
  if (visited_vertices.find(src_vertex) != visited_vertices.end()) return;
  assert (src_vertex < g.vertexCount);

  bool orientation = src_vertex & 1;
  uint32_t src_readid = src_vertex >> 1;
  uint32_t beg = start_pos, end;
  if (remaining_depth_bases > g.readseq[src_readid].length())
  {
    end = g.readseq[src_readid].length();
    remaining_depth_bases = remaining_depth_bases - g.readseq[src_readid].length();
  }
  else
  {
    end = remaining_depth_bases;
    remaining_depth_bases = 0;
  }

  assert (end <= g.readseq[src_readid].length());
  assert (beg < end);

  computeMinimizersFromString(minimizers, g.readseq[src_readid], beg, end, orientation, param);
  visited_vertices.insert (src_vertex);

  if (remaining_depth_bases > 0)
    for (uint32_t j = g.offsets[src_vertex]; j < g.offsets[src_vertex+1]; j++)
      if (g.redundant[j] == false)
        dfs_procedure (g, j, 0, remaining_depth_bases, visited_vertices, mmWalkRead, param);
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
      //confirm if parent reads are still available, otherwise nothing to do
      {
        uint32_t available_parent_count = 0;
        for (uint32_t j = g.containment_offsets[i]; j < g.containment_offsets[i+1]; j++) {
          if (g.redundant[j] == false)
            available_parent_count++;
        }
        if (available_parent_count == 0)
          continue;
      }

      //collect minimizers by starting DFS from contained read
      {
        assert (g.readseq[i].length() > 0);
        uint32_t depth_no_bases = param.depth * g.readseq[i].length();
        std::set<uint32_t> visited_vertices;
        for (uint32_t j = g.offsets[i]; j < g.offsets[i+1]; j++)
          if (g.redundant[j] == false)
            dfs_procedure (g, j, 0, depth_no_bases, visited_vertices, mmWalkRead, param);
      }

      //collect minimizers from parent reads
      {

      }

      //compare minimizer lists and decide
    }
  }
}

#endif
