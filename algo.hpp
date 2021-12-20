#ifndef ASM_ALGO_H
#define ASM_ALGO_H

#include <set>
#include "common.hpp"
#include "graph.hpp"

struct algoParams
{
  uint32_t maxContainmentDegree;  //above this, all contained reads will be considered redundant
  uint32_t depth;                 //how far to traverse while collecting minimizers, factor of read length
  uint32_t k;                     //k-mer length (<=16)
  uint32_t fuzz;                  //fuzz parameter for transitive reduction (Myers 2005)
  float d;                        //[0,1] sampling density
  float cutoff;                   //what fraction of minimizers must match for redundancy
  bool hpc;                       //parse k-mers in homopolymer compressed space

  void printParams ()
  {
    std::cerr << "INFO, printParams(), maxContainmentDegree = " << maxContainmentDegree << "\n";
    std::cerr << "INFO, printParams(), depth = " << depth << "\n";
    std::cerr << "INFO, printParams(), k = " << k << "\n";
    std::cerr << "INFO, printParams(), fuzz = " << fuzz << "\n";
    std::cerr << "INFO, printParams(), d (density) = " << d << "\n";
    std::cerr << "INFO, printParams(), cutoff = " << cutoff << "\n";
    std::cerr << "INFO, printParams(), hpc = " << std::boolalpha << hpc << "\n";
  }
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
  assert (beg >= 0 && beg <= str.length());
  assert (end >= beg && end <= str.length());

  if (end == beg)
    return;

  //'U' is for integer literal of type unsigned int
  uint32_t shift1 = 2*(param.k - 1), mask = (1U<<2*param.k) - 1, kmer = 0;
  if (param.k == 16) mask = UINT32_MAX; //corner case
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

      if (l >= param.k)
      {
        hashkmer = hash32(kmer, mask) * 1.0 / UINT32_MAX;
        if (hashkmer < param.d) container.emplace_back(kmer);
      }
    }
    else {
      l = 0;
    }
  }
}

/**
 * @param[in] src_vertex              starting vertex
 * @param[in] beg                     0-based, leftmost position in the current read string
 * @param[in] end                     1-based, end string processing here in the current read string
 * @param[in] remaining_depth_bases   total count of bases to process during DFS
 */
uint32_t dfs_procedure (graphcontainer &g, uint32_t src_vertex, uint32_t beg, uint32_t end, uint32_t remaining_depth_bases, std::set<uint32_t> &visited_vertices, std::vector<uint32_t> &minimizers, const algoParams &param)
{
  if (visited_vertices.find(src_vertex) != visited_vertices.end()) return 0U; //visited already
  uint32_t src_readid = src_vertex >> 1;

  assert (src_vertex < g.vertexCount);
  assert (beg <= end);
  assert (end <= g.readseq[src_readid].length());
  assert (remaining_depth_bases >= end - beg); //this should be ensured before calling DFS

  uint32_t bases_processed = 0;
  visited_vertices.insert (src_vertex);
  bool rev = src_vertex & 1; //orientation

#ifdef VERBOSE
  std::cerr << "INFO, dfs_procedure(), computing minimizers from read " << src_readid << ", vertex = " << src_vertex << ", offsets = [" << beg << "," << end << "), orientation = '" << "+-"[rev] << "'\n";
#endif

  computeMinimizersFromString(minimizers, g.readseq[src_readid], beg, end, rev, param);
  bases_processed = end - beg;
  remaining_depth_bases = remaining_depth_bases - bases_processed;

  //move to neighbor vertices
  if (remaining_depth_bases > 0) {
    for (uint32_t j = g.offsets[src_vertex]; j < g.offsets[src_vertex+1]; j++) {
      assert (g.edges[j].src == src_vertex);
      uint32_t adjVertexId = g.edges[j].dst;
      uint32_t adjReadId = adjVertexId >> 1;

      if (g.redundant[adjReadId] == false)
      {
        //compute begin and end offsets for next read/vertex
        uint32_t nextReadLen = g.readseq[adjReadId].length();
        rev = adjVertexId & 1;

        if (rev == false)
        {
          beg = g.edges[j].ov_dst; //skip overlapping portion
          end = std::min (nextReadLen, beg + remaining_depth_bases);
        }
        else
        {
          end = nextReadLen - g.edges[j].ov_dst; //skip overlapping portion
          beg = end - std::min (end, remaining_depth_bases);
        }
        bases_processed += dfs_procedure (g, adjVertexId, beg, end, remaining_depth_bases, visited_vertices, minimizers, param);
      }
    }
  }
  return bases_processed;
}

void identifyRedundantReads(graphcontainer &g, const algoParams &param, std::ofstream& log)
{
  std::vector<uint32_t> mmWalkRead; //walk along same orientation as the read, and collect minimizers
  std::vector<uint32_t> mmWalkParentReads; //collect minimizers by walking from parent reads containing the current read
  std::set<uint32_t> visited_vertices;

  //TODO: reconsider if we should process contained reads in a particular order
  for (uint32_t i = 0; i < g.readCount; i++)
  {
    assert (g.readseq[i].length() > 0);
    uint32_t available_parent_count = 0;
    mmWalkRead.clear();
    mmWalkParentReads.clear();

    if (g.contained[i] == true && g.redundant[i] == false)
    {
      //confirm if parent reads are still available, otherwise nothing to do
      {
        for (uint32_t j = g.containment_offsets[i]; j < g.containment_offsets[i+1]; j++)
          if (g.redundant[g.containments[j].dst] == false)
            available_parent_count++;
        if (available_parent_count == 0) continue;
      }

      //user-specified depth during DFS in terms of count of bases
      uint32_t depth_bases = param.depth * g.readseq[i].length();

      //collect minimizers by starting DFS from contained read
      {
        visited_vertices.clear();
        uint32_t bases_processed;
        uint32_t vertexId = i << 1 | 0; //forward orientation
        //should we also walk in opposite orientation? //TODO

        //set begin offset = end offset below to start collecting minimizers from adjacent vertices
        bases_processed = dfs_procedure (g, vertexId, g.readseq[i].length(), g.readseq[i].length(), depth_bases, visited_vertices, mmWalkRead, param);
      }

#ifdef VERBOSE
      std::cerr << "INFO, identifyRedundantReads(), processed readid " << i << ", collected " << mmWalkRead.size() << " minimizers from contained read\n";
#endif

      //collect minimizers from parent reads
      {
        visited_vertices.clear();
        for (uint32_t j = g.containment_offsets[i]; j < g.containment_offsets[i+1]; j++) {

          assert (g.containments[j].src == i);
          uint32_t parentReadId = g.containments[j].dst;
          uint32_t parentReadLen = g.readseq[parentReadId].length();

          if (g.redundant[parentReadId] == false)
          {
            uint32_t revbit = (g.containments[j].rev == true) ? 1U : 0U;
            uint32_t bases_processed = 0, beg, end, parentVertexId = parentReadId << 1 | revbit;

            if  (revbit == 0) {
              //suffix
              beg = g.containments[j].dst_end_offset;
              end = std::min (parentReadLen, beg + depth_bases);
            }
            else {
              //prefix
              end = g.containments[j].dst_start_offset;
              beg = end - std::min (end, depth_bases);
            }
            bases_processed = dfs_procedure (g, parentVertexId, beg, end, depth_bases, visited_vertices, mmWalkParentReads, param);
          }
        }
      }

#ifdef VERBOSE
      std::cerr << "INFO, identifyRedundantReads(), collected " << mmWalkParentReads.size() << " minimizers from " << available_parent_count << " parent reads\n";
#endif

      //compare minimizer lists and decide
      std::sort (mmWalkRead.begin(), mmWalkRead.end());
      std::sort (mmWalkParentReads.begin(), mmWalkParentReads.end());

      auto it1 = mmWalkRead.begin(), it2 = mmWalkParentReads.begin();

      uint32_t countCommon = 0;
      while (it1 != mmWalkRead.end() && it2 != mmWalkParentReads.end())
      {
        if (*it1 == *it2) {
          countCommon++;
          it1++; it2++;
        }
        else if (*it1 < *it2)
          it1++;
        else if (*it2 < *it1)
          it2++;
      }

      assert (countCommon <= mmWalkRead.size());
      assert (countCommon <= mmWalkParentReads.size());


      if (mmWalkRead.size() > 0 && 1.0 * countCommon / mmWalkRead.size() >= param.cutoff) {
        g.redundant[i] = true;

        log << g.umap_inverse[i] << "\tidentifyRedundantReads()\tPARENTS=";
        for (uint32_t j = g.containment_offsets[i]; j < g.containment_offsets[i+1]; j++) {
          uint32_t parentReadId = g.containments[j].dst;
          if (g.redundant[parentReadId] == false)
            log << g.umap_inverse[parentReadId] << ", ";
        }
        log << "\t" << countCommon << "/" << mmWalkRead.size() << ":" << mmWalkParentReads.size() << "\n";
      }
    }
  }

  std::cerr << "INFO, identifyRedundantReads() finished\n";
}

/**
 * we can avoid rescuing contained reads if they have
 * 0 out-degree or 0 in-degree
 * in-degree is approximated from out-degree of reverse read
 */
void pruneEndContainedReads(graphcontainer &g, std::ofstream& log)
{
  for (uint32_t i = 0; i < g.readCount; i++)
  {
    if (g.contained[i] == true && g.redundant[i] == false)
    {
      if (g.getDegree (i << 1 | 1) == 0 || g.getDegree (i << 1 | 0) == 0) {
        g.redundant[i] = true;
        log << g.umap_inverse[i] << "\tpruneEndContainedReads()\n";
      }
    }
  }
  std::cerr << "INFO, pruneEndContainedReads() finished\n";
}

/**
 * @param[in]   removeContainedReads      simply mark all contained reads as redundant
 */
void ovlgraph_simplify (bool removeContainedReads, graphcontainer &g, const algoParams &param, const std::string &logfilename)
{
  std::ofstream logFile (logfilename);

  if (removeContainedReads) g.removeContainedReads (logFile);
  g.index(); //indexing is needed prior to transitive reduction
  if (param.fuzz != INT32_MAX)
    g.transitiveReduction (param.fuzz);
  g.index(); //re-index

  //check redundancy of contained reads
  if (!removeContainedReads)
  {
    g.removeContainedReadsAboveDegree (param.maxContainmentDegree, logFile);
    g.index();

    while (true)
    {
      std::size_t countRedundantReads = std::count(g.redundant.begin(), g.redundant.end(), true);

      g.printGraphStats();
      identifyRedundantReads (g, param, logFile);
      g.index(); //re-index

      g.printGraphStats();
      pruneEndContainedReads(g, logFile);
      g.index(); //re-index

      if (std::count(g.redundant.begin(), g.redundant.end(), true) == countRedundantReads)
        break;
    }
  }

  std::cerr << "INFO, ovlgraph_simplify() finished, printing final stats\n";
  g.printGraphStats();
}

#endif
