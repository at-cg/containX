#ifndef ASM_ALGO_H
#define ASM_ALGO_H

#include <set>
#include <numeric>
#include "common.hpp"
#include "graph.hpp"
#include "param.hpp"
#include "other_algo.hpp"
#include <omp.h>

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
  uint32_t hashThreshold = UINT32_MAX * param.d;

  //adjust boundary by k-1 characters for parsing k-mers at junctions
  if (rev == false) {
    if (beg >= param.k - 1)
      beg -= param.k - 1;
  }
  else {
    if (str.length() - end >= param.k - 1)
      end += param.k - 1;
  }

  assert (beg >= 0 && beg <= str.length());
  assert (end >= beg && end <= str.length());

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
        if (hash32(kmer, mask) <= hashThreshold)
          container.emplace_back(kmer);
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
uint32_t collectMinimizersusingDFS (const graphcontainer &g, uint32_t src_vertex, uint32_t beg, uint32_t end, uint32_t remaining_depth_bases, std::set<uint32_t> &visited_vertices, std::vector<uint32_t> &minimizers, const algoParams &param)
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
  std::cerr << "INFO, collectMinimizersusingDFS(), computing minimizers from read " << src_readid << ", vertex = " << src_vertex << ", offsets = [" << beg << "," << end << "), orientation = '" << "+-"[rev] << "'\n";
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

      if (g.deletedReads[adjReadId] == false)
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
        bases_processed += collectMinimizersusingDFS (g, adjVertexId, beg, end, remaining_depth_bases, visited_vertices, minimizers, param);
      }
    }
  }
  return bases_processed;
}

void identifyRedundantReadsKmerHeuristic(const graphcontainer &g, std::vector<bool> &redundant, const algoParams &param, std::ofstream& log)
{
#pragma omp parallel
  {
    std::vector<uint32_t> mmWalkRead; //walk along same orientation as the read, and collect minimizers
    std::vector<uint32_t> mmWalkParentReads; //collect minimizers by walking from parent reads containing the current read
    std::vector<uint32_t> mmCommon; //common minimizers within the above two
    std::set<uint32_t> visited_vertices;

    //number ids of reads were assigned in the order they were found in paf file
    //static scheduling with chunk size one is preferred as
    //reads with lower ids may have more overlaps
#pragma omp for schedule(static, 1)
    for (uint32_t i = 0; i < g.readCount; i++)
    {
      uint32_t available_parent_count = 0;
      mmWalkRead.clear();
      mmWalkParentReads.clear();
      mmCommon.clear();

      if (g.contained[i] == true && g.deletedReads[i] == false)
      {
        assert (g.readseq[i].length() > 0);

        //user-specified depth during DFS in terms of count of bases
        uint32_t depth_bases = param.depthReadLen * g.readseq[i].length();

        //collect minimizers by starting DFS from contained read
        {
          visited_vertices.clear();
          uint32_t bases_processed;
          uint32_t vertexId = i << 1 | 0; //forward orientation
          //should we also walk in opposite orientation? //TODO

          //set begin offset = end offset below to start collecting minimizers from adjacent vertices
          bases_processed = collectMinimizersusingDFS (g, vertexId, g.readseq[i].length(), g.readseq[i].length(), depth_bases, visited_vertices, mmWalkRead, param);

          vertexId = i << 1 | 1; //reverse orientation
          bases_processed += collectMinimizersusingDFS (g, vertexId, g.readseq[i].length(), g.readseq[i].length(), depth_bases, visited_vertices, mmWalkRead, param);
        }

#ifdef VERBOSE
        std::cerr << "INFO, identifyRedundantReadsKmerHeuristic(), processed readid " << i << ", collected " << mmWalkRead.size() << " minimizers from contained read\n";
#endif

        if (mmWalkRead.size() == 0) continue; //no point going further

        //collect minimizers from parent reads
        {
          visited_vertices.clear();
          for (uint32_t j = g.containment_offsets[i]; j < g.containment_offsets[i+1]; j++) {

            assert (g.containments[j].src == i);
            uint32_t parentReadId = g.containments[j].dst;
            uint32_t parentReadLen = g.readseq[parentReadId].length();

            if (g.deletedReads[parentReadId] == false)
            {
              //walk w.r.t. forward orientation of read string
              uint32_t revbit = (g.containments[j].rev == true) ? 1U : 0U;
              uint32_t bases_processed = 0, beg, end;
              uint32_t parentVertexId = parentReadId << 1 | revbit;

              if  (revbit == 0) {
                //suffix of read string
                beg = g.containments[j].dst_end_offset;
                end = std::min (parentReadLen, beg + depth_bases);
              } else {
                //prefix of read string
                end = g.containments[j].dst_start_offset;
                beg = end - std::min (end, depth_bases);
              }
              bases_processed = collectMinimizersusingDFS (g, parentVertexId, beg, end, depth_bases, visited_vertices, mmWalkParentReads, param);

              //walk w.r.t. reverse orientation of read string
              revbit = (g.containments[j].rev == true) ? 0U : 1U;
              parentVertexId = parentReadId << 1 | revbit;

              if  (revbit == 0) {
                //suffix of read string
                beg = g.containments[j].dst_end_offset;
                end = std::min (parentReadLen, beg + depth_bases);
              } else {
                //prefix of read string
                end = g.containments[j].dst_start_offset;
                beg = end - std::min (end, depth_bases);
              }
              bases_processed += collectMinimizersusingDFS (g, parentVertexId, beg, end, depth_bases, visited_vertices, mmWalkParentReads, param);
            }
          }
        }

#ifdef VERBOSE
        std::cerr << "INFO, identifyRedundantReadsKmerHeuristic(), collected " << mmWalkParentReads.size() << " minimizers from " << available_parent_count << " parent reads\n";
#endif

        //keep unique minimizers only before comparing
        //rationale: there may be duplicate minimizers collected from redundant neighboring reads
        std::sort (mmWalkParentReads.begin(), mmWalkParentReads.end());
        auto last = std::unique (mmWalkParentReads.begin(), mmWalkParentReads.end());
        mmWalkParentReads.erase (last, mmWalkParentReads.end());

        std::sort (mmWalkRead.begin(), mmWalkRead.end());
        last = std::unique (mmWalkRead.begin(), mmWalkRead.end());
        mmWalkRead.erase (last, mmWalkRead.end());

        //common minimizers
        std::set_intersection(mmWalkParentReads.begin(), mmWalkParentReads.end(),
            mmWalkRead.begin(), mmWalkRead.end(),
            std::back_inserter(mmCommon));

        if (mmWalkRead.size() > 0 && 1.0 * mmCommon.size() / mmWalkRead.size() >= param.cutoff) {
          redundant[i] = true;
          if (!param.logFileName.empty()) {
#pragma omp critical
            log << g.umap_inverse.at(i) << "\tidentifyRedundantReadsKmerHeuristic()\tREDUNDANT=T\tPARENTS=";
          }
        }
        else {
          if (!param.logFileName.empty()) {
#pragma omp critical
            log << g.umap_inverse.at(i) << "\tidentifyRedundantReadsKmerHeuristic()\tREDUNDANT=F\tPARENTS=";
          }
        }

        for (uint32_t j = g.containment_offsets[i]; j < g.containment_offsets[i+1]; j++)
        {
          uint32_t parentReadId = g.containments[j].dst;
          if (g.deletedReads[parentReadId] == false)
            if (!param.logFileName.empty()) {
#pragma omp critical
              log << g.umap_inverse.at(parentReadId) << ", ";
            }
        }
        if (!param.logFileName.empty()) {
#pragma omp critical
          log << "\t" << mmCommon.size() << "/" << mmWalkRead.size() << ":" << mmWalkParentReads.size() << "\n";
        }
      }
    }
  }

  std::cerr << "INFO, identifyRedundantReadsKmerHeuristic() finished\n";
}

/**
 * @brief                       check whether one can spell string associated with vid via path from v_src
 *                              to final_dst vertex inside subgraph with allcombined vertices
 * @param[in] g
 * @param[in] offset_vid        zero-based offset of first character in read(vid) yet to be matched
 * @param[in] v_src             starting vertex
 * @param[in] offset_src        zero-based offset of first character in read(v_src) after overlap
 * @param[in] final_dst         destination vertex
 * @param[in] final_offset_dst  zero-based offset of first character in read(final_dst) after overlap
 * @param[in] vid               vertex id of the main read undergoing redundancy check
 * @param[in] subgraph_vertices restrict graph traversal within these vertices
 */
bool checkSequenceViaDFS (const graphcontainer &g, uint32_t offset_vid, uint32_t v_src, uint32_t offset_src, uint32_t final_dst, uint32_t final_offset_dst, uint32_t vid, const std::set<uint32_t> &subgraph_vertices, std::ofstream& log)
{
  assert (subgraph_vertices.find (v_src) != subgraph_vertices.end());
  assert (subgraph_vertices.find (final_dst) != subgraph_vertices.end());
  uint32_t vid_read_len = g.readseq[vid>>1].length();
  uint32_t src_read_len = g.readseq[v_src>>1].length();

  if (offset_vid == vid_read_len) //completed scanning all characters of main read
  {
    if (v_src == final_dst && offset_src == final_offset_dst) return true;
    else return false;
  }
  else
  {
    assert (offset_vid < vid_read_len);   //some characters in main read still pending
    assert (offset_src == src_read_len);  //v_src read is completely checked
    for (uint32_t j = g.offsets[v_src]; j < g.offsets[v_src+1]; j++) //neighbors of v_src
    {
      if (subgraph_vertices.find(g.edges[j].dst) == subgraph_vertices.end()) continue; //not in subgraph
      uint32_t countChars = std::min(g.edges[j].len, vid_read_len - offset_vid);
      assert(countChars > 0);

      std::string adjacentRead = g.readseq[g.edges[j].dst>>1];

      //reverse complement the adjacent read string if needed
      if (g.edges[j].dst & 1)
      {
        std::reverse(adjacentRead.begin(), adjacentRead.end());
        std::transform (adjacentRead.begin(), adjacentRead.end(), adjacentRead.begin(), complement);
      }
      std::string readSubStr = adjacentRead.substr (adjacentRead.length() - g.edges[j].len, countChars);

      std::string refSubStr = g.readseq[vid>>1].substr (offset_vid, countChars);
      if (readSubStr == refSubStr)
        if (checkSequenceViaDFS (g, countChars + offset_vid, g.edges[j].dst, adjacentRead.length() - g.edges[j].len + countChars, final_dst, final_offset_dst, vid, subgraph_vertices, log))
          return true;
    }
    return false;
  }
}

void identifyRedundantReadsDiscardLongerReads(const graphcontainer &g, std::vector<bool> &redundant, const algoParams &param, std::ofstream& log)
{
  //process reads in descending order of their length
  std::vector<uint32_t> idx (g.readCount);
  std::iota (idx.begin(), idx.end(), 0);

  //sort by descending order of read length
  std::sort (idx.begin(), idx.end(), [&](uint32_t i1, uint32_t i2) {return g.readseq[i1].length() > g.readseq[i2].length();});
  assert(g.checkSymmetry());
  uint32_t n_reduced = 0;

  //vertices either adjacent to vid or "contained" inside vid
  std::set<uint32_t> subgraph_vertices;
  std::set<uint32_t> incoming_vertices;
  std::set<uint32_t> outgoing_vertices;
  std::set<uint32_t> contained_vertices;

  for (uint32_t i = 0; i < g.readCount; i++)
  {
    uint32_t rid = idx[i];
    uint32_t vid = rid << 1 | 0; //forward orientation

    subgraph_vertices.clear();
    incoming_vertices.clear();
    outgoing_vertices.clear();
    contained_vertices.clear();

    if (redundant[rid]) continue;

    for (uint32_t j = g.offsets[vid]; j < g.offsets[vid+1]; j++)
      if (redundant[g.edges[j].dst >> 1] == false) {
        subgraph_vertices.insert (g.edges[j].dst);
        outgoing_vertices.insert (g.edges[j].dst);
      }

    for (uint32_t j = g.offsets[vid^1U]; j < g.offsets[(vid^1U) + 1]; j++)
      if (redundant[g.edges[j].dst >> 1] == false) {
        subgraph_vertices.insert (g.edges[j].dst^1U);
        incoming_vertices.insert (g.edges[j].dst^1U);
      }

    for (uint32_t j = g.containment_offsets_children[rid]; j < g.containment_offsets_children[rid+1]; j++) {
      assert (g.containments_copy[j].dst == rid); //"dst" is parent read containing "src"
      if (redundant[g.containments_copy[j].src] == false) {
        uint32_t revbit = (g.containments_copy[j].rev == true) ? 1U : 0U;
        subgraph_vertices.insert (g.containments_copy[j].src << 1 | revbit);
        contained_vertices.insert (g.containments_copy[j].src << 1 | revbit);
      }
    }

    if (subgraph_vertices.find(vid) != subgraph_vertices.end()) continue;
    if (subgraph_vertices.find(vid^1U) != subgraph_vertices.end()) continue;
    if (contained_vertices.size() == 0) continue; //no benefit in pruning reads that don't contain others

    bool checkSucceeded = false;
    bool checkedOnce = false;

    for (uint32_t j = g.offsets[vid^1U]; j < g.offsets[(vid^1U) + 1]; j++)
      for (uint32_t k = g.offsets[vid]; k < g.offsets[vid+1]; k++)
        if ((checkedOnce == false || checkSucceeded == true) && redundant[g.edges[j].dst >> 1] == false && redundant[g.edges[k].dst >> 1] == false) {
          checkSucceeded = checkSequenceViaDFS (g, g.edges[j].ov_src, g.edges[j].dst^1U, g.readseq[g.edges[j].dst>>1].length(), g.edges[k].dst, g.edges[k].ov_src, vid, subgraph_vertices, log);
          checkedOnce = true;
        }

    if (checkSucceeded == true) {
      n_reduced++;
      redundant[rid] = true;
      log << g.umap_inverse.at(rid) << "\tidentifyRedundantReadsDiscardLongerReads()\tREDUNDANT=T\t" << rid << "\tCHILDREN=" << contained_vertices.size() << "\tIN-DEGREE=" << incoming_vertices.size() << "\tOUT-DEGREE=" << outgoing_vertices.size() << std::endl;
    }
    else
      log << g.umap_inverse.at(rid) << "\tidentifyRedundantReadsDiscardLongerReads()\tREDUNDANT=F\t" << rid << "\tCHILDREN=" << contained_vertices.size() << "\tIN-DEGREE=" << incoming_vertices.size() << "\tOUT-DEGREE=" << outgoing_vertices.size() << std::endl;
  }

  std::cerr << "INFO, identifyRedundantReadsDiscardLongerReads(), " << n_reduced << " reads marked for deletion\n";
  std::cerr << "INFO, identifyRedundantReadsDiscardLongerReads() finished\n";
}


void ovlgraph_simplify (graphcontainer &g, const algoParams &param)
{
  std::ofstream logFile (param.logFileName);
  graphCleanup (g, param, logFile);

  uint32_t iter = 0;
  while (iter < param.iter)
  {
    if (param.runHui2016)
      identifyRedundantReadsHuiEtAl (g, param, logFile); //implemented for benchmarking
    else
    {
      std::vector<bool> redundant = g.deletedReads;
      identifyRedundantReadsDiscardLongerReads (g, redundant, param, logFile); //our algorithm

      g.deletedReads = redundant;
      g.index(); //re-index
      g.printGraphStats();
      redundant = g.deletedReads;

      identifyRedundantReadsKmerHeuristic (g, redundant, param, logFile); //our algorithm
      g.deletedReads = redundant;
    }
    g.index(); //re-index
    g.printGraphStats();

    tipCleaning (g, param, logFile);
    g.index(); //re-index
    g.printGraphStats();
    iter++;
  }

  std::cerr << "INFO, ovlgraph_simplify() finished, printing final stats\n";
  g.printGraphStats();
}

#endif
