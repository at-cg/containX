#ifndef ASM_GRAPH_H
#define ASM_GRAPH_H

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fstream>
#include "paf.hpp"
#include "common.hpp"
#include "param.hpp"

/*
                    |< ov_src >|
  src: ------------------------>
                    ||overlap|||
               dst: -------------------------->
                    |< ov_dst >|<---- len --->|

  * if a read x has a suffix-prefix overlap with another read y, arcs x->y and ~y -> ~x are saved
  * vertex ids will be formed as (read id << 1 | orientation); therefore two vertices per read
  * we will maintain edges sorted by key <src, len>
*/
class graphArc
{
  public:
    uint32_t src;
    uint32_t len;
    uint32_t dst;
    uint32_t ov_src;  //strand, length of suffix of src involved in overlap
    uint32_t ov_dst;  //strand, length of prefix of dst involved in overlap
    bool del;         //mark for deletion if needed

    /**
     * The above storage format makes it easy to export in GFA format
     * https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#l-link-line
     */

};

/**
 * Tuples to record which reads are contained in which reads
 * Unlike graph arcs, we will store read ids instead of vertex ids
 */
class containmentTuple
{
  public:
    uint32_t src;              //read which is contained
    uint32_t dst;              //parent read
    uint32_t dst_start_offset; //(0-based; BED-like; closed) leftmost offset in dst's string
    uint32_t dst_end_offset;   //(0-based; BED-like; open)
    uint32_t rev;              //dst orientation

    /**
     * The above storage format makes it easy to export in GFA format
     * https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#c-containment-line
     */
};

class graphcontainer
{
  public:

    uint32_t readCount;
    uint32_t vertexCount;       //twice of readCount

    //hash: read id -> number
    std::unordered_map <std::string, uint32_t> umap;  // size = count of reads
    //inverse hash: number -> read id
    std::unordered_map <uint32_t, std::string> umap_inverse;  // size = count of reads
    //read sequences
    std::vector<std::string> readseq;  //size = count of reads
    std::vector<bool> contained; //size = count of reads
    std::vector<bool> deletedReads;   //size = count of reads
    std::vector<graphArc> edges; //size = 2 x suffix-prefix overlaps
    std::vector<uint32_t> offsets; //for CSR-style indexing
    std::vector<containmentTuple> containments; //size = count of contained overlaps
    std::vector<uint32_t> containment_offsets; //for CSR-style indexing (know id of parent reads)
    std::vector<containmentTuple> containments_copy; //size = count of contained overlaps (same content but ordered differently than containments vector)
    std::vector<uint32_t> containment_offsets_children; //for CSR-style indexing (know id of children reads)

    //save fastq read identifier into hash table, and give it an integer id
    uint32_t addStringToMap(const std::string &str)
    {
      if (umap.find(str) != umap.end())
      {
        return umap[str];
      }
      else
      {
        uint32_t key = (uint32_t) umap.size();
        umap[str] = key;
        return key;
      }
    }

    //initialize basic vectors
    void initVectors()
    {
      readCount = umap.size();
      vertexCount = 2 * readCount;

      if (readCount == 0) return;
      assert (contained.size() == 0);
      assert (deletedReads.size() == 0);
      assert (readseq.size() == 0);

      contained.resize(readCount, false);
      deletedReads.resize(readCount, false);
      readseq.resize(readCount, "");
      inverse_map (umap, umap_inverse); //build inverse

      for (auto &e : containments)
        contained[e.src] = true;
    }

    //parse + save all reads, and mark the contained ones
    void initReadStrings(const char *readfilename)
    {
      //parse reads
      int l;
      gzFile fp;
      kseq_t *seq;
      fp = gzopen(readfilename, "r"); //open the file handler
      seq = kseq_init(fp);
      while ((l = kseq_read(seq)) >= 0) {
        std::string id = seq->name.s;
        std::string readstr = seq->seq.s;
        if (umap.find(id) != umap.end() && deletedReads[umap[id]] == false)
        {
          readseq[umap[id]] = readstr;
        }
      }
      kseq_destroy(seq);
      gzclose(fp);
      std::cerr << "INFO, initReadStrings(), parsed reads from " << readfilename << "\n";
    }

    void printGraphStats ()
    {
      std::cerr << "INFO, printGraphStats(), graph has " << edges.size() << " edges\n";
      std::cerr << "INFO, printGraphStats(), multiedges : " << std::boolalpha << checkMultiEdges() << ", symmetric : " << checkSymmetry() << "\n";
      std::cerr << "INFO, printGraphStats(), graph has " << vertexCount << " vertices from " << readCount << " reads in total\n";
      std::cerr << "INFO, printGraphStats(), " << std::count(contained.begin(), contained.end(), true) << " reads are marked as contained in graph\n";

      uint32_t junctionReads, i;

      for (i = 0, junctionReads = 0; i < readCount; i++)
        if (deletedReads[i] == false)
          if (getDegree(i << 1 | 0) > 1 || getDegree(i << 1 | 1) > 1)
            junctionReads++;
      std::cerr << "INFO, printGraphStats(), " << junctionReads << " reads contribute to junction nodes\n";

      for (i = 0, junctionReads = 0; i < readCount; i++)
        if (contained[i] == true && deletedReads[i] == false)
          if (getDegree(i << 1 | 0) > 1 || getDegree(i << 1 | 1) > 1)
            junctionReads++;
      std::cerr << "INFO, printGraphStats(), " << junctionReads << " contained reads contribute to junction nodes\n";


    }

    //index edges and containments
    void index()
    {
      indexEdges();
      indexContainments();
    }

    //sort edge vectors and index using CSR format
    void indexEdges()
    {
      //get rid of deleted edges as well as edges connecting to deleted vertices
      {
        auto it = std::remove_if (edges.begin(), edges.end(), [this](const graphArc &a) {return deletedReads[a.src >> 1] || deletedReads[a.dst >> 1] || a.del; });
        edges.erase (it, edges.end());
      }

      //order such that edges with same "src" occur next to each other
      //within each bucket, prefer edges of bigger length first to optimize "identifyRedundantReadsDiscardLongerReads()"
      std::sort (edges.begin(), edges.end(), [](const graphArc &a, const graphArc &b) {
          return std::tie (a.src, b.len, a.dst) < std::tie (b.src, a.len, b.dst);});

      //make sure there are no duplicate entries
      auto last = std::unique (edges.begin(), edges.end(), [](const graphArc &a, const graphArc &b) {
          return std::tie (a.src, a.len, a.dst) == std::tie (b.src, b.len, b.dst);});
      edges.erase (last, edges.end());

      /**
       * Build offsets array such that out-edges
       * of vertex i are accessible from edges[offset[i]]
       * inclusive to edges[offset[i] + 1] exclusive
       */
      if (vertexCount == 0) return;
      offsets.resize(vertexCount + 1, 0);

      auto it_b = edges.begin();

      for(uint32_t i = 0; i < vertexCount; i++)
      {
        //Range for adjacency list of vertex i
        auto it_e = std::find_if(it_b, edges.end(), [i](const graphArc &e) { return e.src > i; });
        offsets[i+1] = std::distance(edges.begin(), it_e);
        it_b = it_e;
      }

      assert (it_b == edges.end());
    }

    //index containment relationships, for fast retrieval of all reads containing a particular read
    void indexContainments()
    {
      //get rid of containment relationships involving deleted reads
      {
        auto it = std::remove_if (containments.begin(), containments.end(), [this](const containmentTuple &a) {return deletedReads[a.src] || deletedReads[a.dst]; });
        containments.erase (it, containments.end());

        it = std::remove_if (containments_copy.begin(), containments_copy.end(), [this](const containmentTuple &a){return deletedReads[a.src] || deletedReads[a.dst]; });
        containments_copy.erase (it, containments_copy.end());
      }

      std::sort (containments.begin(), containments.end(), [](const containmentTuple &a, const containmentTuple &b) { return a.src < b.src;});

      std::sort (containments_copy.begin(), containments_copy.end(), [](const containmentTuple &a, const containmentTuple &b) { return a.dst < b.dst;});

      if (readCount == 0) return;
      containment_offsets.resize(readCount + 1, 0);
      containment_offsets_children.resize(readCount + 1, 0);

      auto it_b = containments.begin();

      for(uint32_t i = 0; i < readCount; i++)
      {
        //Range for adjacency list of vertex i (i.e., all parent reads of read i)
        auto it_e = std::find_if(it_b, containments.end(), [i](const containmentTuple &e) { return e.src > i; });
        containment_offsets[i+1] = std::distance(containments.begin(), it_e);
        it_b = it_e;
      }

      assert (it_b == containments.end());

      it_b = containments_copy.begin();

      for(uint32_t i = 0; i < readCount; i++)
      {
        //Range for adjacency list of vertex i (i.e., all child reads of read i)
        auto it_e = std::find_if(it_b, containments_copy.end(), [i](const containmentTuple &e) { return e.dst > i; });
        containment_offsets_children[i+1] = std::distance(containments_copy.begin(), it_e);
        it_b = it_e;
      }

      assert (it_b == containments_copy.end());
      assert (containments_copy.size() == containments.size());
      assert (containment_offsets.size() == containment_offsets_children.size());


      contained.assign (readCount, false);
      for (auto &e : containments)
        contained[e.src] = true;
    }

    //count of out-edges from a graph vertex
    uint32_t getDegree (uint32_t src) const
    {
      assert (offsets.size() == vertexCount + 1);
      assert (src < vertexCount);

      return offsets[src + 1] - offsets[src];
    }

    //return maximum degree value in graph
    uint32_t maxDegree () const
    {
      uint32_t maxdeg = 0;
      for (uint32_t i = 0; i < vertexCount; i++)
        maxdeg = std::max (maxdeg, getDegree(i));

      return maxdeg;
    }

    //count of reads containing the read 'src'
    uint32_t getContaintmentDegree (uint32_t src) const
    {
      assert (containment_offsets.size() == readCount + 1);
      assert (src < readCount);

      return containment_offsets[src + 1] - containment_offsets[src];
    }

    //also see https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
    void outputGFA (const std::string &filename, bool printReadStrings) const
    {
      std::ofstream outstrm (filename);

      //print reads - use prefix 'read' before their id
      for (uint32_t i = 0; i < readCount; i++)
        if (deletedReads[i] == false)
          if (printReadStrings)
            outstrm  << "S\tread" << i << "\t" << readseq[i] << "\n";
          else
            outstrm  << "S\tread" << i << "\t*\n";

      //print arcs
      for (uint32_t i = 0; i < edges.size(); i++) {
        if (edges[i].del == false)
        {
          assert (deletedReads[edges[i].src >> 1] == false);
          assert (deletedReads[edges[i].dst >> 1] == false);
          outstrm << "L\tread" << (edges[i].src >> 1) << "\t" << "+-"[edges[i].src & 1] <<  "\tread" << (edges[i].dst >> 1) << "\t" << "+-"[edges[i].dst & 1] << "\t" << edges[i].ov_src <<  "M\n";
        }
      }

      //print containment lines
      for (uint32_t i = 0; i < containments.size(); i++)
        if (deletedReads[containments[i].src] == false && deletedReads[containments[i].dst] == false)
          outstrm << "C\tread" << containments[i].dst << "\t" << "+-"[containments[i].rev] << "\tread" << containments[i].src << "\t+\t" << containments[i].dst_start_offset << "\t" << containments[i].dst_end_offset - containments[i].dst_start_offset<< "M\n";

      //print summary of reads, may help in debugging
      for (auto &e : umap)
      {
        uint32_t i = e.second; //read id
        //print original read id, length, count of reads containing it, out-degree (fwd), out-degree (rev)
        if (deletedReads[i] == false) {
            outstrm  << "x\tread" << i << "\t" << e.first << "\t" << readseq[i].length() << "\t" << getContaintmentDegree(i) << "\t" << getDegree (i << 1 | 0) << "\t" << getDegree (i << 1 | 1) << "\n";
        }
      }
    }

    //output in graph format similar to what is used by PACE competition
    //https://github.com/PACE-challenge/Treewidth
    //exported graph will be undirected, and will not contain self loops
    //vertex ids in the exported graph are 1-based
    void exportPACE (const std::string &filename) const
    {
      if (filename.empty()) return;
      std::ofstream outstrm (filename);

      std::vector<std::pair<uint32_t, uint32_t>> edgesToPrint;

      for (uint32_t i = 0; i < edges.size(); i++) {
        assert (edges[i].del == false);
        if (edges[i].src < edges[i].dst)
          edgesToPrint.emplace_back(edges[i].src + 1, edges[i].dst + 1);
        else if (edges[i].src > edges[i].dst)
          edgesToPrint.emplace_back(edges[i].dst + 1, edges[i].src + 1);
      }

      std::sort (edgesToPrint.begin(), edgesToPrint.end());
      edgesToPrint.erase(std::unique(edgesToPrint.begin(), edgesToPrint.end()), edgesToPrint.end());

      outstrm  << "p tw " << vertexCount << " " << edgesToPrint.size() << "\n";
      for (auto &e: edgesToPrint)
        outstrm  << e.first << " " << e.second << "\n";
    }

    /**
     * print read ids which are contained but not redundant
     */
    void outputNonRedudantContainedReads (const std::string &filename)
    {
      if (filename.empty()) return;
      std::ofstream outstrm (filename);

      for (auto &e : umap)
      {
        uint32_t i = e.second; //read id
        //print original read id
        if (contained[i] == true && deletedReads[i] == false) {
            outstrm  << e.first << "\n";
        }
      }
    }

    /**
     * print read ids which are contained but not redundant
     */
    void outputNonRedudantReads (const std::string &filename)
    {
      if (filename.empty()) return;
      std::ofstream outstrm (filename);

      for (auto &e : umap)
      {
        uint32_t i = e.second; //read id
        //print original read id
        if (deletedReads[i] == false) {
            outstrm  << e.first << "\n";
        }
      }
    }

    //remove u->v if v'->u' is already removed
    void ensureSymmetry ()
    {
      uint32_t n_reduced = 0;
      for (uint32_t i = 0, j = 0; i < edges.size(); i++)
      {
        if (edges[i].del == true) continue;

        //u->v
        uint32_t u = edges[i].src, v = edges[i].dst;
        uint32_t u_rev = u^1U, v_rev = v^1U;
        //check for v_rev -> u_rev
        for (j = offsets[v_rev]; j < offsets[v_rev +1]; j++)
        {
          assert (edges[j].src == v_rev);
          if (edges[j].dst == u_rev && edges[j].ov_dst == edges[i].ov_src && edges[j].del == false) break; //found
        }
        if (j == offsets[v_rev +1]) {
          edges[i].del = true;
          n_reduced++;
        }
      }

      std::cerr << "INFO, ensureSymmetry() finished, " << n_reduced << " edges marked for deletion\n";
    }

    //assumes indexing is done and edges are sorted
    bool checkMultiEdges ()
    {
      std::vector<bool> markAdj (vertexCount, false);
      for (uint32_t i = 0; i < vertexCount; i++)
      {
        if (getDegree(i) <= 1) continue;
        for (uint32_t j = offsets[i]; j < offsets[i+1]; j++)
        {
          if (markAdj [edges[j].dst] == true) return true; //multi
          markAdj [edges[j].dst] = true;
        }
        for (uint32_t j = offsets[i]; j < offsets[i+1]; j++)
        {
          markAdj [edges[j].dst] = false; //reset for next iteration
        }
      }

      return false;
    }

    //assumes indexing is done and edges are sorted
    bool checkSymmetry () const
    {
      for (uint32_t i = 0, j = 0; i < edges.size(); i++)
      {
        //u->v
        uint32_t u = edges[i].src, v = edges[i].dst;
        uint32_t u_rev = u^1, v_rev = v^1;
        //check for v_rev -> u_rev
        for (j = offsets[v_rev]; j < offsets[v_rev +1]; j++)
        {
          if (edges[j].dst == u_rev && edges[j].ov_dst == edges[i].ov_src) break;
        }
        if (j == offsets[v_rev +1]) return false;
      }
      return true;
    }

};

void processHetReads (const std::unordered_set<std::string> &hetReadsToIgnoreIfContained, graphcontainer &g)
{
  uint32_t n_del = 0;
  for (const auto& readIdStr: hetReadsToIgnoreIfContained) {
    if (g.umap.find(readIdStr) != g.umap.end()) {
      uint32_t id = g.umap[readIdStr];
      if (g.contained[id] == true) {
        if (g.deletedReads[id] == false) {
          n_del++;
          g.deletedReads[id] = true;
        }
      }
    }
  }

  std::cerr << "INFO, processHetReads(), " << n_del << " reads marked for deletion\n";
}

void processHomReads (const std::unordered_set<std::string> &homReadsToIgnoreIfContained, graphcontainer &g)
{
  uint32_t n_del = 0;
  for (auto &e : g.containments)
  {
    if (g.deletedReads[e.src] == true) continue;
    if (g.deletedReads[e.dst] == true) continue;

    assert (g.umap_inverse.find(e.src) != g.umap_inverse.end());
    assert (g.umap_inverse.find(e.dst) != g.umap_inverse.end());
    auto containedReadIdStr = g.umap_inverse[e.src];
    auto parentReadIdStr = g.umap_inverse[e.dst];

    if (homReadsToIgnoreIfContained.find(containedReadIdStr) != homReadsToIgnoreIfContained.end() &&
      homReadsToIgnoreIfContained.find(parentReadIdStr) != homReadsToIgnoreIfContained.end()) {
      if (g.deletedReads[e.src] == false) {
        n_del++;
        g.deletedReads[e.src] = true;
      }
    }
  }

  std::cerr << "INFO, processHomReads(), " << n_del << " reads marked for deletion\n";
}

//note: we assume that overlapper skipped dual mappings, minimap2 overlapping module skips by default
//note: substring overlaps which are not suffix-prefix are ignored, may need to relax this later
void ovlgraph_gen(const char *readfilename, const char *paffilename, const algoParams &param, graphcontainer &g)
{
  //check if het read list from hifiasm is given by user
  std::unordered_set<std::string> hetReadsToIgnoreIfContained;
  if (!param.hetReads.empty())
  {
    std::string str;
    std::ifstream fs(param.hetReads);
    while(getline(fs,str))
      hetReadsToIgnoreIfContained.insert(str);

    std::cerr << "INFO, ovlgraph_gen(), parsed " << hetReadsToIgnoreIfContained.size() << " non-repetitive heterozygous reads from input file " << param.hetReads << "\n";
  }

  //read input paf file
  paf_rec_t r;
  paf_file_t *fp = paf_open(paffilename);
  std::cerr << "INFO, ovlgraph_gen(), reading paf records from " << paffilename << "\n";

  uint64_t totPaf = 0, validPaf = 0, suffPrefPaf = 0, containedPaf = 0;
  while (paf_read(fp, &r) >= 0) {
    totPaf++;
    if (r.ql >= param.min_read_len && r.tl >= param.min_read_len && r.ml * 100.0 / r.bl >= param.min_ovlp_identity)
    {
      if (std::min (r.te - r.ts, r.qe - r.qs) >= param.min_ovlp_len)
      {
        validPaf++;
        std::string qname = r.qn;
        std::string tname = r.tn;

        uint32_t q_readId = g.addStringToMap(qname);
        uint32_t t_readId = g.addStringToMap(tname);

        if (q_readId == t_readId) continue;

        if (r.qe - r.qs == r.ql) //qry is contained in target
        {
          containmentTuple c;
          c.src = q_readId;
          c.dst = t_readId;
          c.rev = r.rev;
          c.dst_start_offset = r.ts;
          c.dst_end_offset = r.te;
          g.containments.emplace_back(c);
          containedPaf++;
        }
        else if (r.te - r.ts == r.tl) //target is contained in qry
        {
          containmentTuple c;
          c.src = t_readId;
          c.dst = q_readId;
          c.rev = r.rev;
          c.dst_start_offset = r.qs;
          c.dst_end_offset = r.qe;
          g.containments.emplace_back(c);
          containedPaf++;
        }

        /**
         * There can be scenarios where read A is contained in read B, and
         * the opposite also holds true, that should be ok
         */

        /**
         * Conditions below check for suffix-prefix overlaps
         * See Insight 13.2 Reverse complements in assembly graphs in MBCT textbook for details
         */

        //a suffix of qry overlaps a prefix of target
        if (r.rev == 0 && r.qs > 0 && r.qe == r.ql && r.ts == 0 && r.te < r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = q_readId << 1 | 0;
            e1.dst = t_readId << 1 | 0;
            e1.ov_src = r.qe - r.qs;
            e1.ov_dst = r.te - r.ts;
            e1.len = r.tl - e1.ov_dst;
          }

          {
            e2.src = t_readId << 1 | 1;
            e2.dst = q_readId << 1 | 1;
            e2.ov_src = r.te - r.ts;
            e2.ov_dst = r.qe - r.qs;
            e2.len = r.ql - e2.ov_dst;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
        }
        //a prefix of qry overlaps a suffix of target
        else if (r.rev == 0 && r.qs == 0 && r.qe < r.ql && r.ts > 0 && r.te == r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = t_readId << 1 | 0;
            e1.dst = q_readId << 1 | 0;
            e1.ov_src = r.te - r.ts;
            e1.ov_dst = r.qe - r.qs;
            e1.len = r.ql - e1.ov_dst;
          }

          {
            e2.src = q_readId << 1 | 1;
            e2.dst = t_readId << 1 | 1;
            e2.ov_src = r.qe - r.qs;
            e2.ov_dst = r.te - r.ts;
            e2.len = r.tl - e2.ov_dst;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
        }
        //a suffix of ~qry overlaps a prefix of target
        else if (r.rev == 1 && r.qs == 0 && r.qe < r.ql && r.ts == 0 && r.te < r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = q_readId << 1 | 1;
            e1.dst = t_readId << 1 | 0;
            e1.ov_src = r.qe - r.qs;
            e1.ov_dst = r.te - r.ts;
            e1.len = r.tl - e1.ov_dst;
          }

          {
            e2.src = t_readId << 1 | 1;
            e2.dst = q_readId << 1 | 0;
            e2.ov_src = r.te - r.ts;
            e2.ov_dst = r.qe - r.qs;
            e2.len = r.ql - e2.ov_dst;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
        }
        //a prefix of ~qry overlaps a suffix of target
        else if (r.rev == 1 && r.qs > 0 && r.qe == r.ql && r.ts > 0 && r.te == r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = t_readId << 1 | 0;
            e1.dst = q_readId << 1 | 1;
            e1.ov_src = r.te - r.ts;
            e1.ov_dst = r.qe - r.qs;
            e1.len = r.ql - e1.ov_dst;
          }

          {
            e2.src = q_readId << 1 | 0;
            e2.dst = t_readId << 1 | 1;
            e2.ov_src = r.qe - r.qs;
            e2.ov_dst = r.te - r.ts;
            e2.len = r.tl - e2.ov_dst;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
        }
      }
    }

    if (totPaf % 10000000 == 0)
      std::cerr << "INFO, ovlgraph_gen(), parsed " << totPaf << " paf records\n";
  }

  paf_close(fp);

  std::cerr << "INFO, ovlgraph_gen(), parsed " << totPaf << " paf records\n";
  std::cerr << "INFO, ovlgraph_gen(), " << validPaf << " records satisfied user-specified cutoffs\n";
  std::cerr << "INFO, ovlgraph_gen(), " << containedPaf << " records belonged to contained overlaps\n";
  std::cerr << "INFO, ovlgraph_gen(), " << suffPrefPaf << " records belonged to proper suffix-prefix overlaps\n";

  std::for_each(g.edges.begin(), g.edges.end(), [](graphArc &e){e.del = false;});
  assert (g.edges.size() <= UINT32_MAX); //otherwise our implementation may not work
  g.containments_copy = g.containments;
  g.initVectors();
  processHetReads (hetReadsToIgnoreIfContained, g);
  g.initReadStrings(readfilename);
  g.index();
  g.printGraphStats();

  std::cerr << "INFO, ovlgraph_gen() finished\n";
}

/**
 * suppose containment degree of a read equals the count of
 * reads containing it; then the following function prints
 * distribution of containment degrees in the graph
 */
void printContainmentDegreeDistribution (graphcontainer &g, const std::string &filename)
{
  uint32_t maxDegree = 0;
  for (uint32_t i = 0; i < g.readCount; i++)
    if (g.deletedReads[i] == false)
      maxDegree = std::max (maxDegree, g.getContaintmentDegree(i));

  std::vector<uint32_t> distribution(maxDegree+1, 0);

  for (uint32_t i = 0; i < g.readCount; i++)
    if (g.deletedReads[i] == false)
      distribution[g.getContaintmentDegree(i)]++;

  //write to file
  std::ofstream outFile(filename);
  for (const auto &e : distribution) outFile << e << "\n";
}

/**
 * the following function prints
 * distribution of vertex out-degree in the graph
 * write to file Degree.txt (overwrite if already exists)
 */
void printDegreeDistribution (graphcontainer &g, const std::string &filename)
{
  uint32_t maxDegree = 0;
  for (uint32_t i = 0; i < g.vertexCount; i++)
    if (g.deletedReads[i>>1] == false)
      maxDegree = std::max (maxDegree, g.getDegree(i));

  std::vector<uint32_t> distribution(maxDegree+1, 0);

  for (uint32_t i = 0; i < g.vertexCount; i++)
    if (g.deletedReads[i>>1] == false)
      distribution[g.getDegree(i)]++;

  //write to file
  std::ofstream outFile(filename);
  for (const auto &e : distribution) outFile << e << "\n";
}

/**
 * the following function prints
 * distribution of vertex out-degree in the graph
 * write to file Degree.txt (overwrite if already exists)
 */
void printDegreeDistributionOnlyContainedVertices (graphcontainer &g, const std::string &filename)
{
  uint32_t maxDegree = 0;
  for (uint32_t i = 0; i < g.vertexCount; i++)
    if (g.contained[i>>1]==true && g.deletedReads[i>>1]==false)
      maxDegree = std::max (maxDegree, g.getDegree(i));

  std::vector<uint32_t> distribution(maxDegree+1, 0);

  for (uint32_t i = 0; i < g.vertexCount; i++)
    if (g.contained[i>>1]==true && g.deletedReads[i>>1]==false)
      distribution[g.getDegree(i)]++;

  //write to file
  std::ofstream outFile(filename);
  for (const auto &e : distribution) outFile << e << "\n";
}

/**
 * the following function prints
 * list of directed edges in the graph
 * write to file edges.DOT (overwrite if already exists)
 */
void printEdgesDOTFormat (graphcontainer &g, const std::string &filename)
{
  std::ofstream outFile(filename);
  outFile << "digraph overlaps {\n";
  for (auto &e : g.edges)
    if (g.deletedReads[e.src >>1] == false && g.deletedReads[e.dst >>1] == false && e.del == false)
      outFile << e.src << " -> " << e.dst << ";\n";
  outFile << "}\n";
}

//consider all contained reads as redundant and remove them from graph
void removeAllContainedReads(graphcontainer &g, const algoParams &param, std::ofstream &log)
{
  for (uint32_t i = 0; i < g.readCount; i++) {
    if (g.deletedReads[i] == true || g.contained[i] == false) continue;
    g.deletedReads[i] = true;
    if (!param.logFileName.empty()) log << g.umap_inverse[i] << "\tremoveAllContainedReads()\n";
  }

  std::vector<graphArc> edges_new;

  for (auto &e : g.edges)
  {
    //check if both end vertices of the edge are still available
    if (g.deletedReads[e.src >> 1] == false && g.deletedReads[e.dst >> 1] == false)
      edges_new.emplace_back(e);
  }

  g.edges = edges_new;
  std::cerr << "INFO, removeAllContainedReads() finished\n";
}

//algorithm motivated from Myers 2005
//assumes indexing is done
//revised to handle multi-graphs
uint32_t transitiveReduction(graphcontainer &g, int fuzz, std::ofstream& log)
{
  std::for_each(g.edges.begin(), g.edges.end(), [](graphArc &e){assert (e.del == false);});
  if (fuzz < 0) return 0U; //stop if reduction is disabled by user

  //save length of edge from vertex being considered
  uint32_t n_reduced = 0;

  for (auto &e : g.edges)
  {
    uint32_t v = e.src;
    if (g.getDegree (v) == 0) continue;

    //neighborhood of v
    for (uint32_t j = g.offsets[v]; j < g.offsets[v+1]; j++) {

      uint32_t w = g.edges[j].dst;

      //neighborhood of w
      for (uint32_t k = g.offsets[w]; k < g.offsets[w+1]; k++) {

        uint32_t sum = g.edges[j].len + g.edges[k].len;

        if (g.edges[k].dst == e.dst) {
          if (sum <= e.len + fuzz && sum + fuzz >= e.len ) { //avoid subtraction with unsigned types
            if (e.del == false) {
              e.del = true;
              n_reduced++;
            }
          }
        }
      }
    }
  }

  std::cerr << "INFO, transitiveReduction(), " << n_reduced << " edges marked for deletion\n";

  return n_reduced;
}

void tipCleaning (graphcontainer &g, const algoParams &param, std::ofstream& log)
{
  std::vector<uint32_t> tipVertexIds;

  for (uint32_t i = 0; i < g.vertexCount; i++)
  {
    tipVertexIds.clear();
    if (g.deletedReads[i >> 1] == true) continue;     //removed already
    if (g.getDegree (i ^ 1) != 0) continue;           //not a tip if there are incoming edges
    if (g.getDegree (i) > 1) continue;                //multiple out-edges
    if (g.getDegree (i) == 0)  {                      //singleton vertex without in and out neighbor
      if (!param.logFileName.empty()) log << g.umap_inverse[i >> 1] << "\ttipCleaning(singleton)\n";
      g.deletedReads[i >> 1] = true;
      continue;
    }

    tipVertexIds.push_back(i);
    uint32_t chainLen = 1; //in term of edge counts
    uint32_t currentVertex = i;
    bool validTip = false;

    for (uint32_t j = 0; j < param.maxTipLen; j++)
    {
      uint32_t next = g.edges [g.offsets[currentVertex]].dst;
      if (g.deletedReads[next >> 1] == false && g.getDegree (next ^ 1) == 1 && g.getDegree (next) == 1) {
        tipVertexIds.push_back (next);
        currentVertex = next;
      }
      else {
        validTip = true; //tip is short enough for pruning
      }
    }

    if (validTip) {
      for (uint32_t &i : tipVertexIds) {
        if (!param.logFileName.empty()) log << g.umap_inverse[i >> 1] << "\ttipCleaning(chain)\n";
        g.deletedReads[i >> 1] = true;
      }
    }
  }

  std::cerr << "INFO, tipCleaning() finished\n";
}

/*
 * Function does the following to clean graph:
 *   Apply minimum contained read length, maximum parent count thresholds
 *   Apply minimum overlap ratio threshold
 *   Transitive reduction
 *   Clean tips
 * Return the graph indexed
 */
void graphCleanup(graphcontainer &g, const algoParams &param, std::ofstream& log)
{
  if (param.removeAllContainedReads) removeAllContainedReads (g, param, log);

  g.index();
  g.printGraphStats();

  transitiveReduction (g, param.fuzz, log);
  g.index(); //re-index
  assert(g.checkSymmetry());
  g.printGraphStats();

  tipCleaning (g, param, log);
  g.index(); //re-index
  g.printGraphStats();
}

#endif
