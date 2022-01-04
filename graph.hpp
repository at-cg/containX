#ifndef ASM_GRAPH_H
#define ASM_GRAPH_H

#include <algorithm>
#include <tuple>
#include <unordered_map>
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
    uint32_t dst;
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
    std::vector<uint32_t> containment_offsets; //for CSR-style indexing

    //save user thresholds
    float min_ovlp_identity;
    int min_ovlp_len;

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

    //parse + save all reads, and mark the contained ones
    void initVectors(const char *readfilename)
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

      //parse reads
      {
        int l;
        gzFile fp;
        kseq_t *seq;
        fp = gzopen(readfilename, "r"); //open the file handler
        seq = kseq_init(fp);
        while ((l = kseq_read(seq)) >= 0) {
          std::string id = seq->name.s;
          std::string readstr = seq->seq.s;
          if (umap.find(id) != umap.end())
          {
            readseq[umap[id]] = readstr;
          }
        }
        kseq_destroy(seq);
        gzclose(fp);
      }
      std::cerr << "INFO, initVectors(), parsed reads from " << readfilename << "\n";
    }

    void printGraphStats ()
    {
      std::cerr << "INFO, printGraphStats(), graph has " << edges.size() << " edges\n";
      std::cerr << "INFO, printGraphStats(), multiedges : " << std::boolalpha << checkMultiEdges() << ", symmetric : " << checkSymmetry() << "\n";
      std::cerr << "INFO, printGraphStats(), graph has " << vertexCount << " vertices from " << readCount << " reads in total\n";
      std::cerr << "INFO, printGraphStats(), " << std::count(contained.begin(), contained.end(), true) << " reads are marked as contained in graph\n";
      uint32_t redundantcontainedReads = 0;
      for (uint32_t i = 0; i < readCount; i++)
        if (contained[i] == true && deletedReads[i] == true)
          redundantcontainedReads++;
      std::cerr << "INFO, printGraphStats(), " << redundantcontainedReads << " reads are marked as deleted and contained in graph\n";

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

      std::sort (edges.begin(), edges.end(), [](const graphArc &a, const graphArc &b) {
          return std::tie (a.src, a.len, a.dst) < std::tie (b.src, b.len, b.dst);});

      //make sure there are no duplicate entries
      auto last = std::unique (edges.begin(), edges.end(), [](const graphArc &a, const graphArc &b) {
          return std::tie (a.src, a.len, a.dst) == std::tie (b.src, b.len, b.dst);});
      assert (last == edges.end());

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
      }

      std::sort (containments.begin(), containments.end(), [](const containmentTuple &a, const containmentTuple &b) {
          return a.src < b.src;});

      if (readCount == 0) return;
      containment_offsets.resize(readCount + 1, 0);

      auto it_b = containments.begin();

      for(uint32_t i = 0; i < readCount; i++)
      {
        //Range for adjacency list of vertex i
        auto it_e = std::find_if(it_b, containments.end(), [i](const containmentTuple &e) { return e.src > i; });
        containment_offsets[i+1] = std::distance(containments.begin(), it_e);
        it_b = it_e;
      }

      assert (it_b == containments.end());
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

    /**
     * print read ids which are contained but not redundant
     */
    void outputNonRedudantContainedReads (const std::string &filename)
    {
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

    //consider all contained reads as redundant and remove them from graph
    void removeContainedReads(const algoParams &param, std::ofstream &log)
    {
      for (uint32_t i = 0; i < readCount; i++) {
        if (deletedReads[i] == false && contained[i] == true) {
          deletedReads[i] = true;
          if (!param.logFileName.empty()) log << umap_inverse[i] << "\tremoveContainedReads()\n";
        }
      }

      std::vector<graphArc> edges_new;

      for (auto &e : edges)
      {
        //check if both end vertices of the edge are still available
        if (deletedReads[e.src >> 1] == false && deletedReads[e.dst >> 1] == false)
          edges_new.emplace_back(e);
      }

      edges = edges_new;
      std::cerr << "INFO, removeContainedReads() finished\n";
    }

    //consider contained reads as redundant if their
    //'containment degree' is > maxDegree
    void removeContainedReadsAboveDegree(const algoParams &param, std::ofstream &log)
    {
      for (uint32_t i = 0; i < readCount; i++) {
        if (getContaintmentDegree (i) > param.maxContainmentDegree && deletedReads[i] == false && contained[i] == true) {
          deletedReads[i] = true;
          if (!param.logFileName.empty()) log << umap_inverse[i] << "\tremoveContainedReadsAboveDegree()\n";
        }
      }


      std::vector<graphArc> edges_new;

      for (auto &e : edges)
      {
        //check if both end vertices of the edge are available
        if (deletedReads[e.src >> 1] == false && deletedReads[e.dst >> 1] == false)
          edges_new.emplace_back(e);
      }

      edges = edges_new;
      std::cerr << "INFO, removeContainedReadsAboveDegree() finished\n";
    }

    //algorithm motivated from Myers 2005
    //assumes indexing is done and edges are sorted
    //TODO: revise this for multi-graphs
    uint32_t transitiveReduction(int fuzz, std::ofstream& log)
    {
      if (fuzz < 0) return 0U; //reduction is disabled by user

      //mark 0 : default, 1 : in-play, 2 : reduced
      std::vector<uint8_t> mark (maxDegree(), 0);

      //save length of edge from vertex being considered
      uint32_t n_reduced = 0;

      for (uint32_t v = 0; v < vertexCount; v++)
      {
        if (getDegree (v) == 0) continue;

        //neighborhood of v
        for (uint32_t j = offsets[v]; j < offsets[v+1]; j++) {
          mark [j - offsets[v]] = 1; //edge in-play
        }

        uint32_t longest = edges[offsets[v+1] - 1].len + fuzz;

        //neighborhood of v
        for (uint32_t j = offsets[v]; j < offsets[v+1]; j++)
        {
          uint32_t w = edges[j].dst;
          if (mark [j - offsets[v]] != 1) continue;

          //neighborhood of w
          for (uint32_t k = offsets[w]; k < offsets[w+1]; k++)
          {
            uint32_t x = edges[k].dst;
            uint32_t sum = edges[j].len + edges[k].len; // v->w + w->x
            if (sum > longest) break;

            for (uint32_t l = offsets[v]; l < offsets[v+1]; l++)
            {
              if (edges[l].dst == x) { //v->x
                if (mark [l - offsets[v]] == 1) { //in-play?
                  if (sum <= edges[l].len + fuzz && sum + fuzz >= edges[l].len) { //length check
                    mark [l - offsets[v]] = 2; //eliminate this edge
                  }
                }
              }
            }
          }
        }

        for (uint32_t j = offsets[v]; j < offsets[v+1]; j++)
        {
          if (mark[j - offsets[v]] == 2)
            edges[j].del = true, n_reduced++;

          mark[j - offsets[v]] = 0; //edge is not in-play anymore
        }
      }

      std::cerr << "INFO, transitiveReduction(), " << n_reduced << " edges marked for deletion\n";

      /**
       * Ideally transitive reduction algorithm should preserve symmetry, but
       * it is not guaranteed if overlapper has <100% recall
       */
      ensureSymmetry (n_reduced);
      std::cerr << "INFO, transitiveReduction(), total " << n_reduced << " edges marked for deletion after ensuring symmetry\n";

      return n_reduced;
    }

    void ensureSymmetry (uint32_t &n_reduced)
    {
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
    bool checkSymmetry ()
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

//note: we assume that overlapper skipped dual mappings, minimap2 overlapping module skips by default
//note: substring overlaps which are not suffix-prefix are ignored, may need to relax this later
void ovlgraph_gen(const char *readfilename, const char *paffilename, float min_ovlp_identity, int min_ovlp_len, graphcontainer &g)
{
  //read input paf file
  paf_rec_t r;
  paf_file_t *fp = paf_open(paffilename);
  if (!fp) {
    std::cerr << "ERROR, ovlgraph_gen(), could not open PAF file\n";
    exit(1);
  }

  std::cerr << "INFO, ovlgraph_gen(), min_ovlp_identity = " << min_ovlp_identity << "\n";
  std::cerr << "INFO, ovlgraph_gen(), min_ovlp_len = " << min_ovlp_len << "\n";
  std::cerr << "INFO, ovlgraph_gen(), reading paf records from " << paffilename << "\n";

  uint64_t totPaf = 0, validPaf = 0, suffPrefPaf = 0, containedPaf = 0;
  while (paf_read(fp, &r) >= 0) {
    totPaf++;
    if (r.ml * 100.0 / r.bl >= min_ovlp_identity)
    {
      if (std::min (r.te - r.ts, r.qe - r.qs) >= min_ovlp_len)
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
        if (r.te - r.ts == r.tl) //target is contained in qry
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
  g.min_ovlp_len = min_ovlp_len;
  g.min_ovlp_identity = min_ovlp_identity;
  g.initVectors(readfilename);
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

#endif
