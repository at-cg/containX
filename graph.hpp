#ifndef ASM_GRAPH_H
#define ASM_GRAPH_H

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <fstream>
#include "paf.hpp"

class graphArc
{
  public:
    uint32_t src;
    uint32_t dst;
    uint32_t rev_src:1, ov_src:31; //strand, length of suffix of src involved in overlap
    uint32_t rev_dst:1, ov_dst:31; //strand, length of prefix of dst involved in overlap
};

class containmentTuple
{
  public:
    uint32_t src;              //read which is contained
    uint32_t dst;
    uint32_t dst_start_offset; //(0-based; BED-like; closed) leftmost pos in forward orientation
    uint32_t dst_end_offset;   //(0-based; BED-like; open)
    uint32_t rev;              //dst strand

    bool operator< (const containmentTuple &r )
    {
      return std::tie(src, rev) < std::tie(r.src, r.rev);
    }
};

class graphcontainer
{
  public:
    std::unordered_map <std::string, uint32_t> umap;  // size = count of reads
    std::vector<std::string> readseq;  //size = count of reads
    std::vector<bool> contained; //size = count of reads
    std::vector<bool> redundant;   //size = count of reads
    std::vector<graphArc> edges; //size = 2 x suffix-prefix overlaps 
    std::vector<containmentTuple> containments; //size = count of contained overlaps 

    //save user thresholds
    float min_ovlp_identity;
    int min_ovlp_len;

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

    void initVectors(const char *readfilename)
    {
      assert (contained.size() == 0);
      assert (redundant.size() == 0);
      assert (readseq.size() == 0);

      int vertexCount = umap.size();
      contained.resize(vertexCount, false);
      redundant.resize(vertexCount, false);
      readseq.resize(vertexCount, "");

      std::sort(containments.begin(), containments.end());

      for (auto &e : containments)
        contained[e.src] = true;


      std::cerr << "INFO, initVectors(), graph has " << vertexCount << " vertices in total\n";
      std::cerr << "INFO, initVectors(), " << std::count(contained.begin(), contained.end(), true) << " reads are marked as contained in graph\n";

      //parse reads
      {
        std::cerr << "INFO, initVectors(), parsing reads from " << readfilename << "\n";
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
    }

    void outputGFA (const std::string &filename) const
    {
      std::ofstream outstrm (filename);

      //print reads - use prefix utg before their id
      for (uint32_t i = 0; i < umap.size(); i++)
        outstrm  << "S\tutg" << i << "\t" << readseq[i] << "\n";

      //print arcs
      for (uint32_t i = 0; i < edges.size(); i++)
        outstrm << "L\tutg" << edges[i].src << "\t" << "+-"[edges[i].rev_src] <<  "\tutg" << edges[i].dst << "\t" << "+-"[edges[i].rev_dst] << "\t" << edges[i].ov_src <<  "M\n"; 

      /**
       * URL: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
       * A link from utg1 to utg2 means that the end of utg1 overlaps with the start of utg2
       *     L <utg1> <utg1-orient> <utg2> <utg2-orient> <Overlap>
       */

      //print containment lines
      for (uint32_t i = 0; i < containments.size(); i++)
        outstrm << "C\tutg" << containments[i].dst << "\t" << "+-"[containments[i].rev] << "\tutg" << containments[i].src << "\t+\t" << containments[i].dst_start_offset << "\t" << containments[i].dst_end_offset - containments[i].dst_start_offset<< "M\n"; 

      /**
       * Example (from GFA1 spec)
       * The following line describes the containment of segment utg2 in the reverse 
       * complement of segment utg1, starting at position 110 of segment utg1 (in its forward orientation).
       *      C  utg1 - utg2 + 110 100M
       */
    }
};

#warning "we assume that overlapper skipped dual mappings, minimap2 overlapping module skips by default"
#warning "substring overlaps which are not suffix-prefix are ignored"
void ovlgraph_gen(const char *readfilename, const char *paffilename, float min_ovlp_identity, int min_ovlp_len, graphcontainer &g)
{
  //read input paf file
  paf_rec_t r;
  paf_file_t *fp = paf_open(paffilename);
  if (!fp) {
    std::cerr << "ERROR, ovlgraph_gen(), could not open PAF file\n";
    exit(1);
  }

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
        uint32_t q_vertexId = g.addStringToMap(qname);
        uint32_t t_vertexId = g.addStringToMap(tname);

        if (q_vertexId == t_vertexId) continue;

        if (r.qe - r.qs == r.ql) //qry is contained in target
        {
          containmentTuple c;
          c.src = q_vertexId;
          c.dst = t_vertexId;
          c.rev = r.rev;  
          c.dst_start_offset = r.ts;
          c.dst_end_offset = r.te;
          g.containments.emplace_back(c);
          containedPaf++;
        }
        if (r.te - r.ts == r.tl) //target is contained in qry
        {
          containmentTuple c;
          c.src = t_vertexId;
          c.dst = q_vertexId;
          c.rev = r.rev;
          c.dst_start_offset = r.qs;
          c.dst_end_offset = r.qe;
          g.containments.emplace_back(c);
          containedPaf++;
        }

        /**
         * Conditions below check for suffix-prefix overlaps
         * See Insight 13.2 Reverse complements in assembly graphs in MBCT textbook for details
         */

        //a suffix of qry overlaps a prefix of target
        if (r.rev == 0 && r.qs > 0 && r.qe == r.ql && r.ts == 0 && r.te < r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = q_vertexId;
            e1.dst = t_vertexId;
            e1.rev_src = 0;
            e1.rev_dst = 0;
            e1.ov_src = r.qe - r.qs;
            e1.ov_dst = r.te - r.ts;
          }

          {
            e2.src = t_vertexId;
            e2.dst = q_vertexId;
            e2.rev_src = 1;
            e2.rev_dst = 1;
            e2.ov_src = r.te - r.ts;
            e2.ov_dst = r.qe - r.qs;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
        }
        //a prefix of qry overlaps a suffix of target
        if (r.rev == 0 && r.qs == 0 && r.qe < r.ql && r.ts > 0 && r.te == r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = t_vertexId;
            e1.dst = q_vertexId;
            e1.rev_src = 0;
            e1.rev_dst = 0;
            e1.ov_src = r.te - r.ts;
            e1.ov_dst = r.qe - r.qs;
          }

          {
            e2.src = q_vertexId;
            e2.dst = t_vertexId;
            e2.rev_src = 1;
            e2.rev_dst = 1;
            e2.ov_src = r.qe - r.qs;
            e2.ov_dst = r.te - r.ts;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
        }
        //a suffix of ~qry overlaps a prefix of target
        if (r.rev == 1 && r.qs == 0 && r.qe < r.ql && r.ts == 0 && r.te < r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = q_vertexId;
            e1.dst = t_vertexId;
            e1.rev_src = 1;
            e1.rev_dst = 0;
            e1.ov_src = r.qe - r.qs;
            e1.ov_dst = r.te - r.ts;
          }

          {
            e2.src = t_vertexId;
            e2.dst = q_vertexId;
            e2.rev_src = 1;
            e2.rev_dst = 0;
            e2.ov_src = r.te - r.ts;
            e2.ov_dst = r.qe - r.qs;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
        }
        //a prefix of ~qry overlaps a suffix of target
        if (r.rev == 1 && r.qs > 0 && r.qe == r.ql && r.ts > 0 && r.te == r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = t_vertexId;
            e1.dst = q_vertexId;
            e1.rev_src = 0;
            e1.rev_dst = 1;
            e1.ov_src = r.te - r.ts;
            e1.ov_dst = r.qe - r.qs;
          }

          {
            e2.src = q_vertexId;
            e2.dst = t_vertexId;
            e2.rev_src = 0;
            e2.rev_dst = 1;
            e2.ov_src = r.qe - r.qs;
            e2.ov_dst = r.te - r.ts;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
        }
      }
    }
  }

  paf_close(fp);

  std::cerr << "INFO, ovlgraph_gen(), parsed " << totPaf << " paf records\n";
  std::cerr << "INFO, ovlgraph_gen(), " << validPaf << " records satisfied user-specified cutoffs\n";
  std::cerr << "INFO, ovlgraph_gen(), " << containedPaf << " records belonged to contained overlaps\n";
  std::cerr << "INFO, ovlgraph_gen(), " << suffPrefPaf << " records belonged to proper suffix-prefix overlaps\n";

  g.min_ovlp_len = min_ovlp_len;
  g.min_ovlp_identity = min_ovlp_identity;

  g.initVectors(readfilename);

}

#endif
