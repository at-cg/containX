#ifndef ASM_GRAPH_H
#define ASM_GRAPH_H

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>
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
    uint32_t dst_start_offset; //0-based
    uint32_t dst_end_offset;   //1-based
    uint32_t rev;              //src strand

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
      std::cerr << "INFO, initVectors(), " << std::count(contained.begin(), contained.end(), true) << " reads are contained in graph\n";

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
};

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
      int ovlp = std::max (r.te - r.ts, r.qe - r.qs);
      if (ovlp >= min_ovlp_len)
      {
        validPaf++;
        std::string qname = r.qn;
        std::string tname = r.tn;
        uint32_t q_vertexId = g.addStringToMap(qname);
        uint32_t t_vertexId = g.addStringToMap(tname);

        if (r.qe - r.qs == r.ql && r.te - r.ts != r.tl) //qry is contained in target
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
        else if (r.qe - r.qs != r.ql && r.te - r.ts == r.tl) //target is contained in qry
        {
          containmentTuple c;
          c.src = t_vertexId;
          c.dst = q_vertexId;
          c.rev = r.rev;
          if (c.rev) {
            c.dst_start_offset = r.ql - r.qe;
            c.dst_end_offset = r.ql - r.qs;
          }
          else
          {
            c.dst_start_offset = r.qs;
            c.dst_end_offset = r.qe;
          }
          g.containments.emplace_back(c);
          containedPaf++;
        }
        //suffix of qry overlaps with prefix of target
        else if (r.qs > 0 && r.qe == r.ql && r.ts == 0 && r.te < r.tl) 
        {
          graphArc e1, e2;

          {
            e1.src = q_vertexId;
            e1.dst = t_vertexId;
            e1.rev_src = r.rev;
            e1.rev_dst = 0;
            e1.ov_src = r.qe - r.qs;
            e1.ov_dst = r.te - r.ts;
          }

          {
            e2.src = t_vertexId;
            e2.dst = q_vertexId;
            e2.rev_src = 1;
            e2.rev_dst = 1 - r.rev;
            e2.ov_src = r.te - r.ts;
            e2.ov_dst = r.qe - r.qs;
          }

          g.edges.emplace_back(e1);
          g.edges.emplace_back(e2);
          suffPrefPaf++;
          
        }
        //prefix of qry overlaps with suffix of target
        else if (r.qs == 0 && r.qe < r.ql && r.ts > 0 && r.te == r.tl)
        {
          graphArc e1, e2;

          {
            e1.src = t_vertexId;
            e1.dst = q_vertexId;
            e1.rev_src = 0;
            e1.rev_dst = r.rev;
            e1.ov_src = r.te - r.ts;
            e1.ov_dst = r.qe - r.qs;
          }

          {
            e2.src = q_vertexId;
            e2.dst = t_vertexId;
            e2.rev_src = 1;
            e2.rev_dst = 1 - r.rev;
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

  g.initVectors(readfilename);

}

#endif
