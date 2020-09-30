#include <vector>
#include <map>

#include "SeqLib/BWAWrapper.h"

#include "ReadTable.h"
#include "OverlapCommon.h" //#include "SuffixArray.h" #include "Util.h"
#include "SGUtil.h" // String Graph

using namespace std;

namespace sv_helper{
  class InMemBwa{
    public:
      InMemBwa(char *ref, bool is_vir);
      ~InMemBwa();
      int is_adjacent(char *seq, int chr_id, int pos, int dist, bool is_reverse_retricted, bool is_reverse,
                      int &fetch_s, int &fetch_e, vector<string> &fetch_cg, int &mapq);
      void bwa_with_contig(const string &contig, const vector<string> &related_refs, const vector<string> &sc_seqs, const vector<string> &su_seqs,
                           int bp_in_contig, bool is_contig_left, int &sc_count, int &su_count);
      void run_get_contig_supt(const string &contig, const vector<string> &related_refs, const map<string, vector<string>> &q2seqs, vector<string> &supt_qnames,
                               int bp1_in_contig, int bp2_in_contig, int flank, int min_map2ctg_qual);
      int align(char *seq, string &aln_seq, int &chr_id, int &pos, int &is_reverse, vector<string> &cg, int &mapq);
    
    private:
      SeqLib::BWAWrapper *bwa;
  };
  
  class LocalSGA{
    public:
      LocalSGA(const string& id, double er, size_t mo, size_t rl) : m_id(id), m_error_rate(er), m_min_overlap(mo), m_readlen(rl) {}
      ~LocalSGA() {}
      void run_asm(vector<string>& seqs, int num_assembly_rounds, vector<string*> &ret_seqs);
      
      void doAssembly(ReadTable *pRT, SeqLib::UnalignedSequenceVector &contigs, int pass);
      void calculateSeedParameters(int read_len, const int minOverlap, int& seed_length, int& seed_stride) const;
      StringGraph* assemble(stringstream& asqg_stream, int minOverlap, int maxEdges, bool bExact, 
                int trimLengthThreshold, bool bPerformTR, bool bValidate, int numTrimRounds, 
                int resolveSmallRepeatLen, int numBubbleRounds, double maxBubbleGapDivergence, 
                double maxBubbleDivergence, int maxIndelLength, int cutoff, string prefix, 
                SeqLib::UnalignedSequenceVector &contigs, bool get_components);
      void remove_exact_dups(SeqLib::UnalignedSequenceVector& cc) const;
      void setToWriteASQG() { m_write_asqg = true; }
      ReadTable* removeDuplicates(ReadTable* pRT);
    
    private:
      string m_id;
      double m_error_rate;
      size_t m_min_overlap;
      int m_readlen;
      
      int maxEdges = 128;
      bool bExact = true;
      bool bPerformTR = false; // transitivie edge reducetion
      bool bValidate = false;
      int resolveSmallRepeatLen = -1;
      size_t numBubbleRounds = 3;
      float gap_divergence = 0.00;
      float divergence = 0.00; //0.05
      int maxIndelLength = 20;
      
      bool m_write_asqg = false;
  };
  
  
  struct AssemblyOptions {
    
    unsigned int verbose = 1;
    string asqgFile;
    string outContigsFile;
    string outVariantsFile;
    string outGraphFile;
    
    unsigned int minOverlap;
    bool bEdgeStats = false;
    bool bSmoothGraph = false;
    int resolveSmallRepeatLen = -1;
    size_t maxEdges = 128;
    
    // Trim parameters
    int numTrimRounds = 10;
    size_t trimLengthThreshold = 300;
    
    // Bubble parameters
    int numBubbleRounds = 3;
    double maxBubbleDivergence = 0.05f;
    double maxBubbleGapDivergence = 0.01f;
    int maxIndelLength = 20;
    
    // 
    bool bValidate;
    bool bExact = true;
    bool bPerformTR = false;
  };
}