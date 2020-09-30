#include <iostream>

#include "sv_helper.h"
#include "SeqLib/BamRecord.h"

#include "SeqLib/UnalignedSequence.h"
#include "svabaOverlapAlgorithm.h"
#include "svabaASQG.h"
#include "SGVisitors.h"

#define MAX_OVERLAPS_PER_ASSEMBLY 20000

using namespace std;

namespace sv_helper{
  InMemBwa::InMemBwa(char *ref, bool is_vir){
    bwa = new SeqLib::BWAWrapper();
    
    if (!is_vir){
      // bwa->SetAScore(2);
      // bwa->SetGapOpen(32);
      // bwa->SetGapExtension(1);
      // bwa->SetMismatchPenalty(18);
      // bwa->SetZDropoff(100);
      // bwa->SetBandwidth(1000);
      // bwa->SetReseedTrigger(1.5);
      // bwa->Set3primeClippingPenalty(5);
      // bwa->Set5primeClippingPenalty(5);
    }
    bwa->LoadIndex(ref);
  }
  
  
  InMemBwa::~InMemBwa(){
    delete bwa;
    bwa = nullptr;
  }
  
  
  int InMemBwa::is_adjacent(char *seq, int chr_id, int pos, int dist, bool is_reverse_retricted, bool is_reverse,
                            int &fetch_s, int &fetch_e, vector<string> &fetch_cg, int &mapq){
    SeqLib::BamRecordVector results;
    bool hardclip = false;
    float secondary_cutoff = 0.6; // secondary alignments must have score >= 0.6*top_score
    int secondary_cap = 20; // max number of secondary alignments to return
    
    bwa->AlignSequence(seq, "test", results, hardclip, secondary_cutoff, secondary_cap);
    if (!results.size()){
      // cerr << "Not mapped: " << seq << endl;
      return 0;
    }
    SeqLib::GenomicRegion grm(chr_id, max(0, pos-dist), pos+dist);
    for (auto& i : results){
      if (grm.GetOverlap(i.AsGenomicRegion()) && (!is_reverse_retricted || i.ReverseFlag() == is_reverse)){
        fetch_s = i.AsGenomicRegion().pos1;
        fetch_e = i.AsGenomicRegion().pos2;
        mapq = i.MapQuality();
        for (auto& c : i.GetCigar()){
          char ctp = c.Type();
          string coplen = to_string(c.Length());
          string ctype(1, ctp);
          fetch_cg.push_back(coplen);
          fetch_cg.push_back(ctype);
        }
        return 1;
      }
    }
    return 0;
  }
  
  int InMemBwa::align(char *seq, string &aln_seq, int &chr_id, int &pos, int &is_reverse, vector<string> &cg, int &mapq){
    SeqLib::BamRecordVector results;
    bool hardclip = false;
    float secondary_cutoff = 0.9; // for vir, use normal value
    int secondary_cap = 10;
    
    bwa->AlignSequence(seq, "test", results, hardclip, secondary_cutoff, secondary_cap);
    for (auto& i : results){
      chr_id = i.AsGenomicRegion().chr;
      pos = i.AsGenomicRegion().pos1;
      is_reverse = i.ReverseFlag();
      mapq = i.MapQuality();
      aln_seq = i.Sequence();
      for (auto& c : i.GetCigar()){
        char ctp = c.Type();
        string coplen = to_string(c.Length());
        string ctype(1, ctp);
        cg.push_back(coplen);
        cg.push_back(ctype);
      }
      return 1;
    }
    return 0;
  }
  
  
  int count_contig_support(const SeqLib::BWAWrapper &bw_ref, const SeqLib::BWAWrapper &bw_c, const string &seq, int bp_in_contig, bool is_contig_left){
    SeqLib::BamRecordVector brv_ref, brv_c;
    bool hardclip = false;
    
    bw_c.AlignSequence(seq, "2c", brv_c, hardclip, 0.60, 10000);
    if (brv_c.size() == 0) 
      return 0;
    
    // get the maximum non-reference alignment score
    int max_as = 0;
    for (auto& r : brv_c) {
      int thisas = 0;
      r.GetIntTag("AS", thisas);
      max_as = max(max_as, thisas);
    }
    
    bw_ref.AlignSequence(seq, "2r", brv_ref, hardclip, 0.60, 10);
    
    // get the maximum reference alignment score
    int max_as_r = 0;
    for (auto& r : brv_ref) {
      int thisas= 0;
      r.GetIntTag("AS", thisas);
      max_as_r = max(max_as_r, thisas);
    }
    
    // reject if better alignment to reference
    if (max_as_r > max_as) return 0;
    
    for (auto &r : brv_c){
      // make sure alignment score is OK
      int thisas = 0;
      r.GetIntTag("AS", thisas);
      if ((double)r.NumMatchBases() * 0.5 > thisas) /* && i.GetZTag("SR").at(0) == 't'*/
        continue;
      
      int r_end=r.PositionEnd(), r_start=r.Position();
      if ((r_end - r_start) < ((double)seq.length() * 0.75))
        continue;
      
      // r_start, r_end is bed format, bp_in_contig is 1-based
      if ((is_contig_left && bp_in_contig >= r_end) || (!is_contig_left && bp_in_contig <= r_start+1))
        continue;
      
      return 1;
    }
    return 0;
  }
  
  
  void InMemBwa::bwa_with_contig(const string &contig, const vector<string> &related_refs, const vector<string> &sc_seqs, const vector<string> &su_seqs,
                                 int bp_in_contig, bool is_contig_left, int &sc_count, int &su_count){
    // make the reference allele BWAWrapper
    SeqLib::BWAWrapper bw_ref, bw_c;
    SeqLib::UnalignedSequenceVector usv_ref, usv_c;
    
    usv_c.push_back({"contig", contig, string()}); // name, seq, qual
    int rind = 0;
    for (auto s : related_refs){
      usv_ref.push_back({to_string(rind++), s, string()});
    }
    bw_c.ConstructIndex(usv_c);
    bw_ref.ConstructIndex(usv_ref);
    
    // // set up custom alignment parameters, mean
    // bw_ref.SetGapOpen(16); // default 6
    // bw_ref.SetMismatchPenalty(9); // default 2
    // bw_c.SetGapOpen(16); // default 6
    // bw_c.SetMismatchPenalty(9); // default 4
    
    for (auto s : sc_seqs) {
      sc_count += count_contig_support(bw_ref, bw_c, s, bp_in_contig, is_contig_left);
    }
    for (auto s : su_seqs) {
      su_count += count_contig_support(bw_ref, bw_c, s, bp_in_contig, is_contig_left);
    }
  }
  

  void add_support_qname(const SeqLib::BWAWrapper &bw_ref, const SeqLib::BWAWrapper &bw_c, const string &qname, const vector<string> &seqs,
                         int bp1_in_contig, int bp2_in_contig, vector<string> &supt_qnames, int flank, int min_map2ctg_qual){
    SeqLib::BamRecordVector brv_ref, brv_c;
    bool hardclip = false;

    for (auto &s : seqs){
      bw_c.AlignSequence(s, "2contig", brv_c, hardclip, 0.60, 20);
      if (brv_c.size() == 0) 
        continue;
      // get the maximum non-reference alignment score
      int max_as = 0;
      for (auto &r : brv_c) {
        int thisas = 0;
        r.GetIntTag("AS", thisas);
        max_as = max(max_as, thisas);
      }
      if (max_as < min_map2ctg_qual) continue;
      
      bw_ref.AlignSequence(s, "2ref", brv_ref, hardclip, 0.60, 10);
      // get the maximum reference alignment score
      int max_as_r = 0;
      for (auto &r : brv_ref) {
        int thisas= 0;
        r.GetIntTag("AS", thisas);
        max_as_r = max(max_as_r, thisas);
      }
      
      // reject if better alignment to reference
      if (max_as_r > max_as) continue;
      
      for (auto &r : brv_c){
        // make sure alignment score is OK
        int thisas = 0;
        r.GetIntTag("AS", thisas);
        if ((double)r.NumMatchBases() * 0.5 > thisas) /* && i.GetZTag("SR").at(0) == 't'*/
          continue;
        
        int r_end=r.PositionEnd(), r_start=r.Position();
        if ((r_end - r_start) < ((double)s.length() * 0.75))
          continue;
        
        // r_start, r_end is bed format, bp_in_contig is 1-based
        if ((r_start > bp2_in_contig-flank) || (r_end < bp1_in_contig+flank))
          continue;
        supt_qnames.push_back(qname);
        return;
      }
      brv_c.clear();
      brv_ref.clear();
    }
  }


  void InMemBwa::run_get_contig_supt(const string &contig, const vector<string> &related_refs, const map<string, vector<string>> &q2seqs, vector<string> &supt_qnames,
                                     int bp1_in_contig, int bp2_in_contig, int flank, int min_map2ctg_qual){
    // make the reference allele BWAWrapper
    SeqLib::BWAWrapper bw_ref, bw_c;
    SeqLib::UnalignedSequenceVector usv_ref, usv_c;
    
    usv_c.push_back({"contig", contig, string()}); // name, seq, qual
    int rind = 0;
    for (auto s : related_refs){
      usv_ref.push_back({to_string(rind++), s, string()});
    }
    bw_c.ConstructIndex(usv_c);
    bw_ref.ConstructIndex(usv_ref);
    
    // // set up custom alignment parameters, mean
    // bw_ref.SetGapOpen(16); // default 6
    // bw_ref.SetMismatchPenalty(9); // default 2
    // bw_c.SetGapOpen(16); // default 6
    // bw_c.SetMismatchPenalty(9); // default 4
    
    map<string, vector<string>>::const_iterator iter;
    for (iter = q2seqs.begin(); iter != q2seqs.end(); iter++) {
      add_support_qname(bw_ref, bw_c, iter->first, iter->second, bp1_in_contig, bp2_in_contig, supt_qnames, flank, min_map2ctg_qual);
    }     
  }

  
  void LocalSGA::run_asm(vector<string>& seqs, int num_assembly_rounds, vector<string*> &ret_seqs){
    ReadTable m_pRT;
    
    int count = 0;
    //! seqs shouldnt contain N
    for (auto& s : seqs) {
      SeqItem si;
      si.id = "read_" + to_string(++count); 
      si.seq = s;
      
      m_pRT.addRead(si);
    }
    
    SeqLib::UnalignedSequenceVector m_contigs;
    doAssembly(&m_pRT, m_contigs, 0);
    int before_size = -1;
    for (int yy = 1; yy != num_assembly_rounds; yy++) {
      if (m_contigs.size() < 2) break; // break because too few contigs to assemle
      //C check contig num, break if not change
      if (m_contigs.size() == before_size) break;
      before_size = m_contigs.size();
      
      // do the second round (on assembled contigs)
      SeqLib::UnalignedSequenceVector tmpc;
      for (auto& j : m_contigs)
        if (j.Seq.length() > (uint) m_readlen)
          tmpc.push_back(j);
      
      ReadTable pRTc0(tmpc);
      m_contigs.clear();
      doAssembly(&pRTc0, m_contigs, yy);
    }
    
    for (auto& j : m_contigs){
      string *jstr = new string(j.Seq);
      ret_seqs.push_back(jstr);
    }
    return;
  }
  
  
  void LocalSGA::doAssembly(ReadTable *pRT, SeqLib::UnalignedSequenceVector &contigs, int pass) {
    if (pRT->getCount() == 0) return;
    
    // clear the hits stream
    stringstream hits_stream, asqg_stream;

    // set the paramters for this run
    double errorRate = m_error_rate;
    int min_overlap = m_min_overlap;
    int cutoff = 0;

    if (pass > 1) {
    //errorRate = 0.015f; // up to one bp mismatch
      errorRate = 0;
      m_error_rate = errorRate;
      min_overlap = m_min_overlap * 2;
    }

    bool exact = errorRate < 0.001f;

    ReadTable * pRT_nd = exact ? removeDuplicates(pRT) : pRT;

    // forward
    SuffixArray* pSAf_nd = new SuffixArray(pRT_nd, 1, false); //1 is num threads. false is silent/no
    RLBWT *pBWT_nd = new RLBWT(pSAf_nd, pRT_nd);
    
    // reverse
    pRT_nd->reverseAll();
    SuffixArray * pSAr_nd = new SuffixArray(pRT_nd, 1, false);
    RLBWT *pRBWT_nd = new RLBWT(pSAr_nd, pRT_nd);
    pRT_nd->reverseAll();
    
    pSAf_nd->writeIndex();
    pSAr_nd->writeIndex();
    
    bool bIrreducibleOnly = true; // default
    int seedLength = 0;
    int seedStride = 0;
    
    if (!exact) calculateSeedParameters(m_readlen, min_overlap, seedLength, seedStride);
    
    svabaOverlapAlgorithm* pOverlapper = new svabaOverlapAlgorithm(pBWT_nd, pRBWT_nd, 
                     errorRate, seedLength,
                     seedStride, bIrreducibleOnly);
    
    pOverlapper->setExactModeOverlap(exact);
    pOverlapper->setExactModeIrreducible(exact);

    pRT_nd->setZero();

    size_t workid = 0;
    SeqItem si;
    
    size_t ocount = 0;
    while (pRT_nd->getRead(si) && (++ocount < MAX_OVERLAPS_PER_ASSEMBLY)) {
      SeqRecord read;
      read.id = si.id;
      read.seq = si.seq;
      OverlapBlockList obl;
      
      OverlapResult rr = pOverlapper->overlapRead(read, min_overlap, &obl);
      pOverlapper->writeOverlapBlocks(hits_stream, workid, rr.isSubstring, &obl);

      svabaASQG::VertexRecord record(read.id, read.seq.toString());
      record.setSubstringTag(rr.isSubstring);
      record.write(asqg_stream);
      
      ++workid;
    }
  
    string line;
    bool bIsSelfCompare = true;
    ReadInfoTable* pQueryRIT = new ReadInfoTable(pRT_nd);

    while(getline(hits_stream, line)) {
      size_t readIdx;
      size_t totalEntries;
      bool isSubstring; 
      OverlapVector ov;
      OverlapCommon::parseHitsString(line, pQueryRIT, pQueryRIT, pSAf_nd, pSAr_nd, bIsSelfCompare, readIdx, totalEntries, ov, isSubstring);

      for(OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter) {
         svabaASQG::EdgeRecord edgeRecord(*iter);
         edgeRecord.write(asqg_stream);
      }
    }
    
    // Get the number of strings in the BWT, this is used to pre-allocated the read table
    delete pOverlapper;
    delete pBWT_nd; 
    delete pRBWT_nd;
    delete pSAf_nd;
    delete pSAr_nd;
    if (exact) delete pRT_nd;
    
    // PERFORM THE ASSMEBLY
    int trimLengthThreshold = 100; 
    int numTrimRounds = 1; 
    StringGraph * oGraph = assemble(asqg_stream, min_overlap, maxEdges, bExact, 
       trimLengthThreshold, bPerformTR, bValidate, numTrimRounds, 
       resolveSmallRepeatLen, numBubbleRounds, gap_divergence, 
       divergence, maxIndelLength, cutoff, m_id + "_", contigs, m_write_asqg);
    
    // optionally output the graph structure
    // if (m_write_asqg)
      // write_asqg(oGraph, asqg_stream, hits_stream, pass);

    // this was allocated in assemble
    delete oGraph;
    delete pQueryRIT;
    
    // remove exact dups
    remove_exact_dups(contigs);
    
    return;
  }
  
  void LocalSGA::remove_exact_dups(SeqLib::UnalignedSequenceVector& cc) const {
  
    set<string> ContigDeDup;
    SeqLib::UnalignedSequenceVector cvec;
    for (auto& i : cc) {
      if (!ContigDeDup.count(i.Seq)) {
        ContigDeDup.insert(i.Seq);
        cvec.push_back(i);
      } else {
        //cerr << "Filtered out a contig for having exact duplicate with another contig" << endl;
      }
    }
    cc = cvec;
  }
  
  void LocalSGA::calculateSeedParameters(int read_len, const int minOverlap, int& seed_length, int& seed_stride) const{
    seed_length = 0;
    
    // The maximum possible number of differences occurs for a fully-aligned read
    int max_diff_high = static_cast<int>(m_error_rate * read_len);

    // Calculate the seed length to use
    // If the error rate is so low that no differences are possible just seed
    // over the entire minOverlap region
    if(max_diff_high > 0)
    {
        // Calculate the maximum number of differences between two sequences that overlap
        // by minOverlap
        int max_diff_low = static_cast<int>(m_error_rate * minOverlap);

         if(max_diff_low == 0)
            max_diff_low = 1;
         
         int seed_region_length = static_cast<int>(ceil(max_diff_low / m_error_rate));
         int num_seeds_low = max_diff_low + 1;
         seed_length = static_cast<int>(seed_region_length / num_seeds_low);
         if(seed_length > static_cast<int>(minOverlap))
            seed_length = minOverlap;
    }
    else
    {
        seed_length = minOverlap;
    }
    seed_stride = seed_length;
  }
  
  StringGraph* LocalSGA::assemble(stringstream& asqg_stream, int minOverlap, int maxEdges, bool bExact, 
          int trimLengthThreshold, bool bPerformTR, bool bValidate, int numTrimRounds, 
          int resolveSmallRepeatLen, int numBubbleRounds, double maxBubbleGapDivergence, 
          double maxBubbleDivergence, int maxIndelLength, int cutoff, string prefix, 
          SeqLib::UnalignedSequenceVector &contigs, bool get_components){

    AssemblyOptions ao;

    StringGraph * pGraph = SGUtil::loadASQG(asqg_stream, minOverlap, true, maxEdges);
    pGraph->m_get_components = get_components;

    if(bExact)
      pGraph->setExactMode(true);
    
    // Visitor functors
    SGTransitiveReductionVisitor trVisit;
    SGGraphStatsVisitor statsVisit;
    SGTrimVisitor trimVisit(trimLengthThreshold);
    SGContainRemoveVisitor containVisit;
    SGValidateStructureVisitor validationVisit;
    
    while(pGraph->hasContainment())
      pGraph->visit(containVisit);

    // Remove any extraneous transitive edges that may remain in the graph
    if(bPerformTR)
      {
        cout << "Removing transitive edges\n";
        pGraph->visit(trVisit);
      }
    
    // Compact together unbranched chains of vertices
    pGraph->simplify();
    
    if(bValidate)
      {
        pGraph->visit(validationVisit);
      }

    // Remove dead-end branches from the graph
    if(numTrimRounds > 0) {
        int numTrims = numTrimRounds;
        while(numTrims-- > 0)
    pGraph->visit(trimVisit);
      }
    
    // Resolve small repeats
    if(resolveSmallRepeatLen > 0 && false) {
        SGSmallRepeatResolveVisitor smallRepeatVisit(resolveSmallRepeatLen);
      }
    
    if(numBubbleRounds > 0)
      {
        SGSmoothingVisitor smoothingVisit(ao.outVariantsFile, maxBubbleGapDivergence, maxBubbleDivergence, maxIndelLength);
        int numSmooth = numBubbleRounds;
        while(numSmooth-- > 0)
    pGraph->visit(smoothingVisit);
        pGraph->simplify(); 
      }

    pGraph->renameVertices(prefix);

    SGVisitorContig av;
    pGraph->visit(av);
    
    SeqLib::UnalignedSequenceVector tmp = av.m_ct;
    for (SeqLib::UnalignedSequenceVector::const_iterator it = tmp.begin(); it != tmp.end(); it++) {
      if ((int)(it->Seq.length()) >= cutoff) {
        contigs.push_back({it->Name + "C", it->Seq, string()}); //postpend with character to distribugish _2 from _22
      }
    }
    
    /*  if (walk_all) {
        SGWalkVector outWalks;
        walkExtra(pGraph, outWalks);
        for (auto& i : outWalks) {
        string seqr = i.getString(SGWT_START_TO_END);
        if ((int)seqr.length() >= cutoff) 
        contigs.push_back(Contig(i.pathSignature(), seqr));
        }
        }*/
    

    return pGraph;
    //    delete pGraph;
     
  }
  
  
  // not totally sure this works...
  ReadTable* LocalSGA::removeDuplicates(ReadTable* pRT) {
    // forward
    SuffixArray* pSAf = new SuffixArray(pRT, 1, false); //1 is num threads. false is silent/no
    RLBWT *pBWT= new RLBWT(pSAf, pRT);

    // reverse
    pRT->reverseAll();
    SuffixArray * pSAr = new SuffixArray(pRT, 1, false);
    RLBWT *pRBWT = new RLBWT(pSAr, pRT);
    pRT->reverseAll();

    svabaOverlapAlgorithm* pRmDupOverlapper = new svabaOverlapAlgorithm(pBWT, pRBWT, 
                                        0, 0, 
                                        0, false);
    
    pRT->setZero();
    ReadTable * pRT_nd = new ReadTable();
    SeqItem sir;
    while (pRT->getRead(sir)) {
      OverlapBlockList OBout;
      SeqRecord read;
      read.id = sir.id;
      read.seq = sir.seq;
      OverlapBlockList obl;
      OverlapResult rr = pRmDupOverlapper->alignReadDuplicate(read, &OBout);

      if (!rr.isSubstring)
        pRT_nd->addRead(sir);
    }

    delete pRmDupOverlapper;
    delete pBWT; 
    delete pRBWT;
    delete pSAf;
    delete pSAr;

    return pRT_nd;
  }
}