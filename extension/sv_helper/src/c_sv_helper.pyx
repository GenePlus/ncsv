from libc.stdlib cimport malloc
from cpython.string cimport PyString_AsString
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from cython.operator cimport dereference as deref

from cpython cimport bool as bool_p
from libcpp cimport bool as bool_c

cdef extern from "sv_helper.h" namespace "sv_helper":
    cdef cppclass InMemBwa:
        InMemBwa(char *, bool_c) except +
        bint is_adjacent(char *seq, int chr_id, int pos, int dist, bool_c is_reverse_retricted, bool_c is_reverse, int &fetch_s, int &fetch_e, vector[string] &fetch_cg, int &mapq)
        void bwa_with_contig(const string &contig, const vector[string] &related_refs, const vector[string] &sc_seqs,
                             const vector[string] &su_seqs, int bp_in_contig, bool_c is_contig_left, int &sc_count, int &su_count)
        void run_get_contig_supt(const string &contig, const vector[string] &related_refs, const map[string, vector[string]] &q2seqs,
                                 vector[string] &supt_qnames, int bp1_in_contig, int bp2_in_contig, int flank, int min_map2ctg_qual)
        bint align(char *seq, string &aln_seq, int &chr_id, int &pos, int &is_reverse, vector[string] &cg, int &mapq)
    
    cdef cppclass LocalSGA:
        LocalSGA(char *, double, size_t, size_t) except +
        void run_asm(vector[string] &seqs, int num_assembly_rounds, vector[string*] &ret_seqs)

cdef class InMemBwaWrapper:
    cdef InMemBwa *c_bwa
    
    def __cinit__(self, char *ref, bool_p is_vir):
        self.c_bwa = new InMemBwa(ref, is_vir)
    
    def run(self, char *seq, int chr_id, int pos, int dist, bool_p is_reverse_retricted, bool_p is_reverse):
        cdef int fetch_s=0, fetch_e=0, mapq=0
        cdef vector[string] fetch_cg
        if self.c_bwa.is_adjacent(seq, chr_id, pos, dist, is_reverse_retricted, is_reverse,
                                  fetch_s, fetch_e, fetch_cg, mapq):
            cg_tuples = [(cg, int(cg_len)) for cg_len, cg in zip(fetch_cg[::2], fetch_cg[1::2])]
            return fetch_s, fetch_e, cg_tuples, mapq
        return False
    
    def align(self, char *seq):
        cdef int chr_id=0, pos=0, mapq=0, is_reverse=0
        cdef vector[string] cg
        cdef string aln_seq
        if self.c_bwa.align(seq, aln_seq, chr_id, pos, is_reverse, cg, mapq):
            return aln_seq, chr_id, pos, cg, mapq, is_reverse>0
        return False
    
    def run_with_contig(self, string contig, vector[string] related_refs, vector[string] sc_seqs, vector[string] su_seqs, int bp_in_contig, bool_p is_contig_left):
        cdef int sc_count=0, su_count=0
        self.c_bwa.bwa_with_contig(contig, related_refs, sc_seqs, su_seqs, bp_in_contig, is_contig_left, sc_count, su_count)
        return sc_count, su_count
    
    def run_get_contig_supt(self, string contig, vector[string] related_refs, const map[string, vector[string]] q2seqs,
                            int bp1_in_contig, int bp2_in_contig, int flank=5, int min_map2ctg_qual=20):
        cdef vector[string] supt_qnames
        self.c_bwa.run_get_contig_supt(contig, related_refs, q2seqs, supt_qnames, bp1_in_contig, bp2_in_contig, flank, min_map2ctg_qual)
        return supt_qnames

    def __dealloc__(self):
        del self.c_bwa

cdef class LocalSGAWrapper:
    cdef LocalSGA *c_sga
    
    def __cinit__(self, char *id, float error_rate, int min_overlap, int readlen):
        self.c_sga = new LocalSGA(id, error_rate, min_overlap, readlen)
    
    def run(self, vector[string] seqs, int num_assembly_rounds):
        cdef vector[string*] ret_seqs
        self.c_sga.run_asm(seqs, num_assembly_rounds, ret_seqs)
        return [deref(x) for x in ret_seqs]
        
    def __dealloc__(self):
        del self.c_sga