from __future__ import division
import __builtin__
import argparse
from collections import Counter, defaultdict, namedtuple
import copy
import logging
import os
import re
import subprocess
import sys
import time
import traceback

from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import numpy as np
import pysam

'''
bamf = '/mnt/NL200/peizhh/Product_Test/standardSamples/output_1012_2/150000316Dnor_160001141D/cancer/5_recal_bam/150000316Dnor_160001141D_cancer_sort_markdup_realign_recal.bam'
ref = '/mnt/NL200/prod/repos/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa' # use hg19 when bam use this reference
bwa = 'bwa'
samtools = 'samtools'
num_threads = 4
'''

class SvConfig(object):
    ''' Arg '''
    FILTER_FLAG = 4+8+256+1024+2048 # unmap, mate unmap, secondary, PCR duplicate, supplementary
    SU_FILTER_FLAG = 256+2048 # used to find single unmap
    ISIZE_MAX = 600 # the length is decided by library construction(like ultrasonic)
                    # ISIZE_MAX is somewhat high, because of the two main peak(cut of 2 nucleosome) for blood sample
    READ_OVERLAP = 200
    MIN_READ_INNER_SPAN = 50
    MIN_MERGE_SC_LEN = 20
    FIND_FLANK = 10 # allow some bp flank for finding, because the fusion regions close to breakpoints are similar
    MIN_SR_REPORT_COUNT = 2
    MIN_DP_REPORT_COUNT = 5
    ''' MSA Arg '''
    MUSCLE_EXE = os.path.join(os.path.dirname(__file__), 'muscle3.8.31_i86linux64')
    # use maxiters=2 to faster
    MUSCLE_CLINE = MuscleCommandline(MUSCLE_EXE, clwstrict=True, diags=True, tree2='-', maxiters=2)
    # use the small threashold to avoid merge pollution
    PHY_MERGE_THRESHOLD = 0.01
    MIN_UNMAP_MERGE_COUNT = 5
    MAPQ_THR = 25
    
    AVAILABLE_CHR = map(str, range(1, 23)) + ['x', 'y']

def run_sv(res_file, bamf, tmp_path, ref, bwa, samtools,
           num_threads=4, hotregion=None, bam_statf=None, ctrl_data=None, ctrl_res_file=None,
           sample_name=None, is_debug=False, start_time=None):
    '''
    main function to do sv analysis
    Suppose Input Bamfile is sorted by location(default by our workflow)
    '''
    
    ''' 0. Check args and prepare data '''
    num_threads = max(1, num_threads)
    is_case = True if res_file else False
    if not sample_name:
        tmp = os.path.basename(bamf).split('_')[:2]
        sample_name = '_'.join(tmp)
    if not os.path.isdir(tmp_path):
        os.makedirs(tmp_path)
    
    sc_tree_path = []
    tmpfile = os.path.join(tmp_path, 'tmp.file'+('_case' if is_case else '_ctrl'))
    sc_dataf = os.path.join(tmp_path, 'sc.data'+('_case' if is_case else '_ctrl'))
    sc_fqf = os.path.join(tmp_path, 'sc.fq'+('_case' if is_case else '_ctrl'))
    sc_bamf = os.path.join(tmp_path, 'sc.bam'+('_case' if is_case else '_ctrl'))
    pe_bamf = os.path.join(tmp_path, 'pe.bam'+('_case' if is_case else '_ctrl'))
    pe_sort_bamf = os.path.join(tmp_path, 'pe.sort.bam'+('_case' if is_case else '_ctrl'))
    unmap_asm_faf = os.path.join(tmp_path, 'unmap_asm.fa'+('_case' if is_case else '_ctrl'))
    unmap_bamf = os.path.join(tmp_path, 'unmap.bam'+('_case' if is_case else '_ctrl'))
    unmap_sort_bamf = os.path.join(tmp_path, 'unmap.sort.bam'+('_case' if is_case else '_ctrl'))
    unmap_asm_bamf = os.path.join(tmp_path, 'unmap_asm.bam'+('_case' if is_case else '_ctrl'))
    
    bam = pysam.AlignmentFile(bamf)
    header = bam.header
    
    cg = 0
    # refer read length
    for i, r in enumerate(bam):
        if i == 0:
            read_len = len(r.query_sequence)
        if r.flag & SvConfig.FILTER_FLAG:
            continue
        for (cg, cg_len) in r.cigartuples:
            if cg != 0:
                break
        if cg == 0:
            # all match
            read_len = len(r.query_sequence)
            break
    
    # Anno Object
    pos_anno = PosAnno()
    pos_anno.load_bam_header(header)
    if hotregion and os.path.isfile(hotregion):
        pos_anno.load_hotregion(hotregion)
        logging.info('HotRegion data have been loaded from {}'.format(hotregion))
    else:
        logging.warning('No HotRegion file given!')
    
    # open a separate pysam.AlignmentFile object and access each time, this is a iterator
    bam.close()
    
    ''' 1. Stat insert size '''
    isize_mean = None
    if bam_statf:
        # use the exist bam stat file to get isize
        #? this stat file not use the filtered reads to do analysis.
        for line in open(bam_statf):
            # pylint: disable=unused-variable
            itm, data, anno = line.strip().split('\t')
            if itm == 'RD_INS':
                r = re.findall(r'^(\d+)\s*SD\(-(\d+)/\+(\d+)\)$', data)
                isize_mean, isize_lstd, isize_rstd = map(float, r[0]) if r else [None]*3
                break
    if not isize_mean:
        # no bam stat file or not found RD_INS
        isize_mean, isize_lstd, isize_rstd = stat_isize(bamf)
    # pylint: disable=unused-variable
    isize_min = isize_mean - 3*isize_lstd
    isize_max = isize_mean + 3*isize_rstd
    isize_max = SvConfig.ISIZE_MAX if isize_max > SvConfig.ISIZE_MAX else isize_max
    
    logging.info('1. Stat insert size well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())
    
    ''' 2. Abstract reads info '''
    sc_store = ScStore(sc_dataf)
    unmap_store = dict()
    pe_store = dict()
    # use header in higher version pysam, template will segment fault when arg is a closed bam
    pe_bam = pysam.AlignmentFile(pe_bamf, 'wb', header=header)
    unmap_bam = pysam.AlignmentFile(unmap_bamf, 'wb', header=header)
    
    for r in pysam.AlignmentFile(bamf):
        tlen = r.template_length
        
        if r.flag & SvConfig.FILTER_FLAG:
            # Single Unmap Read
            if (not r.flag&SvConfig.SU_FILTER_FLAG) and (r.is_unmapped ^ r.mate_is_unmapped):
                qname = r.query_name
                if qname in unmap_store:
                    # only two reads for one qname, all other mapping is dropped
                    mate_r = unmap_store[qname]
                    del unmap_store[qname]
                    r1, r2 = (r, mate_r) if r.is_unmapped else (mate_r, r)
                    if r2.is_duplicate:
                        #! duplicate is not checked if unmap, so skip mate dup data
                        # this may come from sequencing error
                        continue
                    if validate_chr(r2.reference_name):
                        r1.query_name = '{}={}={}={}'.format(r2.reference_name,
                                                             r2.reference_start+1,
                                                             r2.reference_end,
                                                             r2.mapping_quality)
                        unmap_bam.write(r1)
                else:
                    unmap_store[qname] = r
            continue
        elif r.reference_name != r.next_reference_name or abs(tlen) > isize_max:
            # Non-Discordant Require:
            # 1. map to same chr
            # 2. one positive and another is negative (0x10 ^ 0x20)
            # 3. front one is positive, the another is negative (tlen>0 ^ 0x10)
            # 4. isize is less than isize_mean + 3*isize_std
            #! don't use isize_min to trim, because the reads with large sc will be dropped because of low tlen
            qname = r.query_name
            if qname in pe_store:
                mate_r = pe_store[qname]
                del pe_store[qname]
                if validate_chr(r.reference_name) and validate_chr(mate_r.reference_name):
                    #? use tlen is OK, no need to do this
                    r.query_name = '{}={}={}={}'.format(mate_r.reference_name,
                                                        mate_r.reference_start+1,
                                                        mate_r.reference_end,
                                                        mate_r.mapping_quality)
                    mate_r.query_name = '{}={}={}={}'.format(r.reference_name,
                                                             r.reference_start+1,
                                                             r.reference_end,
                                                             r.mapping_quality)
                    pe_bam.write(r)
                    pe_bam.write(mate_r)
            else:
                pe_store[qname] = r
            continue
        elif r.is_reverse == r.mate_is_reverse or (tlen > 0) == r.is_reverse:
            #TODO, add analysis to same strand, may be plalindrome variants.
            continue
        
        for i, (cg, cg_len) in enumerate(r.cigartuples):
            if cg == 4:
                # For bwa mem, don't use hardclip(cigar=5, flag with 0x800). Because read with hardclip is the mapping of sc part
                if i == 0 and validate_chr(r.reference_name):
                    sc_store.add(r, False, cg_len)
                elif i == len(r.cigartuples)-1 and validate_chr(r.reference_name):
                    sc_store.add(r, True, cg_len)
    
    unmap_bam.close()
    pe_bam.close()
    
    logging.info('2. Abstract reads info well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())
    
    ''' 3. Merge sc '''
    with open(sc_fqf, 'w') as outfile:
        for is_sc_right, sc_name, is_ori_read_reverse, sc_seq_data in sc_store.load_all_by_class():
            for seq, qual, count_list, len_list, ori_map_qual in count_leaves(sc_seq_data, reverse=(not is_sc_right)):
                if len(seq) >= SvConfig.MIN_MERGE_SC_LEN and is_sc_right == is_ori_read_reverse:
                    # merge sc to mapping, which are left part in front read and right part in reverse read
                    sc_count = count_list[-1]
                    #! use the tree_ind not tree_path in query name directly, because too long query name fa is not aligned by bwa
                    tree_ind = len(sc_tree_path)
                    sc_tree_path.append([count_list, len_list])
                    outfile.write('@{}\n{}\n+\n{}\n'.format(
                        '{}={}={}={}={}={}'.format(sc_name, is_sc_right, is_ori_read_reverse, sc_count, ori_map_qual, tree_ind), seq, qual))
                else:
                    #TODO indel sc to infer
                    pass
    
    logging.info('3. Merge sc well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())
    
    ''' 4. Analysis single unmap read '''
    # sort the single unmap bam, it has the same reference_start and reference_name with the mate mapped read
    subprocess.call('{} sort {} > {}'.format(samtools, unmap_bamf, unmap_sort_bamf), shell=True)
            
    with open(unmap_asm_faf, 'w') as outfile:
        # single unmap read has the same is_reverse with the mate mapped read
        for is_reverse, unmap_clustered_reads in cluster_reads_by_fr(pysam.AlignmentFile(unmap_sort_bamf)):
            if len(unmap_clustered_reads) > SvConfig.MIN_UNMAP_MERGE_COUNT:
                unmap_chr = unmap_clustered_reads[0].reference_name
                for merged_seq, pair_start, pair_end, ori_map_qual, unmap_count in do_muscle_merge(unmap_clustered_reads, tmpfile):
                    outfile.write('>{}={}={}={}={}={}\n'.format(unmap_chr, pair_start, pair_end, is_reverse, ori_map_qual, unmap_count))
                    outfile.write('{}\n'.format(merged_seq))
    
    logging.info('4. Analysis single unmap read well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())
    
    ''' 5. Run bwa for sc and signle unmap '''
    # use samse to align, not mem(not work well for short sequence(<70bp))
    subprocess.call('{0} aln -t {1} {2} {3} | {0} samse {2} - {3} | {4} view -Sb - > {5}'.format(
        bwa, num_threads, ref, sc_fqf, samtools, sc_bamf), shell=True)
    subprocess.call('{0} mem -t {1} {2} {3} | {4} view -Sb - > {5}'.format(
        bwa, num_threads, ref, unmap_asm_faf, samtools, unmap_asm_bamf), shell=True)
    
    logging.info('5. Run bwa for sc and signle unmap well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())
    
    ''' 6. get annoed split read result with pe support '''
    # sort the pe bam and index
    subprocess.call('{0} sort {1} > {2} && {0} index {2}'.format(samtools, pe_bamf, pe_sort_bamf), shell=True)
    
    merge_flank = MergeFlank(sample_name)
    ResInfo = namedtuple('ResInfo',
                         '''source, hot_status, fusion_type, chr1, bp_pos1, start1, end1, chr2, bp_pos2, start2, end2,
                            is_gene1R, is_reverse1, is_gene2R, is_reverse2,
                            sr_segcigar, su_segcigar, sc_count, unmap_count, pe_count, template_count, mapq''')
    skip_region = defaultdict(list)
    read_inner_span = isize_max-2*read_len
    read_inner_span = SvConfig.MIN_READ_INNER_SPAN if read_inner_span < SvConfig.MIN_READ_INNER_SPAN else read_inner_span
    
    merge_sr = defaultdict(list)
    for sc_data in get_sc_info(sc_bamf):
        # pylint: disable=unused-variable
        (chr1, pos1, is_reverse1, is_gene1R, chr2, pos2, is_reverse2, is_gene2R, 
         is_sc_right, is_same_strand, sc_cigar_md, sc_seq, sc_mapq, sc_count, tree_ind) = sc_data
       
        #! In SR, the reverse1 and reverse2 are the diferent parts of same read
        #      so should use the strand info of mate read(reverse complement to sc read) as is_reverse1
        is_reverse1 ^= 1
        validate_sr = validate_pos('SR', chr1, pos1, pos1, pos1, chr2, pos2, pos2, pos2,
                                   is_gene1R, is_reverse1, is_gene2R, is_reverse2,
                                   pos_anno, isize_max, sc_count, SvConfig.MIN_SR_REPORT_COUNT,
                                   read_inner_span, pe_sort_bamf, skip_region)
        if not validate_sr:
            continue
        hot_status, fusion_type, pe_count, pe_mapq = validate_sr
        if pe_mapq is not None:
            sc_mapq = pe_mapq # use pe as mapq
        
        template_count = count_template_coverage(bamf, chr1, pos1, isize_max)
        
        sc_info_tmp = [sc_cigar_md, sc_seq, sc_count, # used to merge
                       ['SplitRead', hot_status, fusion_type,
                        chr1, pos1, pos1, pos1, chr2, pos2, pos2, pos2,
                        is_gene1R, is_reverse1, is_gene2R, is_reverse2,
                        sc_cigar_md, '*', sc_count, 0, pe_count, template_count, sc_mapq]]
        
        key = '{}={}={}={}={}={}={}={}'.format(
            chr1, pos1, is_reverse1, is_gene1R, chr2, pos2, is_reverse2, is_gene2R)
        #! seq need reverse to be identical to count_list/len_list
        # for current sc_seq get by re-bwa, only think about whether it's reverse, 
        # no relation with seq in sc.fq which may be complement reverse or equal to sc_seq
        leaf = [[sc_seq if is_reverse2 else sc_seq[::-1]] + sc_tree_path[tree_ind], sc_info_tmp]
        merge_sr[key].append(leaf)
    
    # merge the sc count of similar sc sequence(mapping to same position)
    for key, sc_tree_info in merge_sr.items():
        leaves, sc_info_tmps = map(list, zip(*sc_tree_info))
        merge_sc_count = merge_leaves_count(leaves)
        
        sc_info_tmps.sort(key=lambda _: _[2], reverse=True) # sort by sc_count
        merge_seq = []
        for sc_info_tmp in sc_info_tmps:
            sc_cigar_md, sc_seq = sc_info_tmp[:2]
            merge_seq.append('{}({})'.format(sc_seq, sc_cigar_md))
        
        # pe_count and template_count are same, and use the largest count item.
        merge_seq_cigar = '_'.join(merge_seq)
        # remove the sc_cigar_md, sc_seq, sc_count which are used to merge
        sc_res_info = sc_info_tmps[0][3]
        sc_res_info[-7] = merge_seq_cigar
        sc_res_info[-5] = merge_sc_count
        sc_res_info = ResInfo(*sc_res_info)
        
        merge_flank.add(sc_res_info, is_case=is_case)
    
    logging.info('6. get annoed split read result with pe support well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())
    
    ''' 7. get annoed single unmap result '''
    for unmap_data in get_unmap_info(unmap_asm_bamf):
        (chr1, unmap_pos1, start1, end1, is_reverse1, is_gene1R,
         chr2, unmap_pos2, start2, end2, is_reverse2, is_gene2R,
         unmap_cigar, unmap_seq, unmap_mapq, unmap_count) = unmap_data

        validate_su = validate_pos('SU', chr1, unmap_pos1, start1, end1, chr2, unmap_pos2, start2, end2,
                                   is_gene1R, is_reverse1, is_gene2R, is_reverse2,
                                   pos_anno, isize_max, unmap_count, SvConfig.MIN_SR_REPORT_COUNT,
                                   read_inner_span, pe_sort_bamf, skip_region)
        if not validate_su:
            continue
        # pylint: disable=unused-variable
        hot_status, fusion_type, pe_count, pe_mapq = validate_su
        if pe_mapq is not None:
            unmap_mapq = pe_mapq # use pe as mapq
        
        template_count = count_template_coverage(bamf, chr1, unmap_pos1, isize_max)
        
        su_res_info = ResInfo('SingleUnmap', hot_status, fusion_type,
                              chr1, unmap_pos1, start1, end1, chr2, unmap_pos2, start2, end2,
                              is_gene1R, is_reverse1, is_gene2R, is_reverse2,
                              '*', '{}({})'.format(unmap_seq, unmap_cigar),
                              0, unmap_count, pe_count, template_count, unmap_mapq)
        merge_flank.add(su_res_info, is_case=is_case)
    
    #TODO, merge the skip region which is used by sc and unmap
    
    logging.info('7. get annoed single unmap result well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())
    
    ''' 8. get annoed discordant pair result '''
    for pe_merged_data in get_pe_info(bamf, isize_max, pe_sort_bamf, pos_anno, skip_region=skip_region):
        (chr1, pe_pos1, start1, end1, is_reverse1, is_gene1R,
         chr2, pe_pos2, start2, end2, is_reverse2, is_gene2R,
         pe_mapq, pe_count) = pe_merged_data
       
        validate_dp = validate_pos('DP', chr1, pe_pos1, start1, end1, chr2, pe_pos2, start2, end2,
                                   is_gene1R, is_reverse1, is_gene2R, is_reverse2,
                                   pos_anno, isize_max, pe_count, SvConfig.MIN_DP_REPORT_COUNT)
        if not validate_dp:
            continue
        hot_status, fusion_type = validate_dp
        
        template_count = count_template_coverage(bamf, chr1, pe_pos1, isize_max)
        if hot_status.is_hotfusion:
            # do unmap and sc finding for only hot fusion
            #! cost long time, so not for all pe to do. 
            sc_count, sc_in_pe = get_sc_within_pe(chr1, pe_pos1, is_reverse1, sc_bamf)
            unmap_count, unmap_in_pe = get_unmap_within_pe(chr1, pe_pos1, is_reverse1, unmap_asm_bamf)
        else:
            sc_count, sc_in_pe, unmap_count, unmap_in_pe = 0, '*', 0, '*'
        
        dp_res_info = ResInfo('DiscordantPair', hot_status, fusion_type,
                              chr1, pe_pos1, start1, end1, chr2, pe_pos2, start2, end2,
                              is_gene1R, is_reverse1, is_gene2R, is_reverse2,
                              sc_in_pe, unmap_in_pe, sc_count, unmap_count, pe_count, template_count, pe_mapq)
        merge_flank.add(dp_res_info, is_case=is_case)
    
    logging.info('8. get annoed discordant pair result well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())
    
    ''' 9. write result '''
    if is_case:
        if ctrl_data:
            merge_flank.load_ctrl(ctrl_data)
        merge_flank.merge_bp()
        merge_flank.write(res_file, RES_TITLE)
        if ctrl_res_file:
            merge_flank.write(ctrl_res_file, RES_TITLE, only_ctrl=True)
    else:
        return merge_flank
    logging.info('9. write result well-done' + ', Passing {} times 30s'.format(int((time.time()-start_time)/30)) if start_time else '')
    if is_debug:
        logging.warning(hpy().heap())

'''################ Functions and Class ################'''
# -1. same name function from .utils, make this file stand alone
def rm_chr_prefix(data):
    return '\n'.join([s[3:] if s[:3].lower() == 'chr' else s for s in data.split('\n')])

def validate_chr(chr):
    '''
    validate whether the chromsome is sv candidate
    '''
    return True if rm_chr_prefix(chr).lower() in SvConfig.AVAILABLE_CHR else False

# 1.
def stat_isize(bamf):
    '''
    count the mode and rms(root mean square) of bam insert size
    '''
    isize = defaultdict(int)
    for r in pysam.AlignmentFile(bamf):
        #! sc at the both ends will make isize less than the real, but only right std used, so it's ok.
        tlen = r.template_length
        if r.flag & SvConfig.FILTER_FLAG:
            # remove supplementary(in FILTER_FLAG) flag for counting size, isize of 0x800 is 0.
            # bam is not got by bwa mem -a, no need to 256(0x100)
            continue
        if (r.reference_name != r.next_reference_name
                or r.is_reverse == r.mate_is_reverse
                or (tlen > 0) == r.is_reverse):
            # continue if not normal read pair
            continue
        if r.is_read2:
            # only count one read in pe
            continue
        isize[abs(tlen)] += 1
    isize_mode = sorted(isize.items(), key=lambda _: _[1], reverse=True)[0][0]
    # infer the std
    lall = sum([c for v, c in isize.items() if v < isize_mode])
    rall = sum([c for v, c in isize.items() if v > isize_mode])
    mdcnt = isize[isize_mode]
    lcnt, lrms, rcnt, rrms = 0, 0, 0, 0
    for v, c in sorted(isize.items(), key=lambda _: _[0] - isize_mode):
        # from center to outer to filter the abnormal value
        if v < isize_mode and lcnt < lall*0.99:
            lcnt += c
            lrms += (v-isize_mode)**2 * c
        elif v > isize_mode and rcnt < rall*0.99:
            rcnt += c
            rrms += (v-isize_mode)**2 * c
    # try to emulate normal distribution.
    lrms = np.sqrt(lrms*2/(lcnt*2+mdcnt))
    rrms = np.sqrt(rrms*2/(rcnt*2+mdcnt))
    return isize_mode, lrms, rrms

# 3.
def longest_common_string(a, b):
    ''' Deprcated '''
    s = np.zeros([len(a)+1, len(b)+1])
    for ai, a0 in enumerate(a):
        for bi, b0 in enumerate(b):
            if a0 == b0:
                s[ai+1][bi+1] = s[ai][bi] + 1
            else:
                s[ai+1][bi+1] = max(s[ai+1][bi], s[ai][bi+1])
    return s[-1][-1]

def count_leaves(seq_data, reverse=False):
    ''' 
    Count leaves from a Tree(emulated).
    Tree Feat.: 
        1. Each node is a sequence
        2. A node is the children of another node if and only if its's sequence is another node's prefix
        3. A node's count is the occurrence of this node and its parent's count(if exist)
    '''
    if seq_data:
        seq_dict = defaultdict(int)
        qual_dict = defaultdict(str)
        ori_map_qual_list = []
        for seq, qual, ori_map_qual in seq_data:
            if reverse:
                seq = seq[::-1]
            seq_dict[seq] += 1
            #TODO use the average qual?
            qual_dict[seq] = qual
            ori_map_qual_list.append(ori_map_qual)
        seq_list, count_list, len_list = [], [], []
        ori_map_qual = sorted(ori_map_qual_list, reverse=True)[int(3/4*len(seq_data))] # use the Q3 as qual
        # sort by alphabet order to emulate the tree, and preorder traversal
        for seq, count in sorted(seq_dict.items(), key=lambda _: _[0]):
            called = False
            while seq_list:
                s = seq_list[-1] 
                if len(s) >= len(seq) or s != seq[:len(s)]:
                    if not called:
                        # use copy.deepcopy to force return a new object, not refer to the same one
                        yield (s[::-1] if reverse else s), qual_dict[s], copy.deepcopy(count_list), copy.deepcopy(len_list), ori_map_qual
                        called = True
                    seq_list.pop()
                    count_list.pop()
                    len_list.pop()
                else:
                    seq_list.append(seq)
                    # accumulated count
                    count_list.append(count+count_list[-1])
                    len_list.append(len(seq))
                    break
            if not seq_list:
                # exit to the outermoust
                seq_list.append(seq)
                count_list.append(count)
                len_list.append(len(seq))
        # output the last one
        s = seq_list[-1]
        yield (s[::-1] if reverse else s), qual_dict[s], copy.deepcopy(count_list), copy.deepcopy(len_list), ori_map_qual

# 4.
def cluster_reads_by_fr(bam_or_reads, read_overlap=SvConfig.READ_OVERLAP, use_mate=False):
    ''' 
    use sorted bam/reads as input
        to cluster reads whose distance of neighbour is within read_overlap
    '''
    pre_chr, pre_pos = [-1, -1], [-1, -1]
    clustered_reads = [[], []]
    for r in bam_or_reads:
        chr = r.next_reference_name if use_mate else r.reference_name
        start = r.next_reference_start if use_mate else r.reference_start
        rev = (r.mate_is_reverse if use_mate else r.is_reverse)&1
        if chr != pre_chr[rev] or start-pre_pos[rev] > read_overlap:
            if clustered_reads[rev]:
                yield rev, clustered_reads[rev]
            pre_chr[rev] = chr
            clustered_reads[rev] = [r]
        else:
            clustered_reads[rev].append(r)
        pre_pos[rev] = start
    # last one
    for rev, reads in enumerate(clustered_reads):
        if reads:
            # read is Empty if input is empty
            yield rev, reads

def do_muscle_merge(unmap_reads, tmpfile):
    ''' 
    return the merged sequence by clustered unmap_reads and left, right boundary of used reaads start position 
    '''
    child = subprocess.Popen(str(SvConfig.MUSCLE_CLINE), stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             shell=(sys.platform != "win32"))
    input_data, phy = [], []
    for r in unmap_reads:
        seq = r.query_sequence
        # use the unmap data, trim the N (use N to mask Adapter and make up read len)
        while seq.upper().startswith('N'):
            seq = seq[1:]
        while seq.upper().endswith('N'):
            seq = seq[:-1]
        if seq:
            input_data.append('>{}\n{}\n'.format(r.query_name, seq))
    with open(tmpfile, 'w') as outfile:
        write_aln = False
        for line in child.communicate(''.join(input_data))[0].split('\n'):
            if write_aln:
                outfile.write('{}\n'.format(line))
            else:
                phy.append(line)
                if line == ';':
                    # phylip is over
                    write_aln = True
    try:
        align = AlignIO.read(tmpfile, "clustal")
    except ValueError:
        # No records found in handle
        align = [[]]
        phy = [';']
        logging.warning('one unmap record is error in IO')
    seq_len = len(align[0])
    for seqid in parse_phy(phy):
        merged_seq = []
        for i in range(seq_len):
            c = [align[j, i] for j in seqid]
            a = Counter([base for base in c if base != '-'])
            base = a.most_common(1)
            if base:
                # not all gap
                merged_seq.append(base[0][0])
        ori_map_qual_list, starts, ends = [], [], []
        for qname in [align[i].id for i in seqid]:
            # id is the header of fasta
            # pylint: disable=unused-variable
            chr, start, end, ori_map_qual = qname.split('=')
            starts.append(int(start))
            ends.append(int(end))
            ori_map_qual_list.append(int(ori_map_qual))
        
        ori_map_qual = sorted(ori_map_qual_list, reverse=True)[int(3/4*len(seqid))] # use the Q3 as qual
        yield ''.join(merged_seq), min(starts), max(ends), ori_map_qual, len(seqid)

def parse_phy(phy):
    '''
    parse the phylip result
    '''
    #TODO count the accumulated value for compare
    def compare_phy(pre, cur):
        av, groupa = pre
        bv, groupb = cur
        global merge_group
        breaka, breakb = False, False
        merge_value = av < SvConfig.PHY_MERGE_THRESHOLD + (bv < SvConfig.PHY_MERGE_THRESHOLD)*2
        if merge_value&3:
            # could be merge
            merge_group = groupa+groupb
        elif merge_value&2:
            merge_group = groupb
            breaka = True
        elif merge_value&1:
            merge_group = groupa
            breakb = True
        else:
            breaka, breakb = True, True
            merge_group = []
        if breaka and len(groupa) >= SvConfig.MIN_UNMAP_MERGE_COUNT:
            return groupa
        if breakb and len(groupb) >= SvConfig.MIN_UNMAP_MERGE_COUNT:
            return groupb
        
    # remove the stop mark ;
    phy = phy[:-1]
    current_leaf, current_depth = -1, 0
    global merge_group
    upper, pair_depth, merge_group = [], [0], []
    for p in phy:
        p = p.strip()
        if p == '(':
            current_depth += 1
        elif p.startswith(')'):
            current_depth -= 1
            match = re.findall(r'\:(.*)', p)
            if match:
                # not the last )
                value = float(match[0])
                if pair_depth[-1] == current_depth:
                    pair_depth.pop()
                    r = compare_phy(upper.pop(), (value, merge_group))
                    if r:
                        yield r
                else:
                    upper.append((value, merge_group))
        elif p == ',':
            pair_depth.append(current_depth)
        else:
            current_leaf += 1
            value = float(re.findall(r'\:(.*)', p)[0])
            if pair_depth[-1] == current_depth:
                pair_depth.pop()
                r = compare_phy(upper.pop(), (value, [current_leaf]))
                if r:
                    yield r
            else:
                upper.append((value, [current_leaf]))

# 6.
def validate_pos(type, chr1, pos1, start1, end1, chr2, pos2, start2, end2,
                 is_gene1R, is_reverse1, is_gene2R, is_reverse2,
                 pos_anno, isize_max, count, count_limit,
                 read_inner_span=None, pe_sort_bamf=None, skip_region=None):
    if not validate_chr(chr2):
        # chr1 has been validated for SR,SU.
        # both chr1,2 have been validated for DP
        return []
    
    #TODO, fix specific SR
    if type != 'SR' and count < count_limit:
        # judge SR count after pe count
        return []
    elif chr1 == chr2 and (start1-isize_max <= start2 <= end1+isize_max
                           or start1-isize_max <= end2 <= end1+isize_max):
        # for SU/DP, two side overlap
        #TODO, for SR
        #    too close, indel and repeat is not considered now
        #    it also can't find pe support, because min distance of discordant pe is isize_max
        return []
    
    # anno to hot region
    fusion_type = FusionType.parse_from_strand(chr1, pos1, is_reverse1, chr2, pos2, is_reverse2)
    hot_status = pos_anno.anno2hotregion(chr1, pos1, chr2, pos2, is_gene1R, is_gene2R, fusion_type)
    
    if type == 'DP':
        # DP not find pe support
        return hot_status, fusion_type
    
    # find pe support, pe_count and pe_same_strand_count
    pe_count, pe_mapq, block_start1, block_end1, block_start2, block_end2 = get_pe_support(
        chr1, start1, end1, is_reverse1, is_gene1R,
        chr2, start2, end2, is_reverse2, is_gene2R,
        read_inner_span, pe_sort_bamf, pos_anno)
    if pe_count:
        skip_region[chr1].append([block_start1, block_end1, is_reverse1, chr2, block_start2, block_end2, is_reverse2])
        skip_region[chr2].append([block_start2, block_end2, is_reverse2, chr1, block_start1, block_end1, is_reverse1])
    elif type == 'SR' and count < count_limit:
        # only filter the no pe support result
        return []
    return hot_status, fusion_type, pe_count, pe_mapq

def get_sc_info(sc_bamf):
    '''
    got the sc(SR) result from abstracted bam
    '''
    # is_sc_right represents the relative location of each other
    # gene1/2LR represents the relative location of breakpoint
    for r in pysam.AlignmentFile(sc_bamf):
        if r.flag & SvConfig.FILTER_FLAG:
            continue
        
        chr1, pos1, is_sc_right, is_reverse1, sc_count, ori_map_qual, tree_ind = r.query_name.split('=')
        chr2, is_reverse2 = r.reference_name, r.is_reverse
        pos1, is_reverse1, is_reverse2, is_sc_right, sc_count, ori_map_qual, tree_ind = map(
            int, [pos1, is_reverse1, is_reverse2, is_sc_right, sc_count, ori_map_qual, tree_ind])
        #! the sc part has been reverse complement if original is negative, so need to reverse to back to real strand.
        if is_reverse1:
            is_reverse2 = is_reverse2 ^ 1
        is_same_strand = is_reverse1 == is_reverse2
        is_gene1R = 0 if is_sc_right else 1
        if is_same_strand ^ is_sc_right:
            pos2 = r.reference_end
            is_gene2R = 0
        else:
            pos2 = r.reference_start+1
            is_gene2R = 1
        
        sc_mapq = min(r.mapq, ori_map_qual) # use the min one as mapping quality
        #! MD is only for aln?
        yield (chr1, pos1, is_reverse1, is_gene1R, chr2, pos2, is_reverse2, is_gene2R,
               is_sc_right, is_same_strand, '{}/{}'.format(r.cigarstring, r.get_tag('MD')),
               r.query_sequence, sc_mapq, sc_count, tree_ind)

def get_pe_support(chr1, start1, end1, is_reverse1, is_gene1R,
                   chr2, start2, end2, is_reverse2, is_gene2R,
                   read_inner_span, pe_sort_bamf, pos_anno):
    '''
    use the indexed pe_bam to find suppport pe of the given fusion region
    '''
    pe_count = 0
    pe_mapq_all = []
    pe_support_flank = SvConfig.FIND_FLANK
    if is_gene1R:
        block_start1, block_end1 = max(1, start1-pe_support_flank), min(start1+read_inner_span-1, pos_anno.chrlen[chr1])
    else: 
        block_start1, block_end1 = max(1, end1-read_inner_span+1), min(end1+pe_support_flank, pos_anno.chrlen[chr1])
    if is_gene2R:
        block_start2, block_end2 = max(1, start2-pe_support_flank), min(start2+read_inner_span-1, pos_anno.chrlen[chr2])
    else:
        block_start2, block_end2 = max(1, end2-read_inner_span+1), min(end2+pe_support_flank, pos_anno.chrlen[chr2])
    
    for pe_read in pysam.AlignmentFile(pe_sort_bamf).fetch(chr1, block_start1-1, block_end1):
        # pylint: disable=unused-variable
        mate_chr, mate_start, mate_end, mate_mapq = pe_read.query_name.split('=')
        mate_start, mate_end = map(int, [mate_start, mate_end])
        
        if chr2 != mate_chr:
            continue
        if is_gene1R and not block_start1 <= pe_read.reference_start+1 <= block_end1:
            continue
        elif not is_gene1R and not block_start1 <= pe_read.reference_end <= block_end1:
            continue
        elif is_gene2R and not block_start2 <= mate_start <= block_end2:
            continue
        elif not is_gene2R and not block_start2 <= mate_end <= block_end2:
            continue
        
        if pe_read.is_reverse == is_reverse1 and pe_read.mate_is_reverse == is_reverse2:
            pe_mapq_all.append(min(int(mate_mapq), pe_read.mapping_quality))  # use the min part mapq of pe
            pe_count += 1
    if pe_count:
        pe_mapq_all.sort(reverse=True)
        pe_mapq = pe_mapq_all[int(3/4*pe_count)] # use the Q3 as mapq of this dp
    else:
        pe_mapq = None
    return pe_count, pe_mapq, block_start1, block_end1, block_start2, block_end2

def count_template_coverage(bamf, chr, pos, isize_max):
    '''
    Count the `Non-Abnormal Template` coverage of 1-based position 
    Split Reads will also be counted.
    '''
    count = 0
    isize_max = 1 if isize_max < 1 else isize_max
    for r in pysam.AlignmentFile(bamf).fetch(chr, max(0, pos-isize_max), pos):
        if r.flag&SvConfig.FILTER_FLAG:
            continue
        tlen = r.template_length
        if tlen <= 0 or tlen > isize_max:
            # one part of pe is unmap -> tlen=0
            # pe is not at same chromsome -> tlen=0
            # not allow too long tlen, max is isize_max
            continue
        # no need r.is_reverse == r.mate_is_reverse or (tlen>0) == r.is_reverse, allow these reads
        if r.reference_start+1 <= pos <= r.reference_start+tlen:
            count += 1
    return count

def merge_leaves_count(leaves):
    '''
    Merge leaves by count the common node once
    leaves:
        (seq, count_list, len_list) for one leaf
    '''
    if not leaves:
        return 0
    leaves.sort(key=lambda _: _[0])
    seq_prev, count_list, len_list_prev = leaves[0]
    total_count = count_list[-1]
    for seq, count_list, len_list in leaves[1:]:
        max_ind = -1
        for ind, (l1, l2) in enumerate(zip(len_list_prev, len_list)):
            if seq[:l1] == seq_prev[:l2]:
                max_ind = ind
        if max_ind >= 0:
            # has common prefix
            total_count += count_list[-1] - count_list[max_ind]
        else:
            total_count += count_list[-1]
        seq_prev, len_list_prev = seq, len_list
    return total_count

# 7.
def get_unmap_info(unmap_asm_bamf):
    '''
    got the single unmap(SU) result from abstracted bam
    '''
    for r in pysam.AlignmentFile(unmap_asm_bamf):
        if r.is_unmapped:
            continue
        chr1, start1, end1, is_reverse1, ori_map_qual, unmap_count = r.query_name.strip().split('=')
        start1, end1, is_reverse1, ori_map_qual, unmap_count = map(int, [start1, end1, is_reverse1, ori_map_qual, unmap_count])
        chr2, start2, end2 = r.reference_name, r.reference_start+1, r.reference_end
        
        # sequence of unmap is complement reversed if mate is_reverse because of the same is_reverse flag
        is_reverse2 = r.is_reverse&1
        if is_reverse1:
            is_reverse2 ^= 1
        
        # gene1 is mapped mate reads, and gene2 is unmap merged reads
        is_same_strand = is_reverse1 == is_reverse2
        is_gene1R = 1 if is_reverse1 else 0
        is_gene2R = is_gene1R if is_same_strand else is_gene1R^1
        
        unmap_pos1 = start1 if is_gene1R else end1
        unmap_pos2 = start2 if is_gene2R else end2
        
        unmap_mapq = min(r.mapq, ori_map_qual)
        yield (chr1, unmap_pos1, start1, end1, is_reverse1, is_gene1R,
               chr2, unmap_pos2, start2, end2, is_reverse2, is_gene2R,
               r.cigarstring, r.query_sequence, unmap_mapq, unmap_count)

# 8.
def get_pe_cluster_by_dist(pe_reads, rev, rev_mate, read_overlap=SvConfig.READ_OVERLAP):
    is_same_strand = rev == rev_mate
    is_R = 1 if rev else 0
    is_mate_R = is_R if is_same_strand else is_R^1
    
    PeAttr = namedtuple('PeAttr', 'chr1, start1, end1, chr2, start2, end2, mapq')
    # strands can decide bp location
    pes = []
    for r in pe_reads:
        chr, start, end, mapq = r.query_name.split('=')
        pes.append(PeAttr(r.reference_name,
                          r.reference_start,
                          r.reference_end,
                          chr,
                          int(start),
                          int(end),
                          min(int(mapq), r.mapping_quality)))  # use the min part mapq of pe
    
    if is_R:
        bps = [p.start1 for p in pes]
        bp = min(bps)
    else:
        bps = [p.end1 for p in pes]
        bp = max(bps)
    if is_mate_R:
        bp_mates = [p.start2 for p in pes]
        bp_mate = min(bp_mates)
    else:
        bp_mates = [p.end2 for p in pes]
        bp_mate = max(bp_mates)
    
    dists = []
    for i, (current_bp, current_mate_bp) in enumerate(zip(bps, bp_mates)):
        dists.append([i, abs(current_bp-bp) + abs(current_mate_bp-bp_mate)])
    dists.sort(key=lambda _: _[1]) # sort by dist
    if dists[-1][1] - dists[0][1] <= read_overlap:
        # no need to separate
        yield pes
    else:
        pre_dist = dists[0][1] # first one
        pre_group = []
        for i, d in dists:
            if abs(d-pre_dist) > read_overlap:
                yield pre_group
                pre_group = [pes[i]]
            else:
                pre_group.append(pes[i])
            pre_dist = d
        yield pre_group
    
def get_pe_info(bamf, isize_max, pe_sort_bamf, pos_anno, skip_region=None):
    '''
    get pe(DP) result from abnormal pe which are not used to support SR/SU
    '''
    for is_reverse1, pe_clustered_reads in cluster_reads_by_fr(pysam.AlignmentFile(pe_sort_bamf)):
        pe_mate_sorted_reads = sorted(pe_clustered_reads, cmp=lambda a, b: pos_anno.cmp_reads(a, b, use_mate=True))
        for is_reverse2, pe_mate_clustered_reads in cluster_reads_by_fr(pe_mate_sorted_reads, use_mate=True):
            if len(pe_mate_clustered_reads) < SvConfig.MIN_DP_REPORT_COUNT:
                # quick pass
                continue
            for pes in get_pe_cluster_by_dist(pe_mate_clustered_reads, is_reverse1, is_reverse2):
                if len(pes) < SvConfig.MIN_DP_REPORT_COUNT:
                    continue
                
                chr1 = pes[0].chr1
                start1 = min([p.start1 for p in pes]) 
                end1 = max([p.end1 for p in pes])
                chr2 = pes[0].chr2
                start2 = min([p.start2 for p in pes]) 
                end2 = max([p.end2 for p in pes])
                pe_mapq_all = [p.mapq for p in pes] 
                
                if pos_anno.cmp_chrpos((chr1, start1), (chr2, start2)) >= 0:
                    # report one side of pe
                    continue
                
                # skip the sc and unmap report region
                if skip_region:
                    is_skip = False
                    for block_start1, block_end1, block_is_reverse1, block_chr2, block_start2, block_end2, block_is_reverse2 in skip_region[chr1]:
                        if chr2 != block_chr2:
                            continue
                        if ((block_start1 <= start1 <= block_end1 or block_start1 <= end1 <= block_end1)
                                and (block_start2 <= start2 <= block_end2 or block_start2 <= end2 <= block_end2)
                                and (is_reverse1 == block_is_reverse1)
                                and (is_reverse2 == block_is_reverse2)):
                            is_skip = True
                            break
                    if is_skip:
                        continue
                
                #! save the is_reverse1, not change it because it used in a for loop
                cur_rev1, cur_rev2 = is_reverse1, is_reverse2
                
                is_same_strand = cur_rev1 == cur_rev2
                is_gene1R = 1 if cur_rev1 else 0
                is_gene2R = is_gene1R if is_same_strand else is_gene1R^1
                
                pe_pos1 = start1 if is_gene1R else end1
                pe_pos2 = start2 if is_gene2R else end2
                
                # use the more abundant part as main
                sum1 = count_template_coverage(bamf, chr1, pe_pos1, isize_max)
                sum2 = count_template_coverage(bamf, chr2, pe_pos2, isize_max)
                
                if sum2 > sum1:
                    chr1, pe_pos1, start1, end1, cur_rev1, is_gene1R, chr2, pe_pos2, start2, end2, cur_rev2, is_gene2R = (
                        chr2, pe_pos2, start2, end2, cur_rev2, is_gene2R, chr1, pe_pos1, start1, end1, cur_rev1, is_gene1R)
                
                pe_count = len(pes)
                pe_mapq_all.sort(reverse=True)
                pe_mapq = pe_mapq_all[int(3/4*pe_count)] # use the Q3 as mapq of this dp
                yield (chr1, pe_pos1, start1, end1, cur_rev1, is_gene1R,
                       chr2, pe_pos2, start2, end2, cur_rev2, is_gene2R,
                       pe_mapq, pe_count)

def get_sc_within_pe(chr1, pe_pos, is_reverse1, sc_bamf):
    # back trace the sc sequence for pe, not match or match to other place
    flank = SvConfig.FIND_FLANK
    records = defaultdict(int)
    max_count, max_pos, max_info = 0, 0, ''
    for r in pysam.AlignmentFile(sc_bamf):
        if not r.qname.startswith('{}='.format(chr1)):
            # judge chr1==chr, quick pass
            continue
        # pylint: disable=unused-variable
        chr, pos, is_sc_right, is_ori_read_reverse, sc_count, ori_map_qual, _ = r.qname.strip().split('=')
        pos, is_ori_read_reverse, sc_count = map(int, [pos, is_ori_read_reverse, sc_count])
        if pe_pos-flank <= pos <= pe_pos+flank and is_ori_read_reverse^is_reverse1:
            # record the max count of this pos
            if sc_count > max_count:
                max_count, max_pos = sc_count, pos
                max_info = '{}({}={}={})'.format(
                    r.seq, 
                    '*' if r.is_unmapped else r.cigarstring, 
                    '{}:{}'.format(chr, pos),
                    'unmap' if r.is_unmapped else '{}:{}'.format(r.reference_name, r.reference_start+1))
                records[pos] += 1
    if records:
        # sc use prefix tree, so all other results count as 1
        sc_count, sc_in_pe = max_count+records[max_pos]-1, max_info
    else:
        sc_count, sc_in_pe = 0, '*'
    return sc_count, sc_in_pe

def get_unmap_within_pe(chr1, pe_pos, is_reverse1, unmap_asm_bamf):
    # back trace the unmap sequence for pe, not match or match to other place
    flank = SvConfig.FIND_FLANK
    records = []
    record_count = 0
    use_unmap_r = 0 if is_reverse1 else 1 # plus strand use end(start + end .. mate)
    
    for r in pysam.AlignmentFile(unmap_asm_bamf):
        if not r.qname.startswith('{}='.format(chr1)):
            # judge chr1==chr, quick pass
            continue
        chr, start, end, is_reverse, ori_map_qual, unmap_count = r.qname.strip().split('=')
        pos = end if use_unmap_r else start
        pos, is_reverse, unmap_count = map(int, [pos, is_reverse, unmap_count])
        if pe_pos-flank <= pos <= pe_pos+flank and is_reverse == is_reverse1:
            # record the max count of this pos
            record_count += unmap_count
            records.append('{}({}={}={})'.format(
                r.seq,
                '*' if r.is_unmapped else r.cigarstring, 
                '{}:{}'.format(chr, pos),
                'unmap' if r.is_unmapped else '{}:{}'.format(r.reference_name, r.reference_start+1)))
    if records:
        # unmap use all counts, they are not use common prefix
        unmap_count, unmap_in_pe = record_count, '='.join(records)
    else:
        unmap_count, unmap_in_pe = 0, '*'
    return unmap_count, unmap_in_pe

class ScStore(object):
    '''
    class to deal with sc abstracted from original bam
    '''
    def __init__(self, sc_dataf):
        ''' 
        data(list):
            item 0(dict): data for sc_left
                key: sc_name
                    {chr}_{break_point}
                value(dict):
                    key: is_original_read_negative
                        0/1
                    value(list): sequence info
                        [sc_sequence, sc_quality]
            item 1: data for sc_right, structure is same to item 0
        '''
        self.dataf = sc_dataf
        self.outfile = open(sc_dataf, 'w')
    
    def add(self, r, is_sc_right, cg_len):
        '''
        Add the read r to self data
        '''
        query_seq = r.query_sequence
        query_qual = ''.join([__builtin__.chr(_+33) for _ in r.query_qualities])
        sc_seq = query_seq[-cg_len:] if is_sc_right else query_seq[:cg_len]
        sc_qual = query_qual[-cg_len:] if is_sc_right else query_qual[:cg_len]
        # trim N(only in end, generated by adapter mask)
        while sc_seq.upper().startswith('N'):
            sc_seq = sc_seq[1:]
            sc_qual = sc_qual[1:]
        while sc_seq.upper().endswith('N'):
            sc_seq = sc_seq[:-1]
            sc_qual = sc_qual[:-1]
        if not sc_seq:
            # SC part is all N
            return
        ori_map_qual = r.mapping_quality
        sc_name = '{}={}'.format(r.reference_name, 
                                 # trans position to 1-base in reality
                                 # reference_end is the one past the last aligned residue
                                 r.reference_end if is_sc_right else r.reference_start+1)
        # is_sc_right used as index is auto convert to 0/1
        # r.is_reverse as key is Bool, so use bit op & to convert to 0/1
        #TODO ?mate start(1-based) 
        
        self.outfile.write('{}:{}:{}\t{}\t{}\t{}\n'.format(is_sc_right&1, sc_name, r.is_reverse&1, sc_seq, sc_qual, ori_map_qual))
    
    def _flush(self):
        self.outfile.flush()
    
    def load_all_by_class(self):
        self._flush()
        self._sort()
        
        pre_name = ''
        pre_list = []
        for line in open(self.sort_dataf):
            tmp = line.strip().split('\t')
            if tmp[0] != pre_name:
                if pre_list:
                    is_sc_right, sc_name, is_ori_read_reverse = pre_name.split(':')
                    yield int(is_sc_right), sc_name, int(is_ori_read_reverse), pre_list
                pre_list = []
                pre_name = tmp[0]
            sc_seq, sc_qual, ori_map_qual = tmp[1:]
            pre_list.append([sc_seq, sc_qual, int(ori_map_qual)])
        # last one
        if pre_list:
            is_sc_right, sc_name, is_ori_read_reverse = pre_name.split(':')
            yield int(is_sc_right), sc_name, int(is_ori_read_reverse), pre_list
    
    def _sort(self):
        sort_dataf = self.dataf + '.sort'
        # sort by is_sc_right, sc_name, is_reverse
        subprocess.call('sort -k 1,1 -i {} -o {}'.format(self.dataf, sort_dataf), shell=True)
        self.sort_dataf = sort_dataf


class FusionType(object):
    DEL = 0
    DUP = 1
    INV = 2
    CTX1 = 3
    CTX2 = 4
    
    TYPE2STR = {DEL:'DEL/ITX',
                DUP:'DUP/ITX',
                INV:'INV',
                CTX1:'CTX-1',
                CTX2:'CTX-2'}
    
    def __init__(self, _type):
        if _type not in self.TYPE2STR:
            raise ValueError('Wrong type number "{}" received'.format(_type))
        self.type = _type
        self.type_str = self.TYPE2STR[_type]
    
    def __repr__(self):
        return self.type_str
    
    @classmethod
    def parse_from_strand(cls, chr1, pos1, is_reverse1, chr2, pos2, is_reverse2):
        if chr1 != chr2:
            if is_reverse1 == is_reverse2:
                # cross type
                return FusionType(cls.CTX2)
            return FusionType(cls.CTX1)
        if is_reverse1 == is_reverse2:
            return FusionType(cls.INV)
        elif (pos1 < pos2) ^ is_reverse1:
            # for RL and pos1 > pos2, it's a DEL
            return FusionType(cls.DEL)
        return FusionType(cls.DUP)


class HotStatus(object):
    REVHOT = 3 # another side of hotfusion, for inv and ctx-2
    HOTFUSION = 4
    
    def __init__(self, num):
        self.num = num
    
    def __repr__(self):
        if self.num < self.REVHOT:
            return str(self.num)
        elif self.num == self.REVHOT:
            return 'RevHot'
        return 'HotFusion'
    
    @property
    def is_hotfusion(self):
        return True if self.num >= self.REVHOT else False


class FusionAnnoHelper(object):
    
    def __init__(self):
        self.regs = {'1':set(), '2':set()}
    
    def add(self, reg_num, lr, type, anno):
        # use tuple to store in set
        self.regs[reg_num].add((lr, type, anno))
    
    @classmethod
    def union_fusion(cls, fah1, fah2):
        for is_rev, (reg_num1, reg_num2) in enumerate([['1', '2'], ['2', '1']]):
            regs1 = fah1.regs[reg_num1]
            regs2 = fah2.regs[reg_num2]
            union_regs = set.intersection(regs1, regs2)
            yield is_rev, union_regs

class PosAnno(object):
    '''
    class to do hot anno and handle the comparaion
    '''
    def __init__(self):
        self.chr2gpos = defaultdict(int)
        self.hotregion = None
        self.chrlen = defaultdict(int)
    
    def load_bam_header(self, header):
        '''
        get global pos from bam header
        '''
        pre_pos = 0
        chr2gpos = defaultdict(int)
        for d in header['SQ']:
            chr2gpos[d['SN']] = pre_pos + d['LN']
            pre_pos += d['LN']
            self.chrlen[d['SN']] = d['LN']
        self.chr2gpos = chr2gpos
    
    def load_hotregion(self, hotregion):
        '''
        load hot region file, which is generated by sv flow
        '''
        self.hotregion = hotregion
    
    def anno2hotregion(self, chr1, pos1, chr2, pos2, is_gene1R, is_gene2R, fusion_type):
        '''
        anno 1-based region to hot region
        '''
        num_hotgene = 0
        if self.hotregion is None:
            return HotStatus(num_hotgene)
        
        hotregion_anno = pysam.TabixFile(self.hotregion)
        fusion = [FusionAnnoHelper(), FusionAnnoHelper()]
        for ind, (c, p) in enumerate([[chr1, pos1], [chr2, pos2]]):
            c = rm_chr_prefix(c)
            try:
                res = list(hotregion_anno.fetch(c, p-1, p))
            except ValueError:
                # for the situation that c not in bed, show 'could not create iterator for region'
                res = []
            
            in_hotgene = False
            for r in res:
                tmp = r.split('\t')
                _type = tmp[3]
                if _type == 'G':
                    in_hotgene = True
                elif _type == 'F':
                    fusion[ind].add(*tmp[4:])
            if in_hotgene:
                num_hotgene += 1
        
        lr = 'R' if is_gene1R else 'L'
        lr += 'R' if is_gene2R else 'L'
        is_revhot = False
        for is_rev, union_regs in FusionAnnoHelper.union_fusion(*fusion):
            for subr in union_regs:
                # maybe in multi hotfusion region
                sub_type = subr[1]
                
                # for DEL in DEL/ITX, use in to judge
                if sub_type in str(fusion_type):
                    if 'INV' in str(fusion_type) or 'CTX' in str(fusion_type):
                        # for CTX-1, lr order should be changed to the input order
                        subr_lr = subr[0][::-1] if is_rev else subr[0]
                        # judge whether revhot
                        if lr == subr_lr:
                            return HotStatus(HotStatus.HOTFUSION)
                        is_revhot = True
                    else:
                        return HotStatus(HotStatus.HOTFUSION)
        if is_revhot:
            # not a hotfustion, but is rev of a hotfusion
            return HotStatus(HotStatus.REVHOT)
        return HotStatus(num_hotgene)
    
    def cmp_chrpos(self, x, y, is_pe=False):
        '''
        Compare the position
        input type for x/y is (refenece, pos)
        if is_pe is True, input type for x/y is (refenece1, pos1, reference2, pos2)
        '''
        chr2gpos = self.chr2gpos
        r = chr2gpos[x[0]] + x[1] - (chr2gpos[y[0]] + y[1])
        if r or not is_pe:
            return r
        return chr2gpos[x[2]] + x[3] - (chr2gpos[y[2]] + y[3])
    
    def cmp_reads(self, x, y, use_mate=False):
        c1 = x.next_reference_name if use_mate else x.reference_name
        p1 = x.next_reference_start if use_mate else x.reference_start
        c2 = y.next_reference_name if use_mate else y.reference_name
        p2 = y.next_reference_start if use_mate else y.reference_start
        return self.cmp_chrpos((c1, p1), (c2, p2))


class MergeFlank(object):
    def __init__(self, sample_name):
        self.data = defaultdict(list)
        self.flank = SvConfig.FIND_FLANK
        self.group_id = 0
        self.sample_name = sample_name
    
    class _MergeData(object):
        _MergeInfo = namedtuple('MergeInfo', 'show_info, total_vars, is_case')
    
        def __init__(self, source, ordered_bp_pos1, ordered_bp_pos2, sc_count, unmap_count, res_info, is_case, order_rev):
            self.source, self.ordered_bp_pos1, self.ordered_bp_pos2 = source, ordered_bp_pos1, ordered_bp_pos2
            self.sc_count, self.unmap_count, self.res_info, self.is_case, self.order_rev = sc_count, unmap_count, res_info, is_case, order_rev
            self.template_count = res_info.template_count
            
            self.main_bp = False
            self.two_side_sr = [0, 0]
            self.two_side_su = [0, 0]
            self.is_uniq = True
            self.group_id = '.'
        
        def get_show_info(self, sample_name):
            [source, hot_status, fusion_type, chr1, bp_pos1, start1, end1, chr2, bp_pos2, start2, end2,
             is_gene1R, is_reverse1, is_gene2R, is_reverse2,
             sr_segcigar, su_segcigar, ori_sc_count, ori_unmap_count, pe_count, template_count, mapq] = self.res_info
            
            # adjust SR/SU which merged adjacent bps
            sc_count = max(self.two_side_sr)
            unmap_count = max(self.two_side_su)
            
            ori_alt_depth = ori_sc_count + ori_unmap_count + pe_count
            # adjusted count:
            # double SR which is the another side SR, make the assumption that the SR will be symmetrical
            alt_depth = 2*sc_count + 2*unmap_count + pe_count
            
            #! split reads also coverage the breakpoint
            normal_depth = template_count - ori_sc_count
            
            ori_depth = ori_alt_depth + normal_depth
            depth = alt_depth + normal_depth
            
            if start1 == end1:
                posfmt1 = '{}:{}'.format(chr1, start1)
            else:
                posfmt1 = '{}:{}-{}'.format(chr1, start1, end1)
            if start2 == end2:
                posfmt2 = '{}:{}'.format(chr2, start2)
            else:
                posfmt2 = '{}:{}-{}'.format(chr2, start2, end2)
            
            link = 'http://10.0.1.7:4384/cnv-cgi/sv.cgi?s={}&p1={}:{}&p2={}:{}&lr={}{}'.format(
                sample_name, chr1, bp_pos1, chr2, bp_pos2, 'R' if is_gene1R else 'L', 'R' if is_gene2R else 'L')
            
            #TODO exclude other sr count at same location
            
            show_info = [hot_status,
                         fusion_type,
                         self.group_id,
                         'Y' if self.main_bp else 'N',
                         'Y' if self.is_case else 'N',
                         'Y' if self.is_uniq else 'N',
                         mapq,
                         'Y' if mapq >= SvConfig.MAPQ_THR else 'N',
                         '{:.4f}'.format(ori_alt_depth/ori_depth),
                         '{:.4f}'.format(alt_depth/depth),
                         posfmt1, posfmt2,
                         source,
                         '{0}{1}|{2}{3}'.format(*['R' if is_gene1R else 'L', 
                                                  '-' if is_reverse1 else '+',
                                                  'R' if is_gene2R else 'L',
                                                  '-' if is_reverse2 else '+',
                                                 ]), 
                         sr_segcigar,
                         su_segcigar,
                         ori_sc_count,
                         ori_unmap_count,
                         sc_count,
                         unmap_count,
                         pe_count,
                         alt_depth,
                         depth,
                         link,]
            return self._MergeInfo(show_info, alt_depth, self.is_case)
    
    def add(self, res_info, is_case=True):
        ''' add the data by the given order.
        '''
        source, chr1, bp_pos1, chr2, bp_pos2 = res_info.source, res_info.chr1, res_info.bp_pos1, res_info.chr2, res_info.bp_pos2
        is_reverse1, is_reverse2 = res_info.is_reverse1, res_info.is_reverse2
        # if source == 'SplitRead':
            # count = res_info.sc_count
        # elif source == 'SingleUnmap':
            # count = res_info.unmap_count
        # elif source == 'DiscordantPair':
            # count = res_info.pe_count
        # else:
            # raise ValueError("Can't recognize the source {}".format(source))
        sc_count, unmap_count = res_info.sc_count, res_info.unmap_count
        
        if str(chr1) == str(chr2):
            if bp_pos2 < bp_pos1:
                key = '{}_{}_{}_{}'.format(chr2, chr1, is_reverse2, is_reverse1)
                # last two is mainBP and OtherBpCount
                self.data[key].append(self._MergeData(source, bp_pos2, bp_pos1, sc_count, unmap_count, res_info, is_case, True))
            else:
                key = '{}_{}_{}_{}'.format(chr1, chr2, is_reverse1, is_reverse2)
                self.data[key].append(self._MergeData(source, bp_pos1, bp_pos2, sc_count, unmap_count, res_info, is_case, False))
        else:
            if str(chr2) < str(chr1):
                key = '{}_{}_{}_{}'.format(chr2, chr1, is_reverse2, is_reverse1)
                self.data[key].append(self._MergeData(source, bp_pos2, bp_pos1, sc_count, unmap_count, res_info, is_case, True))
            else:
                key = '{}_{}_{}_{}'.format(chr1, chr2, is_reverse1, is_reverse2)
                self.data[key].append(self._MergeData(source, bp_pos1, bp_pos2, sc_count, unmap_count, res_info, is_case, False))
    
    def load_ctrl(self, ctrl_data):
        for k, d in ctrl_data.data.items():
            self.data[k].extend(d)
    
    def merge_bp(self, merge_data=None):
        if merge_data is None:
            for chr_data in self.data.values():
                merge_data = []
                pre_pos = -1
                chr_data.sort(key=lambda _: _.ordered_bp_pos1) # sort by ordered_bp_pos1
                for d in chr_data:
                    if d.ordered_bp_pos1-pre_pos > self.flank:
                        self.merge_bp(merge_data)
                        merge_data = [d]
                        pre_pos = d.ordered_bp_pos1
                    else:
                        merge_data.append(d) # pass the id
                self.merge_bp(merge_data)
        else:
            if not merge_data:
                # no data
                return
            group = []
            merge_data.sort(key=lambda _: _.ordered_bp_pos2) # sort by ordered_bp_pos2
            for d in merge_data:
                for g in group:
                    g0 = g[0]
                    if abs(g0.ordered_bp_pos1-d.ordered_bp_pos1) <= self.flank and abs(g0.ordered_bp_pos2-d.ordered_bp_pos2) <= self.flank:
                        g.append(d)
                        break
                else:
                    # no break
                    group.append([d])
            for g in group:
                case_data = [d for d in g if d.is_case]
                ctrl_data = [d for d in g if not d.is_case]
                is_uniq = False if (case_data and ctrl_data) else True
                for case_ctrl_data in [case_data, ctrl_data]:
                    group_id = '.'
                    if not case_ctrl_data:
                        # no data
                        continue
                    if len(case_ctrl_data) > 1:
                        case_ctrl_data.sort(key=lambda _: _.template_count, reverse=True)
                        group_id = self.group_id
                        self.group_id += 1
                    
                    # count the support at each side
                    two_side_sr = [0, 0]
                    two_side_su = [0, 0]
                    
                    for order_rev in [False, True]:
                        base_sr_pos = None
                        for d in case_ctrl_data:
                            if d.order_rev != order_rev:
                                continue
                            
                            d.is_uniq = is_uniq
                            d.group_id = group_id
                            if d.source == 'SplitRead':
                                if base_sr_pos is None:
                                    base_sr_pos = d.ordered_bp_pos2 if order_rev else d.ordered_bp_pos1
                                    two_side_sr[order_rev] += d.sc_count
                                else:
                                    sr_pos = d.ordered_bp_pos2 if order_rev else d.ordered_bp_pos1
                                    if sr_pos != base_sr_pos:
                                        #TODO, In merge SR same p1, but different p2 will share the same base count
                                        # If allow the same p1, will make count bigger than real, and vice verca
                                        two_side_sr[order_rev] += d.sc_count
                            elif d.source == 'SingleUnmap':
                                two_side_su[order_rev] += d.unmap_count
                            elif d.source == 'DiscordantPair':
                                # add the refetch ones
                                two_side_sr[order_rev] += d.sc_count
                                two_side_su[order_rev] += d.unmap_count
                    
                    for d in case_ctrl_data:
                        d.two_side_sr = two_side_sr
                        d.two_side_su = two_side_su
                    
                    # set the main one as the largest template count
                    case_ctrl_data[0].main_bp = True
    
    def write(self, res_file, RES_TITLE, only_ctrl=False):
        show_infos = []
        for d in self.data.values():
            show_infos.extend([_.get_show_info(self.sample_name) for _ in d])
        show_infos.sort(key=lambda _: _.total_vars, reverse=True) # sort by total_vars
        with open(res_file, 'w') as outfile:
            outfile.write(','.join(RES_TITLE)+'\n')
            for r in show_infos:
                if only_ctrl and r.is_case:
                    # skip the is_case = Y item
                    continue
                outfile.write(','.join(map(str, r.show_info))+'\n')

def _rec_memory(pid, memf):
    import re, time, datetime
    header = []
    have_header = False
    with open(memf, 'w') as outfile:
        while True:
            data_list = []
            for line in open('/proc/{}/status'.format(pid)):
                if not line.startswith('Vm'):
                    continue
                tit, data = line.strip().split(':')
                if not have_header:
                    header.append(tit)
                try:
                    data_list.append(str(int(float(re.findall(r'(\d+)\skB', data)[0]) /1024)) + ' Mb')
                except Exception as e:
                    data_list.append(data.strip())
            if not have_header:
                have_header = True
                outfile.write(','.join(['Time'] + header) + '\n')
            outfile.write(','.join([str(datetime.datetime.now())] + data_list) + '\n')
            outfile.flush()
            time.sleep(30)

def get_memory(res_file):
    from multiprocess import Pool
    memf = os.path.join(os.path.dirname(res_file), 'sv_memory.log')
    pid = os.getpid()
    p = Pool(1)
    p.apply_async(_rec_memory, [pid, memf])
    p.close()
    

############### MAIN #################
if __name__ == '__main__':
    # setup logging
    logging.basicConfig(format="[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)s] %(message)s", level=logging.INFO)
    
    # SplitRead: SR
    # SingleUnmap: SU
    # DiscordantPair: DP
    
    RES_TITLE = ['NumOfHot', 'FusionType', 'MergeBP', 'MainBP', 'IsCase', 'IsUnique', 'Qual', 'PassQual',
                 'OriFreq', 'Freq', 'Pos1', 'Pos2', 'Source', 'BreakPointInfo', 'SRSeqCigar', 'SUSeqCigar',
                 'OriSRCount', 'OriSUCount', 'SRCount', 'SUCount', 'DPCount', 'AltDepth', 'Depth', 'Link']
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tmp', type=str, required=True,
                        help='tmp file path')
    parser.add_argument('-r', '--ref', type=str, required=True,
                        help='reference genome, it should be same with the one used by bam')
    parser.add_argument('-H', '--hotregion', type=str, required=True,
                        help='hot region file, generated by sv flow helper script')
    parser.add_argument('-w', '--bwa', type=str, default='bwa',
                        help='default is bwa')
    parser.add_argument('-s', '--samtools', type=str, default='samtools',
                        help='default is samtools')
    parser.add_argument('-n', '--numthreads', type=int, default=4,
                        help='num of threads, default is 4')
    parser.add_argument('-c', '--ctrlout', type=str, default=None,
                        help='norm sv output file, csv format')
    parser.add_argument('-S', '--sample', type=str, default=None,
                        help='sample name, set None to infer it from bam path')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='open debug logging')
    parser.add_argument('--bamstat', type=str, default=None,
                        help='bam stat file, isize could be abstracted from it')
    parser.add_argument('outfile',
                        help='output file, csv format')
    parser.add_argument('bamt', 
                        help='bam to do sv analysis')
    parser.add_argument('bamn', nargs='?',
                        help='ctrl bam')
    in_args = parser.parse_args()
    res_file, tmp_path, bamn, bamt, ref = in_args.outfile, in_args.tmp, in_args.bamn, in_args.bamt, in_args.ref
    bwa, samtools, hotregion, num_threads, = in_args.bwa, in_args.samtools, in_args.hotregion, in_args.numthreads
    ctrl_res_file, sample_name, is_debug = in_args.ctrlout, in_args.sample, in_args.debug
    bam_statf = in_args.bamstat
    
    logging.info('ncsv Ver 0.2.3')
    logging.info('Run Args:\n{} {} {}'.format(sys.executable,
                                              os.path.abspath(sys.argv[0]),
                                              ' '.join(sys.argv[1:])))
    
    if is_debug:
        from guppy import hpy
        get_memory(res_file)
    start_time = time.time()
    
    try:
        ctrl_data = None
        args = [tmp_path, ref, bwa, samtools]
        kwargs = {'num_threads':num_threads,
                  'hotregion':hotregion,
                  'bam_statf':bam_statf,
                  'sample_name':sample_name,
                  'is_debug':is_debug,
                  'start_time':start_time
                 }
        if bamn:
            ctrl_data = run_sv('', bamn, *args, **kwargs)
        run_sv(res_file, bamt, *args, ctrl_data=ctrl_data, ctrl_res_file=ctrl_res_file, **kwargs)
    except Exception as e:
        logging.error('Run_sv Error: {}\n{}\n'.format(e, traceback.format_exc()))
        sys.exit(1) 