import argparse 
from Bio import SeqIO
import pysam
import numpy as np
import scipy as sp
import re
from pysam import AlignmentFile as BAM
from pysam import VariantFile as VCF
from pysam import index as indexBam
from pysam import TabixFile as BED
import gzip
#from pysam import FastaFile as FASTA
import time
from DCutils.utils import extractDepthSnv,extractDepthIndel

"""
def calculateMismatches(seqs,referenceSequence):
    mismatch = 0
    for seq in seqs:
        sequence = seq.query_sequence
        for ii in range(len(sequence)):
            if sequence[ii] != referenceSequence[ii]:
                mismatch += 1
                break
    return mismatch
"""

def log10(mat):
    return np.log10(np.where(mat>0,mat,np.finfo(float).eps))

def power10(mat):
    return 10**np.where(mat>=np.log10(np.finfo(float).eps),mat,-np.inf)

def log2(mat):
    return np.log2(np.where(mat>0,mat,np.finfo(float).eps))

def power2(mat):
    return np.where(mat>=-100,2**mat,0)

def log(mat):
    return np.log(np.where(mat>0,mat,np.finfo(float).eps))

def prepare_reference_mats(chrom,start,end,reference_seq,germline_bed,noise_bed,n_cov_bed,params):
    ### Define and Initialize
    af_miss = params['mutRate']
    m = end - start
    base2num = {'A':0,'T':1,'C':2,'G':3}
    snp_mask = np.full(m,False,dtype=bool)
    indel_mask = np.full(m,False,dtype=bool)
    noise_mask = np.full(m,False,dtype=bool)
    n_cov_mask = np.full(m,False,dtype=bool)
    
    ### Initialize prior mat as no germline resource
    reference_int = np.array([base2num.get(b,4) for b in reference_seq],dtype=int)
    noise_mask[reference_int == 4] = True
    reference_int[reference_int == 4] = 0
    prior_mat = np.full([m,4],af_miss,dtype=float)
    prior_mat[np.ogrid[:m],reference_int] = 1 - 3 * af_miss

    ### Adjust by germline
    for rec in germline_bed.fetch(chrom,start,end):
        ind = rec.pos - 1 - start
        ref = rec.ref
        afs = rec.info['AF']
        if len(ref) == 1:
            has_snp = False
            for ii,alt in enumerate(rec.alts):
                if len(alt) == 1:
                    prior_mat[ind,base2num[alt]] = afs[ii]
                    has_snp = True
                    if afs[ii] >= 0.01:
                        snp_mask[ind] = True
                elif afs[ii] >= 0.01: 
                    indel_mask[ind] = True
            if has_snp:
                prior_mat[ind,base2num[ref]] = 1 - np.delete(prior_mat[ind,:],base2num[ref]).sum()
        if len(ref) != 1:
            for ii,alt in enumerate(rec.alts):
                if afs[ii] >= 0.01:
                    indel_mask[ind+1:ind+len(alt)] = True
    
    ### Preparep noise mask
    for rec in noise_bed.fetch(chrom,start,end,parser=pysam.asBed()):
        interval_start = max(rec.start,start)
        interval_end = min(rec.end,end)
        interval_len = interval_end - interval_start
        interval_start_ind = interval_start - 1 - start
        noise_mask[interval_start_ind:interval_start_ind + interval_len] = True

    ### Preparep normal coverage mask
    for rec in n_cov_bed.fetch(chrom,start,end,parser=pysam.asTuple()):
        depth = int(rec[3])
        if depth >= params["minNdepth"]:
            continue
        interval_start = max(int(rec[1]),start)
        interval_end = min(int(rec[2]),end)
        interval_len = interval_end - interval_start
        interval_start_ind = interval_start - 1 - start
        n_cov_mask[interval_start_ind:interval_start_ind + interval_len] = True
    return prior_mat,snp_mask,indel_mask,noise_mask,n_cov_mask,reference_int

def determineTrimLength(seq,params,processed_flag):
    if seq.is_forward:
        overlap = 0 # Never mask overlap of forward read
        left = params['trim5']
        right_frag = params['trim5'] - min(params['trim5'],abs(seq.template_length) - seq.reference_length)
        right_read = params['trim3']
        right = max(right_frag,right_read)
    else:
        ### Mask overlap of reverse read
        if processed_flag:
            mate_cigar = seq.get_tag("MC")
            re_cigar =  re.search('(?:(\d+)S)?(\d+)M(?:(\d+)S)?',mate_cigar)
            mate_reference_length = int(re_cigar.group(2))
            overlap = max(0,seq.reference_length + mate_reference_length - abs(seq.template_length))
        else: overlap =0
        right_frag = params['trim5']
        right_read = params['trim3']
        right=max(right_frag,right_read)
        left_frag = params['trim5'] - min(params['trim5'],abs(seq.template_length) - seq.reference_length)
        left=max(left_frag,overlap,params['trim3'])
    return left,right


def genotypeDSSnv(seqs,reference_int,prior_mat,antimask,params):
    params["print_flag"]=0
    np.set_printoptions(edgeitems=30, linewidth=10000,suppress=True,)
    prob_damage = params['amperr']
    prob_damage_mat = np.ones([4,4],dtype=float)*prob_damage/3
    np.fill_diagonal(prob_damage_mat,1-prob_damage)
    log_prob_damage_mat = log10(prob_damage_mat)
    prob_amp = params['amperr']
    prob_amp_mat = np.ones([4,4],dtype=float)*prob_amp/3
    np.fill_diagonal(prob_amp_mat,1-prob_amp)
    log_prob_amp_mat = log10(prob_amp_mat)

    #seq_mat = 1-power10(-qual_mat_merged/10).reshape(4,n,1)

    ### Assign reads to each strand
    F1R2 = []
    F2R1 = []
    for seq in seqs:
        if (seq.is_read1 and seq.is_forward) or (seq.is_read2 and seq.is_reverse):
            F1R2.append(seq)
        if (seq.is_read2 and seq.is_forward) or (seq.is_read1 and seq.is_reverse):
            F2R1.append(seq)   
    

    ### Determine match length

    cigar = seqs[0].cigarstring
    re_cigar =  re.search('(?:(\d+)S)?(\d+)M(?:(\d+)S)?',cigar)
    leftS = re_cigar.group(1) or 0
    matches = re_cigar.group(2)
    rightS = re_cigar.group(3) or 0

    n = min([seq.reference_length for seq in seqs])
    m_F1R2 = len(F1R2)
    m_F2R1 = len(F2R1)

    ### Prepare sequence matrix and quality matrix for each strand
    base2num = {'A':0,'T':1,'C':2,'G':3}
    F1R2_seq_mat = np.zeros([m_F1R2,n],dtype=int) # Base(ATCG) x reads x pos
    F1R2_qual_mat = np.zeros([m_F1R2,n])
    F2R1_seq_mat = np.zeros([m_F2R1,n],dtype=int) # Base(ATCG) x reads x pos
    F2R1_qual_mat = np.zeros([m_F2R1,n])
    
    
    for mm,seq in enumerate(F1R2):
        qualities = seq.query_alignment_qualities
        sequence = seq.query_alignment_sequence
        for nn in range(n):
            base = sequence[nn]
            if base == 'N':
                F1R2_seq_mat[mm,nn] = 0
                F1R2_qual_mat[mm,nn] = np.log10(0.25)
            else:
                F1R2_seq_mat[mm,nn] = base2num[base]
                F1R2_qual_mat[mm,nn] = qualities[nn]
    
    F1R2_qual_mat_merged = np.zeros([4,n])
    F1R2_count_mat = np.zeros([4,n],dtype=int)
    for nn in range(0,4):
        F1R2_qual_mat_merged[nn,:] = F1R2_qual_mat.sum(axis=0,where=(F1R2_seq_mat==nn))
        F1R2_count_mat[nn,:] = (F1R2_seq_mat==nn).sum(axis=0) 

    for mm,seq in enumerate(F2R1):
        qualities = seq.query_alignment_qualities
        sequence = seq.query_alignment_sequence
        for nn in range(n):
            base = sequence[nn]
            if base == 'N':
                F2R1_seq_mat[mm,nn] = 0
                F2R1_qual_mat[mm,nn] = np.log10(0.25)
            else:
                F2R1_seq_mat[mm,nn] = base2num[base]
                F2R1_qual_mat[mm,nn] = qualities[nn]
    F2R1_qual_mat_merged = np.zeros([4,n])
    F2R1_count_mat = np.zeros([4,n],dtype=int)
    for nn in range(0,4):
        F2R1_qual_mat_merged[nn,:] = F2R1_qual_mat.sum(axis=0,where=(F2R1_seq_mat==nn))
        F2R1_count_mat[nn,:] = (F2R1_seq_mat==nn).sum(axis=0)

    total_count_mat = F1R2_count_mat + F2R1_count_mat
    total_count_mat[reference_int,np.ogrid[:n]] = -1
    antimask[(total_count_mat >= 1).sum(axis=0) > 1] = False
    alt_int = np.argmax(total_count_mat,axis=0)


    F1R2_masked_qual_mat = F1R2_qual_mat_merged[:,antimask]
    F2R1_masked_qual_mat = F2R1_qual_mat_merged[:,antimask]
    reference_int_masked = reference_int[antimask]
    alt_int_masked = alt_int[antimask]
    Pamp = log_prob_damage_mat[reference_int_masked,alt_int_masked]
    
    mmm = alt_int_masked.size
    prior_ref = log10(prior_mat[antimask][np.ogrid[:mmm],reference_int_masked])
    prior_alt = log10(prior_mat[antimask][np.ogrid[:mmm],alt_int_masked])
    
    F1R2_ref_prob_mat = -F1R2_masked_qual_mat[reference_int_masked,np.ogrid[:reference_int_masked.size]]/10
    F1R2_alt_prob_mat = -F1R2_masked_qual_mat[alt_int_masked,np.ogrid[:reference_int_masked.size]]/10
    F2R1_ref_prob_mat = -F2R1_masked_qual_mat[reference_int_masked,np.ogrid[:reference_int_masked.size]]/10
    F2R1_alt_prob_mat = -F2R1_masked_qual_mat[alt_int_masked,np.ogrid[:reference_int_masked.size]]/10


    def calculatePosterior(Pamp,Pref,Palt,prior_ref,prior_alt):
        log_e_10 = np.log(10)
        Pamp = power10(Pamp)
        Pref = power10(Pref)
        Palt = power10(Palt)

        first_term = log10(1 + 2*Pamp*Pref - Pamp - Pref)
        second_term = log10(Pamp + Palt - Pamp*Palt)
        ref_prob = prior_ref + first_term + second_term

        first_term = log10(1 + 2*Pamp*Palt - Pamp - Palt)
        second_term = log10(Pamp + Pref - Pamp*Pref)
        alt_prob = prior_alt + first_term + second_term
        return ref_prob,alt_prob
    #print(prior_mat)
    F1R2_ref_prob,F1R2_alt_prob = calculatePosterior(Pamp,F1R2_ref_prob_mat,F1R2_alt_prob_mat,prior_ref,prior_alt)
    F2R1_ref_prob,F2R1_alt_prob = calculatePosterior(Pamp,F2R1_ref_prob_mat,F2R1_alt_prob_mat,prior_ref,prior_alt)

    #TS_gt = np.ones([4,n]) * 0.25
    #BS_gt = np.ones([4,n]) * 0.25
    F1R2_ARLR_masked = F1R2_alt_prob - F1R2_ref_prob
    F2R1_ARLR_masked = F2R1_alt_prob - F2R1_ref_prob

    F1R2_ARLR = np.zeros(n)
    F2R1_ARLR = np.zeros(n)
    F1R2_ARLR[antimask] = F1R2_ARLR_masked
    F2R1_ARLR[antimask] = F2R1_ARLR_masked
    #print(F1R2_ref_qual_mat)
    #print(F1R2_alt_qual_mat)
    return F1R2_ARLR,F2R1_ARLR,reference_int,alt_int,antimask,F1R2_count_mat,F2R1_count_mat

def genotypeSSSnv(seq_mat,qual_mat,prior_mat,params,A_log_prob,T_log_prob,C_log_prob,G_log_prob):
    np.seterr(all="raise")
    #qual_mat[qual_mat<=15] = 0
    ### Define and Initialize
    np.fill_diagonal(prob_amp_mat,1-prob_amp)
    m,n = seq_mat.shape
    base2num = {'A':0,'T':1,'C':2,'G':3}
    ### Merge quality 
    qual_mat_merged = np.zeros([4,n])
    count_mat = np.zeros([4,n],dtype=int)
    for nn in range(0,4):
        qual_mat_merged[nn,:] = qual_mat.sum(axis=0,where=(seq_mat==nn))
        count_mat[nn,:] = (seq_mat==nn).sum(axis=0) + 1
    qual_mat_merged = qual_mat_merged - np.log(0.25)*10
    #confident = ((count_mat.sum(axis=0) - count_mat.max(axis=0)) == 3)
    #qual_mat_merged_ambiguous = qual_mat_merged[~confident]
    A_probs = np.zeros([16,n])
    T_probs = np.zeros([16,n])
    C_probs = np.zeros([16,n])
    G_probs = np.zeros([16,n])

    #A_log_prob_readA = 
    #A_log_prob_A = A_log_prob[[0,4,8,12],:] + qual_mat_merged



def genotypeDSIndel(seqs,bam,params):
    np.set_printoptions(edgeitems=30, linewidth=10000,suppress=True)
    prob_damage = params['damage']/100
    log_prob_damage = np.log10(prob_damage)
    log_prob_nodamage = np.log10(1 - prob_damage)
    #log_prior = log10(prior_mat)
    ### Assign reads to each strand
    F1R2 = []
    F2R1 = []
    for seq in seqs:
        if (seq.is_read1 and seq.is_forward) or (seq.is_read2 and seq.is_reverse):
            F1R2.append(seq)
        if (seq.is_read2 and seq.is_forward) or (seq.is_read1 and seq.is_reverse):
            F2R1.append(seq)   
    chrom = seqs[0].reference_name
    start = seqs[0].reference_start
    end = seqs[0].reference_end
    f1r2_names_dict = {f1r2.query_name:nn for nn,f1r2 in enumerate(F1R2)}
    f2r1_names_dict = {f2r1.query_name:nn for nn,f2r1 in enumerate(F2R1)}
    f1r2_genotype = np.zeros([len(F1R2),1])
    f2r1_genotype = np.zeros([len(F2R1),1])
    for nn,col in enumerate(bam.pileup(chrom,start,end,flag_filter=2816)):
        indel_flag = 0
        for read in col.pileups:
            if read.alignment.query_name in f1r2_names_dict:
                if read.indel != 0:
                    indel_flag = 1
                    f1r2_genotype[f1r2_names_dict[read.alignment.query_name],0] = read.indel
                print(read.indel,read.alignment.query_name,f1r2_names_dict)
        if indel_flag == 1:
            break
    for nn,col in enumerate(bam.pileup(chrom,start,end,flag_filter=2816)):
        indel_flag = 0
        for read in col.pileups:
            if read.alignment.query_name in f1r2_names_dict:
                if read.indel != 0:
                    indel_flag = 1
                    f1r2_genotype[f1r2_names_dict[read.alignment.query_name],0] = read.indel
                print(read.indel,read.alignment.query_name,f1r2_names_dict)
        if indel_flag == 1:
            break         
    """
    for mm,seq in enumerate(F1R2):
        qualities = seq.query_alignment_qualities
        sequence = seq.query_alignment_sequence
        for nn in range(n):
            base = sequence[nn]
            if base == 'N':
                F1R2_seq_mat[mm,nn] = 0
                F1R2_qual_mat[mm,nn] = np.log10(0.25)
            else:
                F1R2_seq_mat[mm,nn] = base2num[base]
                F1R2_qual_mat[mm,nn] = qualities[nn]
    
    for mm,seq in enumerate(F2R1):
        qualities = seq.query_alignment_qualities
        sequence = seq.query_alignment_sequence
        for nn in range(n):
            base = sequence[nn]
            if base == 'N':
                F2R1_seq_mat[mm,nn] = 0
                F2R1_qual_mat[mm,nn] = np.log10(0.25)
            else:
                F2R1_seq_mat[mm,nn] = base2num[base]
                F2R1_qual_mat[mm,nn] = qualities[nn]
    """

    log_gt_mat = np.ones([n,4]) * np.log(0.25)
    F1R2_gt = genotypeSSSnv(F1R2_seq_mat[:,antimask],F1R2_qual_mat[:,antimask],prior_mat[antimask,:],params)
    F2R1_gt = genotypeSSSnv(F2R1_seq_mat[:,antimask],F2R1_qual_mat[:,antimask],prior_mat[antimask,:],params)
    m,n = F1R2_gt.shape
    #pairwise = np.repeat(F1R2_gt.reshape((m,n,1)),4,axis=2) + np.repeat(F2R1_gt.reshape((m,1,n)),4,axis=1)
    log_gt_mat_masked = np.zeros([m,4])
    consensus_gt = np.logical_and(F1R2_gt >= 0.95,F2R1_gt >= 0.95,dtype=int)
    
    """
    for nn in range(0,4):
        ind_no_nn = [0,1,2,3]
        ind_no_nn.remove(nn)
        damage_prob_mat = np.zeros([4,4])
        damage_prob_mat[nn,nn] = 2 * log_prob_nodamage 
        damage_prob_mat[ind_no_nn,nn] = log_prob_nodamage + log_prob_damage
        damage_prob_mat[nn,ind_no_nn] = log_prob_nodamage + log_prob_damage
        damage_prob_mat[ind_no_nn,ind_no_nn] = 2 * log_prob_damage 
        pairwise_nn = pairwise + np.repeat(damage_prob_mat.reshape([1,4,4]),m,axis=0)
        log_per_base_prob = log10(((power10(pairwise_nn)).sum(axis=2).sum(axis=1)))
        log_gt_mat_masked[:,nn] = log_per_base_prob + log_prior[antimask,nn]
    """
    #gt_mat_masked = power10(log_gt_mat_masked)
    #log_gt_mat_masked = log10(gt_mat_masked / gt_mat_masked.sum(axis=1,keepdims=True))
    gt_mat[antimask,:] = consensus_gt
    return gt_mat

def bamIterateMultipleRegion(bamObject,regions):
    for region in regions:
        for rec in bamObject.fetch(*region):
            if len(region)  >= 2:
                if rec.reference_start < region[1]:
                    continue
            yield rec,region

def callBam(params,processNo,chunkSize):
    ### Get parameters
    bam = params['tumorBam']
    nbam = params['normalBam']
    regions = params['regions']
    germline = VCF(params['germline'],'rb')
    ncoverage = BED(params['ncoverage'])
    reference = params['reference']
    minMapq = params["mapq"]
    mutRate = params["mutRate"]
    pcut = params["pcutoff"]
    nn = processNo
    output = 'tmp/' + params["output"] + '_' + str(nn)
    noise= BED(params['noise'])
    base2num = {'A':0,'T':1,'C':2,'G':3}
    num2base = 'ATCG'
    muts = []
    muts_dict = dict()
    duplex_read_num_dict = dict()

    print("Process" + str(processNo)+": Initiated")
    ### Initialize
    total_coverage = 0
    starttime = time.time()
    tumorBam = BAM(bam,'rb')
    normalBam = BAM(nbam,'rb')
    currentReadSet = []
    currentBc = ''
    currentStart = -1
    currentReadDict = {}
    #chromPast = 'startChrom'
    fasta = reference
    recCount = 0
    currentCheckPoint = 100000
    duplex_count = 0
    reference_mat_chrom = 'anyChrom'
    reference_mat_start = 0
    locus_bed = gzip.open(output+"_coverage.bed.gz",'wb')
    processed_read_names = set()

    for rec,region in bamIterateMultipleRegion(tumorBam,regions):
        recCount += 1
        if recCount == currentCheckPoint:
            print("Process"+str(processNo)+": processed "+str(recCount)+" reads in "+\
                str((time.time()-starttime)/60)+" minutes" )
            currentCheckPoint += 100000
        if rec.mapping_quality <= minMapq or \
        rec.is_supplementary or \
        rec.is_secondary or \
        rec.has_tag('DT') or \
        not rec.is_proper_pair or \
        rec.is_qcfail:
            continue
        if rec.get_tag('AS') - rec.get_tag('XS') <= 50: continue
        #if rec.cigartuples[0][0] == 4: continue
        if rec.cigarstring.count('I') + rec.cigarstring.count('D') >= 2: continue
        start = rec.reference_start
        bc = rec.query_name.split('_')[1]
        bcsplit = bc.split('+')
        bc1 = bcsplit[0]
        bc2 = bcsplit[1]     
        chrom = tumorBam.get_reference_name(rec.reference_id) 
        if currentStart == -1: currentStart = start 
        #print(start,currentStart)
        if start == currentStart:
            #print(currentReadDict,1)
            if currentReadDict.get(bc1+'+'+bc2) != None:
                currentReadDict[bc1+'+'+bc2]["seqs"].append(rec)
                if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
                    currentReadDict[bc1+'+'+bc2]["F1R2"] += 1
                else:
                    currentReadDict[bc1+'+'+bc2]["F2R1"] += 1           
            elif currentReadDict.get(bc2+'+'+bc1) != None:
                currentReadDict[bc2+'+'+bc1]["seqs"].append(rec)
                if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
                    currentReadDict[bc2+'+'+bc1]["F1R2"] += 1
                else:
                    currentReadDict[bc2+'+'+bc1]["F2R1"] += 1    
            else:
                currentReadDict.update({bc:{"seqs":[rec],"F1R2":0,"F2R1":0}})
                if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
                    currentReadDict[bc1+'+'+bc2]["F1R2"] += 1
                else:
                    currentReadDict[bc1+'+'+bc2]["F2R1"] += 1    

        else:
            """
            Calling block starts
            """
            for key in currentReadDict.keys():
                readSet = currentReadDict[key]["seqs"]
                setBc = key.split('+')
                setBc1 = setBc[0]
                setBc2 = setBc[1]
                F2R1 = currentReadDict[key]["F2R1"]
                F1R2 = currentReadDict[key]["F1R2"]
                duplex_no = f"{min([F1R2,F2R1])}+{max([F1R2,F2R1])}"
                if duplex_read_num_dict.get(duplex_no):
                    duplex_read_num_dict[duplex_no] += 1
                else:
                    duplex_read_num_dict[duplex_no] = 1
                if setBc1 != setBc2 and F2R1 > 1 and F1R2 > 1:
                    #if dupPast == 'startChrom':
                        #chromPast = tumorBam.get_reference_name(rec.reference_id)
                    #else:
                        #chromPast = chrom
                    #chrom = tumorBam.get_reference_name(rec.reference_id)
                    ### Determine if there is any indels in reads
                    indel_bool = [('I' in seq.cigarstring or 'D' in seq.cigarstring) for seq in readSet]
                    if any(indel_bool): 
                        mlmlml=1
                        #genotypeDSIndel(readSet,tumorBam,params)

                    





                    else:
                        rs_reference_end = max([r.reference_end for r in readSet])
                        chromNow = readSet[0].reference_name
                        if chromNow != reference_mat_chrom or rs_reference_end >= reference_mat_end:
                            ### Output coverage 
                            if "coverage" in locals():
                                non_zero_positions = np.nonzero(coverage)
                                for pos in non_zero_positions[0].tolist():
                                    locus_bed.write(('\t'.join([reference_mat_chrom,str(pos),str(pos+1),str(coverage[pos])])+"\n").encode("utf-8"))
                                total_coverage += np.sum(coverage)
                            reference_mat_chrom = chromNow
                            reference_mat_start = readSet[0].reference_start
                            try: region_end = region[2]
                            except: region_end = 10E10
                            contig_end = tumorBam.get_reference_length(chromNow)
                            reference_mat_end = min(readSet[0].reference_start + 100000,max(region_end,readSet[0].reference_end),contig_end)
                            try:prior_mat,snp_mask,indel_mask,noise_mask,n_cov_mask,ref_np = prepare_reference_mats(reference_mat_chrom,\
                                                                                                reference_mat_start,\
                                                                                                reference_mat_end,fasta[reference_mat_chrom][reference_mat_start:reference_mat_end],\
                                                                                                germline,\
                                                                                                noise,\
                                                                                                ncoverage,\
                                                                                                params)
                            except: print(readSet[0].reference_start + 100000,region_end,readSet[0].reference_end,contig_end,reference_mat_start,chromNow,rec.reference_name,readSet[0].reference_name)
                            coverage = np.zeros(100000,dtype=int)
                        ### Record read names to check if mate has been processed
                        processed_flag = 0
                        for seq in readSet:
                            if seq.query_name in processed_read_names:
                                processed_read_names.remove(seq.query_name)
                                processed_flag = 1
                                break
                        if processed_flag == 0:
                            processed_read_names.add(readSet[0].query_name)

                        start_ind = readSet[0].reference_start - reference_mat_start
                        reference_length = min([read.reference_length for read in readSet])
                        end_ind = readSet[0].reference_start + reference_length - reference_mat_start
                        masks = np.zeros([4,end_ind-start_ind],dtype=bool)
                        try:masks[0,:] = snp_mask[start_ind:end_ind]
                        except: print(start_ind,end_ind)
                        masks[1,:] = noise_mask[start_ind:end_ind]
                        masks[2,:] = n_cov_mask[start_ind:end_ind]
                        left,right = determineTrimLength(readSet[0],params=params,processed_flag=processed_flag)
                        masks[3,:left] = True
                        masks[3,-right:] = True 
                        antimask = np.all(~masks,axis=0)
                        ### If the whole reads are masked:
                        if not np.any(antimask): continue
                        ### Calculate genotype probability
                        F1R2_ARLR,F2R1_ARLR,ref_int,alt_int,antimask,F1R2_count,F2R1_count =  genotypeDSSnv(readSet,ref_np[start_ind:end_ind],prior_mat[start_ind:end_ind,:],antimask,params)
                        ref_int = ref_np[start_ind:end_ind]
                        refs_ind = np.nonzero(np.logical_and(F1R2_ARLR <= -params["pcutoff"], F2R1_ARLR <= -params["pcutoff"]))[0].tolist()
                        muts_ind = np.nonzero(np.logical_and(F1R2_ARLR >= params["pcutoff"], F2R1_ARLR >= params["pcutoff"]))[0].tolist()
                        pass_bool = np.full(F1R2_ARLR.size,False,dtype=bool)
                        pass_bool[refs_ind] = True
                        pass_bool[muts_ind] = True
                        #pass_bool = np.logical_and(antimask,llt_pass)
                        coverage[start_ind:end_ind][pass_bool] += 1
                        #mut_bool = (DS_consensus != ref_np[start_ind:end_ind])
                        #muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                        pos = [mut_ind + start_ind + reference_mat_start for mut_ind in muts_ind]
                        #muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                        mut_positions = [mut_ind + start_ind + reference_mat_start + 1 for mut_ind in muts_ind]
                        if len(mut_positions) > 1:
                            break_flag = 0
                            for nn in range(1,len(mut_positions)):
                                if mut_positions[nn] - mut_positions[nn-1] != 1:
                                    break_flag = 1
                            if break_flag:
                                continue
                        NMs = [seq.get_tag("NM") for seq in readSet]
                        averageNM = sum(NMs)/len(NMs)
                        if averageNM - len(mut_positions) > 1: continue
                        for nn in range(len(mut_positions)):
                            mut_chrom = reference_mat_chrom
                            mut_pos = mut_positions[nn]
                            mut_ref = num2base[ref_int[muts_ind[nn]]]
                            mut_alt = num2base[alt_int[muts_ind[nn]]]
                            if muts_dict.get('_'.join([mut_chrom,str(mut_pos),mut_ref,mut_alt])): continue
                            na,nr,ndp = extractDepthSnv(normalBam,mut_chrom,mut_pos,mut_ref,mut_alt)
                            if na > params["maxAltCount"]: continue
                            ta,tr,tdp = extractDepthSnv(tumorBam,mut_chrom,mut_pos,mut_ref,mut_alt)
                            if readSet[0].is_forward:
                                readPos5p = min(muts_ind[nn] + 1,abs(readSet[0].template_length)-muts_ind[nn])
                                readPos3p = abs(readSet[0].reference_length)-muts_ind[nn]
                            else:
                                readPos5p = min(readSet[0].reference_length - muts_ind[nn], abs(readSet[0].template_length)-readSet[0].reference_length+muts_ind[nn])
                                readPos3p = readSet[0].reference_length - muts_ind[nn]
                            mut = {"chrom":mut_chrom, \
                                "pos":mut_pos, \
                                "ref":mut_ref, \
                                "alt":mut_alt, \
                                "infos":{"RP5":readPos5p,"RP3":readPos3p,"F1R2":F1R2,"F2R1":F2R1, \
                                "TG":F1R2_ARLR[muts_ind[nn]], \
                                "BG":F2R1_ARLR[muts_ind[nn]], \
                                "TC":",".join(F1R2_count[:,muts_ind[nn]].tolist()), \
                                "BC":",".join(F2R1_count[:,muts_ind[nn]].tolist())}, \
                                "formats":['AC',"RC","DP"],"samples":[[ta,tr,tdp],[na,nr,ndp]]}
                            muts_dict['_'.join([mut_chrom,str(mut_pos),mut_ref,mut_alt])] = 0
                            muts.append(mut) 
                    duplex_count += 1   
            """
            Calling block ends
            """
            currentReadDict = {bc:{"seqs":[rec],"F1R2":0,"F2R1":0}}
            if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
                currentReadDict[bc1+'+'+bc2]["F1R2"] += 1
            else:
                currentReadDict[bc1+'+'+bc2]["F2R1"] += 1    
            currentStart = start      

    for key in currentReadDict.keys():
        """
        Calling block starts
        """
        for key in currentReadDict.keys():
            readSet = currentReadDict[key]["seqs"]
            setBc = key.split('+')
            setBc1 = setBc[0]
            setBc2 = setBc[1]
            F2R1 = currentReadDict[key]["F2R1"]
            F1R2 = currentReadDict[key]["F1R2"]
            readSet = currentReadDict[key]["seqs"]
            setBc = key.split('+')
            setBc1 = setBc[0]
            setBc2 = setBc[1]
            F2R1 = currentReadDict[key]["F2R1"]
            F1R2 = currentReadDict[key]["F1R2"]
            duplex_no = f"{min([F1R2,F2R1])}+{max([F1R2,F2R1])}"
            if duplex_read_num_dict.get(duplex_no):
                duplex_read_num_dict[duplex_no] += 1
            else:
                duplex_read_num_dict[duplex_no] = 1
            if setBc1 != setBc2 and F2R1 > 1 and F1R2 > 1:
                #if chromPast == 'startChrom':
                    #chromPast = tumorBam.get_reference_name(rec.reference_id)
                #else:
                    #chromPast = chrom
                #chrom = tumorBam.get_reference_name(rec.reference_id)
                ### Determine if there is any indels in reads
                indel_bool = [('I' in seq.cigarstring or 'D' in seq.cigarstring) for seq in readSet]
                if any(indel_bool): 
                    mlmlml=1
                    #genotypeDSIndel(readSet,tumorBam,params)

                





                else:
                    rs_reference_end = max([r.reference_end for r in readSet])
                    chromNow = readSet[0].reference_name
                    if chromNow != reference_mat_chrom or rs_reference_end >= reference_mat_end:
                        ### Output coverage 
                        if "coverage" in locals():
                            non_zero_positions = np.nonzero(coverage)
                            for pos in non_zero_positions[0].tolist():
                                locus_bed.write(('\t'.join([reference_mat_chrom,str(pos),str(pos+1),str(coverage[pos])])+"\n").encode("utf-8"))
                            total_coverage += np.sum(coverage)
                        reference_mat_chrom = chromNow
                        reference_mat_start = readSet[0].reference_start
                        try: region_end = region[2]
                        except: region_end = 10E10
                        contig_end = tumorBam.get_reference_length(chromNow)
                        reference_mat_end = min(readSet[0].reference_start + 100000,max(region_end,readSet[0].reference_end),contig_end)
                        try:prior_mat,snp_mask,indel_mask,noise_mask,n_cov_mask,ref_np = prepare_reference_mats(reference_mat_chrom,\
                                                                                            reference_mat_start,\
                                                                                            reference_mat_end,fasta[reference_mat_chrom][reference_mat_start:reference_mat_end],\
                                                                                            germline,\
                                                                                            noise,\
                                                                                            ncoverage,\
                                                                                            params)
                        except: print(readSet[0].reference_start + 100000,region_end,readSet[0].reference_end,contig_end,reference_mat_start,chromNow,rec.reference_name,readSet[0].reference_name)
                        coverage = np.zeros(100000,dtype=int)
                    ### Record read names to check if mate has been processed
                    processed_flag = 0
                    for seq in readSet:
                        if seq.query_name in processed_read_names:
                            processed_read_names.remove(seq.query_name)
                            processed_flag = 1
                            break
                    if processed_flag == 0:
                        processed_read_names.add(readSet[0].query_name)

                    start_ind = readSet[0].reference_start - reference_mat_start
                    reference_length = min([read.reference_length for read in readSet])
                    end_ind = readSet[0].reference_start + reference_length - reference_mat_start
                    masks = np.zeros([4,end_ind-start_ind],dtype=bool)
                    try:masks[0,:] = snp_mask[start_ind:end_ind]
                    except: print(start_ind,end_ind)
                    masks[1,:] = noise_mask[start_ind:end_ind]
                    masks[2,:] = n_cov_mask[start_ind:end_ind]
                    left,right = determineTrimLength(readSet[0],params=params,processed_flag=processed_flag)
                    masks[3,:left] = True
                    masks[3,-right:] = True 
                    antimask = np.all(~masks,axis=0)
                    ### If the whole reads are masked:
                    if not np.any(antimask): continue
                    ### Calculate genotype probability
                    F1R2_ARLR,F2R1_ARLR,ref_int,alt_int,antimask,F1R2_count,F2R1_count =  genotypeDSSnv(readSet,ref_np[start_ind:end_ind],prior_mat[start_ind:end_ind,:],antimask,params)
                    refs_ind = np.nonzero(np.logical_and(F1R2_ARLR <= -params["pcutoff"], F2R1_ARLR <= -params["pcutoff"]))[0].tolist()
                    muts_ind = np.nonzero(np.logical_and(F1R2_ARLR >= params["pcutoff"], F2R1_ARLR >= params["pcutoff"]))[0].tolist()
                    pass_bool = np.full(F1R2_ARLR.size,False,dtype=bool)
                    pass_bool[refs_ind] = True
                    pass_bool[muts_ind] = True
                    #pass_bool = np.logical_and(antimask,llt_pass)
                    coverage[start_ind:end_ind][pass_bool] += 1
                    #mut_bool = (DS_consensus != ref_np[start_ind:end_ind])
                    #muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                    pos = [mut_ind + start_ind + reference_mat_start for mut_ind in muts_ind]
                    #muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                    mut_positions = [mut_ind + start_ind + reference_mat_start + 1 for mut_ind in muts_ind]
                    if len(mut_positions) > 1:
                        break_flag = 0
                        for nn in range(1,len(mut_positions)):
                            if mut_positions[nn] - mut_positions[nn-1] != 1:
                                break_flag = 1
                        if break_flag:
                            break
                    for nn in range(len(mut_positions)):
                        mut_chrom = reference_mat_chrom
                        mut_pos = mut_positions[nn]
                        mut_ref = num2base[ref_int[muts_ind[nn]]]
                        mut_alt = num2base[alt_int[muts_ind[nn]]]
                        if muts_dict.get('_'.join([mut_chrom,str(mut_pos),mut_ref,mut_alt])): continue
                        na,nr,ndp = extractDepthSnv(normalBam,mut_chrom,mut_pos,mut_ref,mut_alt)
                        if na > params["maxAltCount"]: continue
                        ta,tr,tdp = extractDepthSnv(tumorBam,mut_chrom,mut_pos,mut_ref,mut_alt)
                        if readSet[0].is_forward:
                            readPos5p = min(muts_ind[nn] + 1,abs(readSet[0].template_length)-muts_ind[nn])
                            readPos3p = abs(readSet[0].reference_length)-muts_ind[nn]
                        else:
                            readPos5p = min(readSet[0].reference_length - muts_ind[nn], abs(readSet[0].template_length)-readSet[0].reference_length+muts_ind[nn])
                            readPos3p = readSet[0].reference_length - muts_ind[nn]
                        mut = {"chrom":mut_chrom, \
                            "pos":mut_pos, \
                            "ref":mut_ref, \
                            "alt":mut_alt, \
                            "infos":{"RP5":readPos5p,"RP3":readPos3p,"F1R2":F1R2,"F2R1":F2R1, \
                            "TG":F1R2_ARLR[muts_ind[nn]], \
                            "BG":F2R1_ARLR[muts_ind[nn]], \
                            "TC":",".join(F1R2_count[:,muts_ind[nn]].tolist()), \
                            "BC":",".join(F2R1_count[:,muts_ind[nn]].tolist())}, \
                            "formats":['AC',"RC","DP"],"samples":[[ta,tr,tdp],[na,nr,ndp]]}
                        muts_dict['_'.join([mut_chrom,str(mut_pos),mut_ref,mut_alt])] = 0
                        muts.append(mut) 
                duplex_count += 1   
        """
        Calling block ends
        """
    return muts,total_coverage,recCount,duplex_count,duplex_read_num_dict

