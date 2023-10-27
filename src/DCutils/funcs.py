from subprocess import check_output
import subprocess
from scipy.stats import binom
from Bio import SeqIO
import pysam
import numpy as np
import re
from pysam import AlignmentFile as BAM
from pysam import VariantFile as VCF
from pysam import index as indexBam
from pysam import TabixFile as BED
import gzip


"""
def extractDepthSnv(bam,chrom,pos,ref,alt):
    allele = []
    for pileupcolumn in bam.pileup(chrom,pos-1,pos,ignore_orphans=False,min_base_quality=0):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                #if pileupread.
                if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:
                    allele.append(pileupread.alignment.query_sequence[pileupread.query_position])
            break
    alt_count = allele.count(alt)
    ref_count = allele.count(ref)
    alts_count = len(allele) - ref_count
    depth = len(allele)
    return alt_count,alts_count,ref_count,depth
"""


"""
def extractDepthSnv(bam,chrom,pos,ref,alt):
    pileupLine = check_output(["samtools","mpileup","-Q","1",bam,"-r",chrom+":"+str(pos)+"-"+str(pos)],stderr=subprocess.DEVNULL).decode('ascii')
    infos = pileupLine.strip('\n').split('\t')
    if len(infos) <= 1: return 0,0,0
    depth = int(infos[3])
    alleles = infos[4].upper()
    altCount = alleles.count(alt)
    refCount = alleles.count(ref)
    return altCount,refCount,depth

"""


def findIndels(seq):
    refPos = seq.reference_start
    readPos = 0
    indels = {}
    for cigar in seq.cigartuples:
        if cigar[0] == 0:
            refPos += cigar[1]
            readPos += cigar[1]
        if cigar[0] in [3, 4]:
            readPos += cigar[1]
        if cigar[0] == 1:
            pos = refPos
            sequence = seq.query_sequence[readPos : readPos + cigar[1]]
            quals = seq.query_qualities[readPos : readPos + cigar[1]]
            indels.update({str(pos) + "I": [sequence, quals, readPos]})
            readPos += cigar[1]
        if cigar[0] == 2:
            pos = refPos
            indels.update({str(pos) + "D": [cigar[1], readPos]})
            refPos += cigar[1]
    return indels


def extractDepthSnv(bam, chrom, pos, ref, alt):
    has_barcode = True
    for seq in bam.fetch(chrom, pos - 1, pos):
        if len(seq.query_name.split("_")[-1].split("+")) != 2:
            has_barcode = False
        break
    if has_barcode:
        refReadBarcodes = {}
        altReadBarcodes = {}
        for pileupcolumn in bam.pileup(chrom, pos - 1, pos, min_base_quality=0):
            if pileupcolumn.pos == pos - 1:
                for pileupread in pileupcolumn.pileups:
                    if (
                        pileupread.is_del
                        or pileupread.is_refskip
                        or pileupread.alignment.is_secondary
                        or pileupread.alignment.is_supplementary
                    ):
                        continue
                    if (
                        pileupread.alignment.query_sequence[pileupread.query_position]
                        == alt
                    ):
                        barcodes = pileupread.alignment.query_name.split("_")[-1].split(
                            "+"
                        )
                        if not altReadBarcodes.get(
                            barcodes[0] + "+" + barcodes[1]
                        ) or not altReadBarcodes.get(barcodes[1] + "+" + barcodes[0]):
                            altReadBarcodes.update({barcodes[0] + "+" + barcodes[1]: 1})
                    else:
                        barcodes = pileupread.alignment.query_name.split("_")[-1].split(
                            "+"
                        )
                        if not refReadBarcodes.get(
                            barcodes[0] + "+" + barcodes[1]
                        ) or not refReadBarcodes.get(barcodes[1] + "+" + barcodes[0]):
                            refReadBarcodes.update({barcodes[0] + "+" + barcodes[1]: 1})
        altAlleleCount = len(altReadBarcodes.keys())
        refAlleleCount = len(refReadBarcodes.keys())
        depth = altAlleleCount + refAlleleCount
    else:
        altAlleleCount = 0
        refAlleleCount = 0
        for pileupcolumn in bam.pileup(chrom, pos - 1, pos, min_base_quality=0):
            if pileupcolumn.pos == pos - 1:
                for pileupread in pileupcolumn.pileups:
                    if (
                        pileupread.is_del
                        or pileupread.is_refskip
                        or pileupread.alignment.is_secondary
                        or pileupread.alignment.is_supplementary
                        or pileupread.alignment.is_duplicate
                    ):
                        continue
                    if (
                        pileupread.alignment.query_sequence[pileupread.query_position]
                        == alt
                    ):
                        altAlleleCount += 1
                    else:
                        refAlleleCount += 1
        depth = altAlleleCount + refAlleleCount
    return altAlleleCount, refAlleleCount, depth


def extractDepthIndel(bam, chrom, pos, ref, alt):
    if len(ref) == 1:
        indelPos = str(pos) + "I"
    else:
        indelPos = str(pos) + "D"
    refReadBarcodes = {}
    altReadBarcodes = {}
    for read in bam.fetch(chrom, pos - 1, pos):
        matchFlag = 0
        if read.mapping_quality <= 30:
            continue
        if "I" in read.cigarstring or "D" in read.cigarstring:
            indels = findIndels(read)
            if indels.get(indelPos):
                matchFlag = 1
        if matchFlag == 1:
            barcodes = read.query_name.split("_")[-1].split("+")
            if not altReadBarcodes.get(
                barcodes[0] + "+" + barcodes[1]
            ) or not altReadBarcodes.get(barcodes[1] + "+" + barcodes[0]):
                altReadBarcodes.update({barcodes[0] + "+" + barcodes[1]: 1})
        else:
            barcodes = read.query_name.split("_")[-1].split("+")
            if not refReadBarcodes.get(
                barcodes[0] + "+" + barcodes[1]
            ) or not refReadBarcodes.get(barcodes[1] + "+" + barcodes[0]):
                refReadBarcodes.update({barcodes[0] + "+" + barcodes[1]: 1})
    altAlleleCount = len(altReadBarcodes.keys())
    refAlleleCount = len(refReadBarcodes.keys())
    depth = altAlleleCount + refAlleleCount
    return altAlleleCount, refAlleleCount, depth


def createVcfStrings(chromDict, infoDict, formatDict, filterDict, recs):
    lines = ["##fileformat=VCFv4.2"]
    for filterr in filterDict.keys():
        lines.append(
            '##FILTER=<ID={},Description="{}">'.format(filterr, filterDict[filterr])
        )
    for info in infoDict.keys():
        lines.append(
            '##INFO=<ID={},Number={},Type={},Description="{}">'.format(
                info, infoDict[info][0], infoDict[info][1], infoDict[info][2]
            )
        )
    for form in formatDict.keys():
        lines.append(
            '##FORMAT=<ID={},Number={},Type={},Description="{}">'.format(
                form, formatDict[form][0], formatDict[form][1], formatDict[form][2]
            )
        )
    for chrom in chromDict.keys():
        lines.append("##contig=<ID={},length={}>".format(chrom, chromDict[chrom]))
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL")
    for rec in recs:
        chrom = rec["chrom"]
        pos = rec["pos"]
        alt = rec["alt"]
        ref = rec["ref"]
        infos = rec["infos"]
        formats = rec["formats"]
        samples = rec["samples"]
        lineEntries = [
            chrom,
            str(pos),
            ".",
            ref,
            alt,
            ".",
            "PASS",
            ";".join(
                [
                    "{info}={val}".format(info=info, val=infos[info])
                    for info in infoDict.keys()
                ]
            ),
            ":".join(formats),
            ":".join([str(s) for s in samples[0]]),
            ":".join([str(s) for s in samples[1]]),
        ]
        lines.append("\t".join(lineEntries))
    return "\n".join(lines) + "\n"


def checkCallablity(F1R2, F2R1, F1R2BQ, F2R1BQ, mutRate, perror, pcut):
    mockFBaseProb = genotypeProbSnv(
        ["A"] * F1R2, F1R2BQ, [mutRate, mutRate, mutRate, 1 - 3 * mutRate], perror
    )
    mockRBaseProb = genotypeProbSnv(
        ["A"] * F2R1, F2R1BQ, [mutRate, mutRate, mutRate, 1 - 3 * mutRate], perror
    )
    if mockFBaseProb[3] <= pcut and mockRBaseProb[3] <= pcut:
        return True
    return False

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
    return np.log10(np.where(mat > 0, mat, np.finfo(float).eps))


def power10(mat):
    return 10 ** np.where(mat >= np.log10(np.finfo(float).eps), mat, -np.inf)


def log2(mat):
    return np.log2(np.where(mat > 0, mat, np.finfo(float).eps))


def power2(mat):
    return np.where(mat >= -100, 2**mat, 0)


def log(mat):
    return np.log(np.where(mat > 0, mat, np.finfo(float).eps))

def calculatePosterior(Pamp, Pref, Palt, prior_ref, prior_alt):
    log_e_10 = np.log(10)
    Pamp = power10(Pamp)
    Pref = power10(Pref)
    Palt = power10(Palt)

    first_term = log10(1 + 2 * Pamp * Pref - Pamp - Pref)
    second_term = log10(Pamp + Palt - Pamp * Palt)
    ref_prob = prior_ref + first_term + second_term

    first_term = log10(1 + 2 * Pamp * Palt - Pamp - Palt)
    second_term = log10(Pamp + Pref - Pamp * Pref)
    alt_prob = prior_alt + first_term + second_term
    return ref_prob, alt_prob


def prepare_reference_mats(
    chrom, start, end, reference_seq, germline_bed, noise_bed, n_cov_bed, params
):
    ### Define and Initialize
    af_miss = params["mutRate"]
    m = end - start
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    snp_mask = np.full(m, False, dtype=bool)
    indel_mask = np.full(m, False, dtype=bool)
    noise_mask = np.full(m, False, dtype=bool)
    n_cov_mask = np.full(m, False, dtype=bool)

    ### Initialize prior mat as no germline resource
    reference_int = np.array([base2num.get(b, 4) for b in reference_seq], dtype=int)
    noise_mask[reference_int == 4] = True
    reference_int[reference_int == 4] = 0
    prior_mat = np.full([m, 4], af_miss, dtype=float)
    prior_mat[np.ogrid[:m], reference_int] = 1 - 3 * af_miss

    ### Adjust by germline
    for rec in germline_bed.fetch(chrom, start, end):
        ind = rec.pos - 1 - start
        ref = rec.ref
        afs = rec.info["AF"]
        if len(ref) == 1:
            has_snp = False
            for ii, alt in enumerate(rec.alts):
                if len(alt) == 1:
                    prior_mat[ind, base2num[alt]] = afs[ii]
                    has_snp = True
                    if afs[ii] >= 0.01:
                        snp_mask[ind] = True
                elif afs[ii] >= 0.01:
                    indel_mask[ind] = True
            if has_snp:
                prior_mat[ind, base2num[ref]] = (
                    1 - np.delete(prior_mat[ind, :], base2num[ref]).sum()
                )
        if len(ref) != 1:
            for ii, alt in enumerate(rec.alts):
                if afs[ii] >= 0.01:
                    indel_mask[ind + 1 : ind + len(alt)] = True

    ### Preparep noise mask
    for rec in noise_bed.fetch(chrom, start, end, parser=pysam.asBed()):
        interval_start = max(rec.start, start)
        interval_end = min(rec.end, end)
        interval_len = interval_end - interval_start
        interval_start_ind = interval_start - 1 - start
        noise_mask[interval_start_ind : interval_start_ind + interval_len] = True

    ### Preparep normal coverage mask
    for rec in n_cov_bed.fetch(chrom, start, end, parser=pysam.asTuple()):
        depth = int(rec[3])
        if depth >= params["minNdepth"]:
            continue
        interval_start = max(int(rec[1]), start)
        interval_end = min(int(rec[2]), end)
        interval_len = interval_end - interval_start
        interval_start_ind = interval_start - 1 - start
        n_cov_mask[interval_start_ind : interval_start_ind + interval_len] = True
    return prior_mat, snp_mask, indel_mask, noise_mask, n_cov_mask, reference_int


def determineTrimLength(seq, params, processed_flag):
    if seq.is_forward:
        overlap = 0  # Never mask overlap of forward read
        left = params["trim5"]
        right_frag = params["trim5"] - min(
            params["trim5"], abs(seq.template_length) - seq.reference_length
        )
        right_read = params["trim3"]
        right = max(right_frag, right_read)
    else:
        ### Mask overlap of reverse read
        if processed_flag:
            mate_cigar = seq.get_tag("MC")
            re_cigar = re.search("(?:(\d+)S)?(\d+)M(?:(\d+)S)?", mate_cigar)
            mate_reference_length = int(re_cigar.group(2))
            overlap = max(
                0,
                seq.reference_length + mate_reference_length - abs(seq.template_length),
            )
        else:
            overlap = 0
        right_frag = params["trim5"]
        right_read = params["trim3"]
        right = max(right_frag, right_read)
        left_frag = params["trim5"] - min(
            params["trim5"], abs(seq.template_length) - seq.reference_length
        )
        left = max(left_frag, overlap, params["trim3"])
    return left, right


def genotypeDSSnv(seqs, reference_int, prior_mat, antimask, params):
    """debug
    params["print_flag"] = 0
    np.set_printoptions(
        edgeitems=30,
        linewidth=10000,
        suppress=True,
    )
    """
    #prob_damage = params["amperr"]
    #prob_damage_mat = np.ones([4, 4], dtype=float) * prob_damage / 3
    #np.fill_diagonal(prob_damage_mat, 1 - prob_damage)
    #log_prob_damage_mat = log10(prob_damage_mat)
    prob_amp = params["amperr"]
    prob_amp_mat = np.ones([4, 4], dtype=float) * prob_amp / 3
    np.fill_diagonal(prob_amp_mat, 1 - prob_amp)
    log_prob_amp_mat = log10(prob_amp_mat)

    # seq_mat = 1-power10(-qual_mat_merged/10).reshape(4,n,1)

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
    re_cigar = re.search("(?:(\d+)S)?(\d+)M(?:(\d+)S)?", cigar)
    leftS = re_cigar.group(1) or 0
    matches = re_cigar.group(2)
    rightS = re_cigar.group(3) or 0

    n = min([seq.reference_length for seq in seqs])
    m_F1R2 = len(F1R2)
    m_F2R1 = len(F2R1)

    ### Prepare sequence matrix and quality matrix for each strand
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    F1R2_seq_mat = np.zeros([m_F1R2, n], dtype=int)  # Base(ATCG) x reads x pos
    F1R2_qual_mat = np.zeros([m_F1R2, n])
    F2R1_seq_mat = np.zeros([m_F2R1, n], dtype=int)  # Base(ATCG) x reads x pos
    F2R1_qual_mat = np.zeros([m_F2R1, n])

    for mm, seq in enumerate(F1R2):
        qualities = seq.query_alignment_qualities
        sequence = seq.query_alignment_sequence
        for nn in range(n):
            base = sequence[nn]
            if base == "N":
                F1R2_seq_mat[mm, nn] = 0
                F1R2_qual_mat[mm, nn] = np.log10(0.25)
            else:
                F1R2_seq_mat[mm, nn] = base2num[base]
                F1R2_qual_mat[mm, nn] = qualities[nn]

    F1R2_qual_mat_merged = np.zeros([4, n])
    F1R2_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F1R2_qual_mat_merged[nn, :] = F1R2_qual_mat.sum(
            axis=0, where=(F1R2_seq_mat == nn)
        )
        F1R2_count_mat[nn, :] = (F1R2_seq_mat == nn).sum(axis=0)

    for mm, seq in enumerate(F2R1):
        qualities = seq.query_alignment_qualities
        sequence = seq.query_alignment_sequence
        for nn in range(n):
            base = sequence[nn]
            if base == "N":
                F2R1_seq_mat[mm, nn] = 0
                F2R1_qual_mat[mm, nn] = np.log10(0.25)
            else:
                F2R1_seq_mat[mm, nn] = base2num[base]
                F2R1_qual_mat[mm, nn] = qualities[nn]
    F2R1_qual_mat_merged = np.zeros([4, n])
    F2R1_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F2R1_qual_mat_merged[nn, :] = F2R1_qual_mat.sum(
            axis=0, where=(F2R1_seq_mat == nn)
        )
        F2R1_count_mat[nn, :] = (F2R1_seq_mat == nn).sum(axis=0)

    total_count_mat = F1R2_count_mat + F2R1_count_mat
    total_count_mat[reference_int, np.ogrid[:n]] = -1
    antimask[(total_count_mat >= 1).sum(axis=0) > 1] = False
    alt_int = np.argmax(total_count_mat, axis=0)

    F1R2_masked_qual_mat = F1R2_qual_mat_merged[:, antimask]
    F2R1_masked_qual_mat = F2R1_qual_mat_merged[:, antimask]
    reference_int_masked = reference_int[antimask]
    alt_int_masked = alt_int[antimask]
    Pamp = log_prob_amp_mat[reference_int_masked, alt_int_masked]

    mmm = alt_int_masked.size
    prior_ref = log10(prior_mat[antimask][np.ogrid[:mmm], reference_int_masked])
    prior_alt = log10(prior_mat[antimask][np.ogrid[:mmm], alt_int_masked])

    F1R2_ref_prob_mat = (
        -F1R2_masked_qual_mat[
            reference_int_masked, np.ogrid[: reference_int_masked.size]
        ]
        / 10
    )
    F1R2_alt_prob_mat = (
        -F1R2_masked_qual_mat[alt_int_masked, np.ogrid[: reference_int_masked.size]]
        / 10
    )
    F2R1_ref_prob_mat = (
        -F2R1_masked_qual_mat[
            reference_int_masked, np.ogrid[: reference_int_masked.size]
        ]
        / 10
    )
    F2R1_alt_prob_mat = (
        -F2R1_masked_qual_mat[alt_int_masked, np.ogrid[: reference_int_masked.size]]
        / 10
    )

    # print(prior_mat)
    F1R2_ref_prob, F1R2_alt_prob = calculatePosterior(
        Pamp, F1R2_ref_prob_mat, F1R2_alt_prob_mat, prior_ref, prior_alt
    )
    F2R1_ref_prob, F2R1_alt_prob = calculatePosterior(
        Pamp, F2R1_ref_prob_mat, F2R1_alt_prob_mat, prior_ref, prior_alt
    )

    # TS_gt = np.ones([4,n]) * 0.25
    # BS_gt = np.ones([4,n]) * 0.25
    F1R2_ARLR_masked = F1R2_alt_prob - F1R2_ref_prob
    F2R1_ARLR_masked = F2R1_alt_prob - F2R1_ref_prob

    F1R2_ARLR = np.zeros(n)
    F2R1_ARLR = np.zeros(n)
    F1R2_ARLR[antimask] = F1R2_ARLR_masked
    F2R1_ARLR[antimask] = F2R1_ARLR_masked
    # print(F1R2_ref_qual_mat)
    # print(F1R2_alt_qual_mat)
    return (
        F1R2_ARLR,
        F2R1_ARLR,
        reference_int,
        alt_int,
        antimask,
        F1R2_count_mat,
        F2R1_count_mat,
    )

def genotypeDSIndel(seqs, bam, params):
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
    f1r2_names_dict = {f1r2.query_name: nn for nn, f1r2 in enumerate(F1R2)}
    f2r1_names_dict = {f2r1.query_name: nn for nn, f2r1 in enumerate(F2R1)}
    f1r2_genotype = np.zeros([len(F1R2), 1])
    f2r1_genotype = np.zeros([len(F2R1), 1])
    for nn, col in enumerate(bam.pileup(chrom, start, end, flag_filter=2816)):
        indel_flag = 0
        for read in col.pileups:
            if read.alignment.query_name in f1r2_names_dict:
                if read.indel != 0:
                    indel_flag = 1
                    f1r2_genotype[
                        f1r2_names_dict[read.alignment.query_name], 0
                    ] = read.indel
                print(read.indel, read.alignment.query_name, f1r2_names_dict)
        if indel_flag == 1:
            break
    for nn, col in enumerate(bam.pileup(chrom, start, end, flag_filter=2816)):
        indel_flag = 0
        for read in col.pileups:
            if read.alignment.query_name in f1r2_names_dict:
                if read.indel != 0:
                    indel_flag = 1
                    f1r2_genotype[
                        f1r2_names_dict[read.alignment.query_name], 0
                    ] = read.indel
                print(read.indel, read.alignment.query_name, f1r2_names_dict)
        if indel_flag == 1:
            break

    log_gt_mat = np.ones([n, 4]) * np.log(0.25)
    F1R2_gt = genotypeSSSnv(
        F1R2_seq_mat[:, antimask],
        F1R2_qual_mat[:, antimask],
        prior_mat[antimask, :],
        params,
    )
    F2R1_gt = genotypeSSSnv(
        F2R1_seq_mat[:, antimask],
        F2R1_qual_mat[:, antimask],
        prior_mat[antimask, :],
        params,
    )
    m, n = F1R2_gt.shape
    # pairwise = np.repeat(F1R2_gt.reshape((m,n,1)),4,axis=2) + np.repeat(F2R1_gt.reshape((m,1,n)),4,axis=1)
    log_gt_mat_masked = np.zeros([m, 4])
    consensus_gt = np.logical_and(F1R2_gt >= 0.95, F2R1_gt >= 0.95, dtype=int)

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
    # gt_mat_masked = power10(log_gt_mat_masked)
    # log_gt_mat_masked = log10(gt_mat_masked / gt_mat_masked.sum(axis=1,keepdims=True))
    gt_mat[antimask, :] = consensus_gt
    return gt_mat

def profileTriNucMismatches(seqs,reference_int,chrom,start,end,params):
    fasta = params["reference"]
    reverse_comp = {"A":"T",'T':"A","C":"G","G":"C"}
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2base = "ATCG"
    base_changes = ['C>A','C>G','C>T','T>A','T>C','T>G']

    F1R2 = []
    F2R1 = []
    for seq in seqs:
        if (seq.is_read1 and seq.is_forward) or (seq.is_read2 and seq.is_reverse):
            F1R2.append(seq)
        if (seq.is_read2 and seq.is_forward) or (seq.is_read1 and seq.is_reverse):
            F2R1.append(seq)

    ### Determine match length

    cigar = seqs[0].cigarstring
    re_cigar = re.search("(?:(\d+)S)?(\d+)M(?:(\d+)S)?", cigar)
    leftS = re_cigar.group(1) or 0
    matches = re_cigar.group(2)
    rightS = re_cigar.group(3) or 0

    n = min([seq.reference_length for seq in seqs])
    F1R2_antimask = np.full(nn,True)
    F2R1_antimask = np.full(nn,True)
    m_F1R2 = len(F1R2)
    m_F2R1 = len(F2R1)

    ### Prepare sequence matrix and quality matrix for each strand
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    F1R2_seq_mat = np.zeros([m_F1R2, n], dtype=int)  # Base(ATCG) x reads x pos
    F1R2_qual_mat = np.zeros([m_F1R2, n])
    F2R1_seq_mat = np.zeros([m_F2R1, n], dtype=int)  # Base(ATCG) x reads x pos
    F2R1_qual_mat = np.zeros([m_F2R1, n])

    for mm, seq in enumerate(F1R2):
        qualities = seq.query_alignment_qualities
        sequence = seq.query_alignment_sequence
        for nn in range(n):
            base = sequence[nn]
            if base == "N":
                F1R2_seq_mat[mm, nn] = 0
                F1R2_qual_mat[mm, nn] = np.log10(0.25)
            else:
                F1R2_seq_mat[mm, nn] = base2num[base]
                F1R2_qual_mat[mm, nn] = qualities[nn]

    F1R2_qual_mat_merged = np.zeros([4, n])
    F1R2_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F1R2_qual_mat_merged[nn, :] = F1R2_qual_mat.sum(
            axis=0, where=(F1R2_seq_mat == nn)
        )
        F1R2_count_mat[nn, :] = (F1R2_seq_mat == nn).sum(axis=0)

    for mm, seq in enumerate(F2R1):
        qualities = seq.query_alignment_qualities
        sequence = seq.query_alignment_sequence
        for nn in range(n):
            base = sequence[nn]
            if base == "N":
                F2R1_seq_mat[mm, nn] = 0
                F2R1_qual_mat[mm, nn] = np.log10(0.25)
            else:
                F2R1_seq_mat[mm, nn] = base2num[base]
                F2R1_qual_mat[mm, nn] = qualities[nn]
    F2R1_qual_mat_merged = np.zeros([4, n])
    F2R1_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F2R1_qual_mat_merged[nn, :] = F2R1_qual_mat.sum(
            axis=0, where=(F2R1_seq_mat == nn)
        )
        F2R1_count_mat[nn, :] = (F2R1_seq_mat == nn).sum(axis=0)

    F1R2_count_mat[reference_int, np.ogrid[:n]] = -1
    F1R2_alt_int = np.argmax(F1R2_count_mat, axis=0)
    F2R1_count_mat[reference_int, np.ogrid[:n]] = -1
    F2R1_alt_int = np.argmax(F1R2_count_mat, axis=0)
    
    #total_qual_mat = F1R2_qual_mat_merged + F2R1_qual_mat_merged
    F1R2_ref_qual = F1R2_qual_mat_merged[reference_int, np.ogrid[:n]] 
    F1R2_alt_qual = F1R2_qual_mat_merged[F1R2_alt_int, np.ogrid[:n]] 
    F1R2_antimask[(F1R2_count_mat >= 1).sum(axis=0) > 1] = False
    F1R2_antimask[F1R2_ref_qual<=60] = False
    F1R2_antimask[np.logical_and(F1R2_alt_qual<=60,F1R2_alt_qual!=0)]  = False
    F1R2_antimask[ref_int == 4] = False
    #alt_ind = np.where(alt_qual!=0)[0]
    
    F2R1_ref_qual = F2R1_qual_mat_merged[reference_int, np.ogrid[:n]] 
    F2R1_alt_qual = F2R1_qual_mat_merged[F2R1_alt_int, np.ogrid[:n]] 
    F2R1_antimask[(F2R1_count_mat >= 1).sum(axis=0) > 1] = False
    F2R1_antimask[F2R1_ref_qual<=60] = False
    F2R1_antimask[np.logical_and(F2R1_alt_qual<=60,F2R1_alt_qual!=0)]  = False
    F2R1_antimask[ref_int == 4] = False

    mismatch_dict = dict()
    for minus_base in ['A','T','C','G']:
        for ref_base in ['C','T']:
                for plus_base in ['A','T','C','G']:
                    mismatch_dict[minus_base+ref_base + plus_base] = [0,0,0,0]
    for nn in range(n):
        if F1R2_antimask[nn]:
            ref_base = num2base[reference_int[nn]]
            alt_base = num2base[F1R2_alt_int[nn]]
            minus_base = fasta[chrom][start+nn-1]
            plus_base = fasta[chrom][start+nn+1]
            if ref_base == 'C' or ref_base == 'T':
                mismatch_dict[minus_base + ref_base + plus_base][F1R2_alt_int[nn]] += 1
            else:
                mismatch[reverse_comp[plus_base] + reverse_comp[ref_base] + reverse_comp[minus_base]][base2num[reverse_comp[alt_base]]] += 1
        if F2R1_antimask[nn]:
            ref_base = num2base[reference_int[nn]]
            alt_base = num2base[F2R1_alt_int[nn]]
            minus_base = fasta[chrom][start+nn-1]
            plus_base = fasta[chrom][start+nn+1]
            if ref_base == 'C' or ref_base == 'T':
                mismatch_dict[minus_base + ref_base + plus_base][F2R1_alt_int[nn]] += 1
            else:
                mismatch[reverse_comp[plus_base] + reverse_comp[ref_base] + reverse_comp[minus_base]][base2num[reverse_comp[alt_base]]] += 1
    print(mismatch)
    return mismatch
    

    

    

    
