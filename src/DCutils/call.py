#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO
import pysam
import numpy as np
import re
from pysam import AlignmentFile as BAM
from pysam import VariantFile as VCF
from pysam import index as indexBam
from pysam import TabixFile as BED
import gzip

# from pysam import FastaFile as FASTA
import time
from DCutils.funcs import *



def bamIterateMultipleRegion(bamObject, regions):
    for region in regions:
        for rec in bamObject.fetch(*region):
            if len(region) >= 2:
                if rec.reference_start < region[1]:
                    continue
            yield rec, region


def callBam(params, processNo, chunkSize):
    ### Get parameters
    bam = params["tumorBam"]
    nbam = params["normalBam"]
    regions = params["regions"]
    germline = VCF(params["germline"], "rb")
    ncoverage = BED(params["ncoverage"])
    reference = params["reference"]
    minMapq = params["mapq"]
    mutRate = params["mutRate"]
    pcut = params["pcutoff"]
    isLearn = params.get("isLearn",False)
    nn = processNo
    output = "tmp/" + params["output"] + "_" + str(nn)
    noise = BED(params["noise"])
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2base = "ATCG"
    muts = []
    muts_dict = dict()
    duplex_read_num_dict = dict()
    FPs = []
    RPs = []
    mismatch_dict = dict()
    for minus_base in ['A','T','C','G']:
        for ref_base in ['C','T']:
                for plus_base in ['A','T','C','G']:
                    mismatch_dict[minus_base+ref_base + plus_base] = [0,0,0,0]

    print("Process" + str(processNo) + ": Initiated")
    ### Initialize
    total_coverage = 0
    starttime = time.time()
    tumorBam = BAM(bam, "rb")
    normalBam = BAM(nbam, "rb")
    currentReadSet = []
    currentBc = ""
    currentStart = -1
    currentReadDict = {}
    # chromPast = 'startChrom'
    fasta = reference
    recCount = 0
    currentCheckPoint = 100000
    duplex_count = 0
    reference_mat_chrom = "anyChrom"
    reference_mat_start = 0
    locus_bed = gzip.open(output + "_coverage.bed.gz", "wb")
    processed_read_names = set()

    for rec, region in bamIterateMultipleRegion(tumorBam, regions):
        recCount += 1
        if recCount == currentCheckPoint:
            print(
                "Process"
                + str(processNo)
                + ": processed "
                + str(recCount)
                + " reads in "
                + str((time.time() - starttime) / 60)
                + " minutes"
            )
            currentCheckPoint += 100000
        if (
            rec.mapping_quality <= minMapq
            or rec.is_supplementary
            or rec.is_secondary
            or rec.has_tag("DT")
            or not rec.is_proper_pair
            or rec.is_qcfail
        ):
            continue
        if rec.get_tag("AS") - rec.get_tag("XS") <= 50:
            continue
        # if rec.cigartuples[0][0] == 4: continue
        #if rec.cigarstring.count("I") + rec.cigarstring.count("D") >= 2:
            #continue
        start = rec.reference_start
        bc = rec.query_name.split("_")[1]
        bcsplit = bc.split("+")
        bc1 = bcsplit[0]
        bc2 = bcsplit[1]
        chrom = tumorBam.get_reference_name(rec.reference_id)
        if currentStart == -1:
            currentStart = start
        # print(start,currentStart)
        if start == currentStart:
            # print(currentReadDict,1)
            if currentReadDict.get(bc1 + "+" + bc2) != None:
                currentReadDict[bc1 + "+" + bc2]["seqs"].append(rec)
                if (rec.is_forward and rec.is_read1) or (
                    rec.is_reverse and rec.is_read2
                ):
                    currentReadDict[bc1 + "+" + bc2]["F1R2"] += 1
                else:
                    currentReadDict[bc1 + "+" + bc2]["F2R1"] += 1
            elif currentReadDict.get(bc2 + "+" + bc1) != None:
                currentReadDict[bc2 + "+" + bc1]["seqs"].append(rec)
                if (rec.is_forward and rec.is_read1) or (
                    rec.is_reverse and rec.is_read2
                ):
                    currentReadDict[bc2 + "+" + bc1]["F1R2"] += 1
                else:
                    currentReadDict[bc2 + "+" + bc1]["F2R1"] += 1
            else:
                currentReadDict.update({bc: {"seqs": [rec], "F1R2": 0, "F2R1": 0}})
                if (rec.is_forward and rec.is_read1) or (
                    rec.is_reverse and rec.is_read2
                ):
                    currentReadDict[bc1 + "+" + bc2]["F1R2"] += 1
                else:
                    currentReadDict[bc1 + "+" + bc2]["F2R1"] += 1

        else:
            """
            Calling block starts
            """
            for key in currentReadDict.keys():
                readSet = currentReadDict[key]["seqs"]
                setBc = key.split("+")
                setBc1 = setBc[0]
                setBc2 = setBc[1]
                F2R1 = currentReadDict[key]["F2R1"]
                F1R2 = currentReadDict[key]["F1R2"]
                duplex_no = f"{min([F1R2,F2R1])}+{max([F1R2,F2R1])}"
                if duplex_read_num_dict.get(duplex_no) is None:
                    duplex_read_num_dict[duplex_no] = [0,0]
                if setBc1 != setBc2 and F2R1 > 1 and F1R2 > 1:
                    indel_bool = [
                        ("I" in seq.cigarstring or "D" in seq.cigarstring)
                        for seq in readSet
                    ]
                    if any(indel_bool):
                        genotypeDSIndel(readSet,tumorBam,params)

                    else:
                        if isLearn:
                            if sum([int(seq.get_tag('NM')) for seq in readSet]) == 0:
                                continue
                        rs_reference_end = max([r.reference_end for r in readSet])
                        chromNow = readSet[0].reference_name
                        if (
                            chromNow != reference_mat_chrom
                            or rs_reference_end >= reference_mat_end
                        ):
                            ### Output coverage
                            if "coverage" in locals():
                                non_zero_positions = np.nonzero(coverage)
                                for pos in non_zero_positions[0].tolist():
                                    locus_bed.write(
                                        (
                                            "\t".join(
                                                [
                                                    reference_mat_chrom,
                                                    str(pos),
                                                    str(pos + 1),
                                                    str(coverage[pos]),
                                                ]
                                            )
                                            + "\n"
                                        ).encode("utf-8")
                                    )
                                total_coverage += np.sum(coverage)
                            reference_mat_chrom = chromNow
                            reference_mat_start = readSet[0].reference_start
                            try:
                                region_end = region[2]
                            except:
                                region_end = 10e10
                            contig_end = tumorBam.get_reference_length(chromNow)
                            reference_mat_end = min(
                                readSet[0].reference_start + 100000,
                                max(region_end, readSet[0].reference_end),
                                contig_end,
                            )
                            (
                                prior_mat,
                                snp_mask,
                                indel_mask,
                                noise_mask,
                                n_cov_mask,
                                ref_np,
                            ) = prepare_reference_mats(
                                reference_mat_chrom,
                                reference_mat_start,
                                reference_mat_end,
                                fasta[reference_mat_chrom][
                                    reference_mat_start:reference_mat_end
                                ],
                                germline,
                                noise,
                                ncoverage,
                                params,
                                )
                            coverage = np.zeros(100000, dtype=int)
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
                        reference_length = min(
                            [read.reference_length for read in readSet]
                        )
                        end_ind = (
                            readSet[0].reference_start
                            + reference_length
                            - reference_mat_start
                        )
                        masks = np.zeros([4, end_ind - start_ind], dtype=bool)
                        masks[0, :] = snp_mask[start_ind:end_ind]
                        masks[1, :] = noise_mask[start_ind:end_ind]
                        masks[2, :] = n_cov_mask[start_ind:end_ind]
                        left, right = determineTrimLength(
                            readSet[0], params=params, processed_flag=processed_flag
                        )
                        masks[3, :left] = True
                        masks[3, -right:] = True
                        antimask = np.all(~masks, axis=0)
                        ### If the whole reads are masked:
                        if not np.any(antimask):
                            continue
                        ### Calculate genotype probability
                        (
                            F1R2_ARLR,
                            F2R1_ARLR,
                            ref_int,
                            alt_int,
                            antimask,
                            F1R2_count,
                            F2R1_count,
                        ) = genotypeDSSnv(
                            readSet,
                            ref_np[start_ind:end_ind],
                            prior_mat[start_ind:end_ind, :],
                            antimask,
                            params,
                        )
                        if isLearn:
                            mismatch_now = profileTriNucMismatches(
                                readSet,
                                ref_np[start_ind:end_ind],
                                params,
                            )
                            for trinuc in mismatch_dict.keys():
                                for nnn in range(4):
                                    mismatch_dict[trinuc][nnn] += mismatch_now[trinuc][nnn]
                        ref_int = ref_np[start_ind:end_ind]
                        refs_ind = np.nonzero(
                            np.logical_and(
                                F1R2_ARLR <= -params["pcutoff"],
                                F2R1_ARLR <= -params["pcutoff"],
                            )
                        )[0].tolist()
                        muts_ind = np.nonzero(
                            np.logical_and(
                                F1R2_ARLR >= params["pcutoff"],
                                F2R1_ARLR >= params["pcutoff"],
                            )
                        )[0].tolist()
                        pass_bool = np.full(F1R2_ARLR.size, False, dtype=bool)
                        pass_bool[refs_ind] = True
                        pass_bool[muts_ind] = True
                        pos = [
                            mut_ind + start_ind + reference_mat_start
                            for mut_ind in muts_ind
                        ]
                        # muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                        mut_positions = [
                            mut_ind + start_ind + reference_mat_start + 1
                            for mut_ind in muts_ind
                        ]
                        if len(mut_positions) > 1:
                            break_flag = 0
                            for nn in range(1, len(mut_positions)):
                                if mut_positions[nn] - mut_positions[nn - 1] != 1:
                                    break_flag = 1
                            if break_flag:
                                continue
                        NMs = [seq.get_tag("NM") for seq in readSet]
                        averageNM = sum(NMs) / len(NMs)
                        if averageNM - len(mut_positions) > 1:
                            continue
                        for nn in range(len(mut_positions)):
                            mut_chrom = reference_mat_chrom
                            mut_pos = mut_positions[nn]
                            mut_ref = num2base[ref_int[muts_ind[nn]]]
                            mut_alt = num2base[alt_int[muts_ind[nn]]]
                            if muts_dict.get(
                                "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                            ):
                                continue
                            na, nr, ndp = extractDepthSnv(
                                normalBam, mut_chrom, mut_pos, mut_ref, mut_alt
                            )
                            if na > params["maxAltCount"]:
                                continue
                            ta, tr, tdp = extractDepthSnv(
                                tumorBam, mut_chrom, mut_pos, mut_ref, mut_alt
                            )
                            if readSet[0].is_forward:
                                readPos5p = min(
                                    muts_ind[nn] + 1,
                                    abs(readSet[0].template_length) - muts_ind[nn],
                                )
                                readPos3p = min(
                                    abs(readSet[0].reference_length) - muts_ind[nn],
                                    muts_ind[nn] + 1
                                )
                            else:
                                readPos5p = min(
                                    readSet[0].reference_length - muts_ind[nn],
                                    abs(readSet[0].template_length)
                                    - readSet[0].reference_length
                                    + muts_ind[nn],
                                )
                                readPos3p = min(
                                    abs(readSet[0].reference_length) - muts_ind[nn],
                                    muts_ind[nn] + 1
                                )
                            FPs.append(readPos5p)
                            RPs.append(readPos3p)

                            mut = {
                                "chrom": mut_chrom,
                                "pos": mut_pos,
                                "ref": mut_ref,
                                "alt": mut_alt,
                                "infos": {
                                    "RP5": readPos5p,
                                    "RP3": readPos3p,
                                    "F1R2": F1R2,
                                    "F2R1": F2R1,
                                    "TG": F1R2_ARLR[muts_ind[nn]],
                                    "BG": F2R1_ARLR[muts_ind[nn]],
                                    "TC": ",".join(
                                        [
                                            str(_)
                                            for _ in F1R2_count[
                                                :, muts_ind[nn]
                                            ].tolist()
                                        ]
                                    ),
                                    "BC": ",".join(
                                        [
                                            str(_)
                                            for _ in F2R1_count[
                                                :, muts_ind[nn]
                                            ].tolist()
                                        ]
                                    ),
                                },
                                "formats": ["AC", "RC", "DP"],
                                "samples": [[ta, tr, tdp], [na, nr, ndp]],
                            }
                            muts_dict[
                                "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                            ] = 0
                            muts.append(mut)
                        if isLearn:
                            continue
                        coverage[start_ind:end_ind][pass_bool] += 1
                        duplex_read_num_dict[duplex_no][1] += np.count_nonzero(pass_bool)
                        duplex_read_num_dict[duplex_no][0] += 1
                        duplex_count += 1
            """
            Calling block ends
            """
            currentReadDict = {bc: {"seqs": [rec], "F1R2": 0, "F2R1": 0}}
            if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
                currentReadDict[bc1 + "+" + bc2]["F1R2"] += 1
            else:
                currentReadDict[bc1 + "+" + bc2]["F2R1"] += 1
            currentStart = start

    """
    Calling block starts
    """
    for key in currentReadDict.keys():
        readSet = currentReadDict[key]["seqs"]
        setBc = key.split("+")
        setBc1 = setBc[0]
        setBc2 = setBc[1]
        F2R1 = currentReadDict[key]["F2R1"]
        F1R2 = currentReadDict[key]["F1R2"]
        duplex_no = f"{min([F1R2,F2R1])}+{max([F1R2,F2R1])}"
        if duplex_read_num_dict.get(duplex_no) is None:
            duplex_read_num_dict[duplex_no] = [0,0]
        if setBc1 != setBc2 and F2R1 > 1 and F1R2 > 1:
            indel_bool = [
                ("I" in seq.cigarstring or "D" in seq.cigarstring)
                for seq in readSet
            ]
            if any(indel_bool):
                mlmlml = 1
                # genotypeDSIndel(readSet,tumorBam,params)

            else:
                if isLearn:
                    if sum([int(seq.get_tag('NM')) for seq in readSet]) == 0:
                        continue
                rs_reference_end = max([r.reference_end for r in readSet])
                chromNow = readSet[0].reference_name
                if (
                    chromNow != reference_mat_chrom
                    or rs_reference_end >= reference_mat_end
                ):
                    ### Output coverage
                    if "coverage" in locals():
                        non_zero_positions = np.nonzero(coverage)
                        for pos in non_zero_positions[0].tolist():
                            locus_bed.write(
                                (
                                    "\t".join(
                                        [
                                            reference_mat_chrom,
                                            str(pos),
                                            str(pos + 1),
                                            str(coverage[pos]),
                                        ]
                                    )
                                    + "\n"
                                ).encode("utf-8")
                            )
                        total_coverage += np.sum(coverage)
                    reference_mat_chrom = chromNow
                    reference_mat_start = readSet[0].reference_start
                    try:
                        region_end = region[2]
                    except:
                        region_end = 10e10
                    contig_end = tumorBam.get_reference_length(chromNow)
                    reference_mat_end = min(
                        readSet[0].reference_start + 100000,
                        max(region_end, readSet[0].reference_end),
                        contig_end,
                    )
                    (
                        prior_mat,
                        snp_mask,
                        indel_mask,
                        noise_mask,
                        n_cov_mask,
                        ref_np,
                    ) = prepare_reference_mats(
                        reference_mat_chrom,
                        reference_mat_start,
                        reference_mat_end,
                        fasta[reference_mat_chrom][
                            reference_mat_start:reference_mat_end
                        ],
                        germline,
                        noise,
                        ncoverage,
                        params,
                        )
                    coverage = np.zeros(100000, dtype=int)
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
                reference_length = min(
                    [read.reference_length for read in readSet]
                )
                end_ind = (
                    readSet[0].reference_start
                    + reference_length
                    - reference_mat_start
                )
                masks = np.zeros([4, end_ind - start_ind], dtype=bool)
                masks[0, :] = snp_mask[start_ind:end_ind]
                masks[1, :] = noise_mask[start_ind:end_ind]
                masks[2, :] = n_cov_mask[start_ind:end_ind]
                left, right = determineTrimLength(
                    readSet[0], params=params, processed_flag=processed_flag
                )
                masks[3, :left] = True
                masks[3, -right:] = True
                antimask = np.all(~masks, axis=0)
                ### If the whole reads are masked:
                if not np.any(antimask):
                    continue
                ### Calculate genotype probability
                (
                    F1R2_ARLR,
                    F2R1_ARLR,
                    ref_int,
                    alt_int,
                    antimask,
                    F1R2_count,
                    F2R1_count,
                ) = genotypeDSSnv(
                    readSet,
                    ref_np[start_ind:end_ind],
                    prior_mat[start_ind:end_ind, :],
                    antimask,
                    params,
                )
                if isLearn:
                    mismatch_now = profileTriNucMismatches(
                        readSet,
                        ref_np[start_ind:end_ind],
                        params,
                    )
                    for trinuc in mismatch_dict.keys():
                        for nnn in range(4):
                            mismatch_dict[trinuc][nnn] += mismatch_now[trinuc][nnn]
                ref_int = ref_np[start_ind:end_ind]
                refs_ind = np.nonzero(
                    np.logical_and(
                        F1R2_ARLR <= -params["pcutoff"],
                        F2R1_ARLR <= -params["pcutoff"],
                    )
                )[0].tolist()
                muts_ind = np.nonzero(
                    np.logical_and(
                        F1R2_ARLR >= params["pcutoff"],
                        F2R1_ARLR >= params["pcutoff"],
                    )
                )[0].tolist()
                pass_bool = np.full(F1R2_ARLR.size, False, dtype=bool)
                pass_bool[refs_ind] = True
                pass_bool[muts_ind] = True
                pos = [
                    mut_ind + start_ind + reference_mat_start
                    for mut_ind in muts_ind
                ]
                # muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                mut_positions = [
                    mut_ind + start_ind + reference_mat_start + 1
                    for mut_ind in muts_ind
                ]
                if len(mut_positions) > 1:
                    break_flag = 0
                    for nn in range(1, len(mut_positions)):
                        if mut_positions[nn] - mut_positions[nn - 1] != 1:
                            break_flag = 1
                    if break_flag:
                        continue
                NMs = [seq.get_tag("NM") for seq in readSet]
                averageNM = sum(NMs) / len(NMs)
                if averageNM - len(mut_positions) > 1:
                    continue
                for nn in range(len(mut_positions)):
                    mut_chrom = reference_mat_chrom
                    mut_pos = mut_positions[nn]
                    mut_ref = num2base[ref_int[muts_ind[nn]]]
                    mut_alt = num2base[alt_int[muts_ind[nn]]]
                    if muts_dict.get(
                        "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                    ):
                        continue
                    na, nr, ndp = extractDepthSnv(
                        normalBam, mut_chrom, mut_pos, mut_ref, mut_alt
                    )
                    if na > params["maxAltCount"]:
                        continue
                    ta, tr, tdp = extractDepthSnv(
                        tumorBam, mut_chrom, mut_pos, mut_ref, mut_alt
                    )
                    if readSet[0].is_forward:
                        readPos5p = min(
                            muts_ind[nn] + 1,
                            abs(readSet[0].template_length) - muts_ind[nn],
                        )
                        readPos3p = min(
                            abs(readSet[0].reference_length) - muts_ind[nn],
                            muts_ind[nn] + 1
                        )
                    else:
                        readPos5p = min(
                            readSet[0].reference_length - muts_ind[nn],
                            abs(readSet[0].template_length)
                            - readSet[0].reference_length
                            + muts_ind[nn],
                        )
                        readPos3p = min(
                            abs(readSet[0].reference_length) - muts_ind[nn],
                            muts_ind[nn] + 1
                        )
                    FPs.append(readPos5p)
                    RPs.append(readPos3p)

                    mut = {
                        "chrom": mut_chrom,
                        "pos": mut_pos,
                        "ref": mut_ref,
                        "alt": mut_alt,
                        "infos": {
                            "RP5": readPos5p,
                            "RP3": readPos3p,
                            "F1R2": F1R2,
                            "F2R1": F2R1,
                            "TG": F1R2_ARLR[muts_ind[nn]],
                            "BG": F2R1_ARLR[muts_ind[nn]],
                            "TC": ",".join(
                                [
                                    str(_)
                                    for _ in F1R2_count[
                                        :, muts_ind[nn]
                                    ].tolist()
                                ]
                            ),
                            "BC": ",".join(
                                [
                                    str(_)
                                    for _ in F2R1_count[
                                        :, muts_ind[nn]
                                    ].tolist()
                                ]
                            ),
                        },
                        "formats": ["AC", "RC", "DP"],
                        "samples": [[ta, tr, tdp], [na, nr, ndp]],
                    }
                    muts_dict[
                        "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                    ] = 0
                    muts.append(mut)
                if isLearn:
                    continue
                coverage[start_ind:end_ind][pass_bool] += 1
                duplex_read_num_dict[duplex_no][1] += np.count_nonzero(pass_bool)
                duplex_read_num_dict[duplex_no][0] += 1
                duplex_count += 1
    """
    Calling block ends
    """
    if isLearn:
        print(FPs,RPs,muts)
        return mismatch_dict,FPs,RPs
    return muts, total_coverage, recCount, duplex_count, duplex_read_num_dict
