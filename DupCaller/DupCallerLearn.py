#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from DCutils.splitBamRegions import splitBamRegions
from DCutils.call import callBam
from DCutils.funcs import createVcfStrings

import argparse
from collections import OrderedDict
from Bio import SeqIO
import os
from multiprocessing import Pool
from pysam import AlignmentFile as BAM
import pysam
import time
from matplotlib import pyplot as plt


if __name__ == "__main__":
    """
    Parse Arguments
    """
    parser = argparse.ArgumentParser(description="Call any mutations from a bam file")
    parser.add_argument("-b", "--bam", type=str, help="bam file")
    parser.add_argument("-g", "--germline", type=str, help="indexed germline af vcf")
    parser.add_argument("-f", "--reference", type=str, help="reference fasta file")
    parser.add_argument("-o", "--output", type=str, help="prefix of the output files")
    parser.add_argument(
        "-r",
        "--regions",
        nargs="+",
        type=str,
        help="contigs to consider for variant calling",
        default=["chr" + str(_) for _ in range(1, 23, 1)] + ["chrX", "chrY"],
    )
    parser.add_argument(
        "-p", "--threads", type=int, help="prefix of the output files", default=1
    )
    parser.add_argument(
        "-ae",
        "--amperr",
        type=float,
        help="estimated polymerase error rate",
        default=2e-4,
    )
    parser.add_argument(
        "-mr",
        "--mutRate",
        type=float,
        help="estimated somatic mutation rate per base",
        default=2.5e-7,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        help="log likelihood ratio threshold of making a mutation call",
        default=2,
    )
    parser.add_argument(
        "-mq",
        "--mapq",
        type=float,
        help="minumum mapq for an alignment to be considered",
        default=30,
    )
    # parser.add_argument('-da','--damage',type=float,default=5E-7)

    parser.add_argument("-n", "--normalBam", type=str, help="deduplicated normal bam")
    parser.add_argument("-m", "--noise", type=str, help="noise mask")
    parser.add_argument("-c", "--ncoverage", type=str, help="coverage bed of ")
    parser.add_argument(
        "-tt",
        "--trimF",
        type=int,
        help="ignore mutation if it is less than n bps from ends of template",
        default=0,
    )
    parser.add_argument(
        "-tr",
        "--trimR",
        type=int,
        help="ignore mutation if it is less than n bps from ends of read",
        default=15,
    )
    parser.add_argument(
        "-d",
        "--minNdepth",
        type=int,
        help="minumum coverage in normal for called variants",
        default=10,
    )
    parser.add_argument("-nad", "--maxAltCount", type=int, default=0)
    parser.add_argument("-maf", "--maxAF", type=int, default=0.1)
    parser.add_argument("-rl", "--readLen", type=int, default=134)
    args = parser.parse_args()

    """
    Store Parameters
    """
    params = {
        "tumorBam": args.bam,
        "normalBam": args.normalBam,
        "germline": args.germline,
        "reference": args.reference,
        "output": args.output,
        "regions": args.regions,
        "threads": args.threads,
        "amperr": args.amperr,
        "mutRate": args.mutRate,
        "pcutoff": args.threshold,
        "mapq": args.mapq,
        "noise": args.noise,
        "trim5": args.trimF,
        "trim3": args.trimR,  # "trim5DBS": args.trim5DBS,\
        # "trim3DBS": args.trim3DBS,\
        "minNdepth": args.minNdepth,
        "maxAltCount": args.maxAltCount,
        "ncoverage": args.ncoverage,  # "damage": args.damage
    }
    """
    Initialze run
    """
    print("..............Loading reference genome.....................")
    fasta = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    startTime = time.time()
    if not os.path.exists("tmp"):
        os.mkdir("tmp")
    bamObject = BAM(args.bam, "rb")

    """
    Execulte variant calling
    """
    if args.threads == 1:
        """
        Single-thread execution
        """
        print(".........Starting variant calling..............")
        # contigs = [(r.strip('\n'),) for r in open(args.regions,'r').readlines()] # Only process contigs in region file
        paramsNow = params
        paramsNow["reference"] = fasta
        paramsNow["isLearn"] = True
        regions = params["regions"]
        paramsNow["regions"] = [
            (chrom, 0, bamObject.get_reference_length(chrom) - 1) for chrom in regions
        ]
        mismatch,FPs,RPs = callBam(
            paramsNow, 0, 1000000
        )
    else:
        """
        Multi-thread execution
        """
        # contigs = [r.strip('\n') for r in open(args.regions,'r').readlines()] # Only process contigs in region file
        contigs = args.regions
        contigLengths = [bamObject.get_reference_length(contig) for contig in contigs]
        print(
            "...........Spliting genomic regions for parallel execution................"
        )
        cutSites, chunkSize = splitBamRegions(
            args.bam, args.threads, contigs
        )  # Split the whole genome for parallel execution
        regionSequence = []
        currentContigIndex = 0

        """
        Determine regions for each process
        """

        for nn, site in enumerate(cutSites[1:]):
            pSite = cutSites[nn]
            if site[0] == pSite[0]:
                regionSequence.append((contigs[site[0]], pSite[1], site[1]))
            else:
                if pSite[1] != 0:
                    regionSequence.append((contigs[pSite[0]], pSite[1]))
                else:
                    regionSequence.append((contigs[pSite[0]],))
                for ii in range(pSite[0] + 1, site[0]):
                    regionSequence.append((contigs[ii],))
                regionSequence.append((contigs[site[0]], 0, site[1]))
        regionSequence.append((contigs[site[0]], site[1]))
        for ii in range(site[0] + 1, len(contigs)):
            regionSequence.append((contigs[ii],))
        print(
            "............Completed region splitting in "
            + str((time.time() - startTime) / 60)
            + " minutes,loading reference genome.................."
        )

        """
        Start variant calling 
        """

        callArguments = []
        startTime2 = time.time()
        print(".........Starting variant calling..............")
        pool = Pool()
        for nn in range(args.threads):
            regions = []
            while len(regionSequence) != 0:
                if len(regionSequence[0]) != 3:
                    regions.append(regionSequence.pop(0))
                else:
                    regions.append(regionSequence.pop(0))
                    break
            # print(regions)
            chroms = [region[0] for region in regions]
            fastaNow = {
                chrom: fasta[chrom] for chrom in chroms
            }  # Takes partial fasta to reduce memory usage
            paramsNow = params
            paramsNow["reference"] = fastaNow
            paramsNow["regions"] = regions
            paramsNow["isLearn"] = True
            callArgument = (paramsNow.copy(), nn, chunkSize)
            callArguments.append(callArgument)
            regions = []
        results = pool.starmap(
            callBam, callArguments
        )  # each result return three list: number of duplex reads, effective lengths, list of mutations
        print(
            "..............Completed bam calling in "
            + str((time.time() - startTime2) / 60)
            + " minutes,merging results................."
        )
        print(results)
        pool.close()
        pool.terminate()
        pool.join()
        mismatch_dicts = [_[0] for _ in results]
        FPs = sum([_[1] for _ in results],[])
        RPs = sum([_[2] for _ in results],[])
        mismatch_dict = dict()
        #print(mismatch_dicts)
        for minus_base in ['A','T','C','G']:
            for ref_base in ['C','T']:
                    for plus_base in ['A','T','C','G']:
                        mismatch_dict[minus_base+ref_base + plus_base] = [0,0,0,0]
        for trinuc in mismatch_dicts[0].keys():
            for nnn in range(4):
                mismatch_dict[trinuc][nnn] = sum(_[trinuc][nnn] for _ in mismatch_dicts)
    with open(params["output"] + "/"+args.output+"_mismatch_trinuc_profile.txt", "w") as f:
        f.write("\tA\tT\tC\tG\n")
        for trinuc in mismatch_dict.keys():
            f.write(f"{trinuc}\t{mismatch_dict[trinuc][0]}\t{mismatch_dict[trinuc][1]}\t{mismatch_dict[trinuc][2]}\t{mismatch_dict[trinuc][3]}\n")
    FPs_count = [0 for _ in range(args.readLen)]
    RPs_count = [0 for _ in range(args.readLen)] 
    for nn in range(args.readLen):
        FPs_count[nn] = FPs.count(nn+1)
        RPs_count[nn] = RPs.count(nn+1)
    with open(params["output"] + "/"+args.output+"_DBS_end_profile.txt", "w") as f:
        f.write("Distance\tMutations_fragment_end\tMutations_read_end\n")
        for nn in range(args.readLen):
            f.write(f"{nn+1}\t{FPs_count[nn]}\t{RPs_count[nn]}\n")   
    print(
        "..............Completed variant calling "
        + str((time.time() - startTime) / 60)
        + " minutes..............."
    )
