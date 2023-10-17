#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from DCutils.splitBamRegions import splitBamRegions
from DCutils.callBamRegionFilter import callBam
from DCutils.utils import createVcfStrings

import argparse
from collections import OrderedDict
from Bio import SeqIO
import os
from multiprocessing import Pool
from pysam import AlignmentFile as BAM
import pysam
import time


if __name__ == "__main__":

    """
    Parse Arguments
    """
    parser = argparse.ArgumentParser(description='Call any mutations from a bam file')
    parser.add_argument('-b','--bam',type=str,help='bam file')
    parser.add_argument('-g','--germline',type=str,help='indexed germline af vcf')
    parser.add_argument('-f','--reference',type=str,help='reference fasta file')
    parser.add_argument('-o','--output',type=str,help='prefix of the output files')
    parser.add_argument('-r','--regions',nargs='+',type=str,help='contigs to consider for variant calling',default = ['chr'+str(_) for _ in range(1,23,1)] + ['chrX','chrY'])
    parser.add_argument('-p','--threads',type=int,help='prefix of the output files',default = 1)
    parser.add_argument('-ae','--amperr',type=float,help='estimated polymerase error rate',default = 1E-5)
    parser.add_argument('-mr','--mutRate',type=float,help='estimated somatic mutation rate per base',default = 2.5E-7)
    parser.add_argument('-t','--threshold',type=float,help='log likelihood ratio threshold of making a mutation call',default = 2)
    parser.add_argument('-mq','--mapq',type=float,help='minumum mapq for an alignment to be considered',default = 30)
    #parser.add_argument('-da','--damage',type=float,default=5E-7)


    parser.add_argument('-n','--normalBam',type=str,help='deduplicated normal bam')
    parser.add_argument('-m','--noise',type=str,help='noise mask')
    parser.add_argument('-c','--ncoverage',type=str,help='coverage bed of ') 
    parser.add_argument('-tt','--trimF',type=int,help='ignore mutation if it is less than n bps from ends of template',default=30)
    parser.add_argument('-tr','--trimR',type=int,help='ignore mutation if it is less than n bps from ends of read',default=8)
    parser.add_argument('-d','--minNdepth',type=int,help='minumum coverage in normal for called variants',default=10)
    parser.add_argument('-nad','--maxAltCount',type=int,default=0)
    parser.add_argument('-maf','--maxAF',type=int,default=0.1)

    args = parser.parse_args()

    """
    Store Parameters
    """
    params = {"tumorBam": args.bam,\
            "normalBam": args.normalBam,\
            "germline": args.germline,\
            "reference": args.reference,\
            "output": args.output,\
            "regions": args.regions,\
            "threads": args.threads,\
            "amperr": args.amperr,\
            "mutRate": args.mutRate,\
            "pcutoff": args.threshold,\
            "mapq": args.mapq,\
            "noise": args.noise,\
            "trim5": args.trimF,\
            "trim3": args.trimR,\
            #"trim5DBS": args.trim5DBS,\
            #"trim3DBS": args.trim3DBS,\
            "minNdepth": args.minNdepth,\
            "maxAltCount": args.maxAltCount,\
            "ncoverage": args.ncoverage,\
            #"damage": args.damage
             }
    """
    Initialze run
    """
    print("..............Loading reference genome.....................")
    fasta = SeqIO.to_dict(SeqIO.parse(args.reference,'fasta'))
    startTime = time.time()
    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    bamObject = BAM(args.bam,'rb')
    
    """
    Execulte variant calling
    """
    if args.threads == 1:
        """
        Single-thread execution
        """
        print(".........Starting variant calling..............")
        #contigs = [(r.strip('\n'),) for r in open(args.regions,'r').readlines()] # Only process contigs in region file
        paramsNow = params
        paramsNow['reference'] = fasta
        regions = params['regions']
        paramsNow['regions'] = [(chrom,0,bamObject.get_reference_length(chrom)-1) for chrom in regions]
        mutsAll,coverage,rec_num,duplex_num,duplex_read_num_single = callBam(paramsNow,0,1000000)
        muts_num = len(mutsAll)
        duplex_combinations = list(duplex_read_num_single.keys())
        duplex_combinations.sort()
        duplex_read_num = OrderedDict({duplex_read_num_single[num] for num in duplex_combinations})
    else:
        """
        Multi-thread execution
        """
        #contigs = [r.strip('\n') for r in open(args.regions,'r').readlines()] # Only process contigs in region file
        contigs = args.regions
        contigLengths = [bamObject.get_reference_length(contig) for contig in contigs]
        print("...........Spliting genomic regions for parallel execution................")
        cutSites,chunkSize = splitBamRegions(args.bam,args.threads,contigs) #Split the whole genome for parallel execution
        regionSequence = []
        currentContigIndex = 0

        """
        Determine regions for each process
        """

        for nn,site in enumerate(cutSites[1:]):
            pSite = cutSites[nn]
            if site[0] == pSite[0]:
                regionSequence.append((contigs[site[0]],pSite[1],site[1]))
            else:
                if pSite[1] != 0:
                    regionSequence.append((contigs[pSite[0]],pSite[1]))
                else:
                    regionSequence.append((contigs[pSite[0]],))
                for ii in range(pSite[0]+1,site[0]):
                    regionSequence.append((contigs[ii],))
                regionSequence.append((contigs[site[0]],0,site[1]))
        regionSequence.append((contigs[site[0]],site[1]))
        for ii in range(site[0]+1,len(contigs)):
            regionSequence.append((contigs[ii],))
        print("............Completed region splitting in " + str((time.time()-startTime)/60) + " minutes,loading reference genome..................")

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
            #print(regions)
            chroms = [region[0] for region in regions]
            fastaNow = {chrom:fasta[chrom] for chrom in chroms} # Takes partial fasta to reduce memory usage
            paramsNow = params
            paramsNow['reference'] = fastaNow
            paramsNow['regions'] = regions
            callArgument = (paramsNow.copy(),nn,chunkSize)
            callArguments.append(callArgument)
            regions = []
        results = pool.starmap(callBam,callArguments)# each result return three list: number of duplex reads, effective lengths, list of mutations
        muts = [_[0] for _ in results]
        coverages = [_[1] for _ in results]
        rec_nums = [_[2] for _ in results]
        duplex_nums = [_[3] for _ in results]
        duplex_read_nums = [_[4] for _ in results]
        print("..............Completed bam calling in " + str((time.time()-startTime2)/60) + " minutes,merging results................." )
        pool.close()
        pool.terminate()
        pool.join()
        mutsAll = sum(muts,[])
        muts_num = len(mutsAll)
        coverage = sum(coverages)
        rec_num = sum(rec_nums)
        duplex_num = sum(duplex_nums)

        duplex_combinations = list(set.union(*[set(d.keys()) for d in duplex_read_nums]))
        duplex_combinations.sort()
        duplex_read_num = OrderedDict({num:sum([d.get(num,0) for d in duplex_read_nums]) for num in duplex_combinations})

    tBam = BAM(args.bam,'rb')
    contigs = tBam.references
    #print(contigs)
    chromDict = {contig:tBam.get_reference_length(contig) for contig in contigs}

    """
    infoDict = {"SS":[1,"Float","Somatic Score of variant"],"NM":[1,"Float","mean NM of alt read group"],"AS":[1,"Float","mean AS of the read group"],"XS":[1,"Float","mean XS of the read group"],\
    "RP5":[1,"Integer","read position"],"RP3":[1,"Integer","distance from 3p"]}
    formatDict = {"AC":[1,"Integer","Count of alt allele"],"RC":[1,"Integer","Count of ref allele"],"DP":[1,"Integer","Depth at the location"],\
    "DC":[1,"Integer","Number of reads in the duplex group"]}
    filterDict = {"PASS":"All filter Passed"}
    """
    infoDict = {"F1R2":[1,"Integer","Number of F1R2 read(s) in the read bundle"], "F2R1":[1,"Integer","Number of F2R1 read(s) in the read bundle"],\
    "RP5":[1,"Integer","read position"],"RP3":[1,"Integer","distance from 3p"], \
    "TG":[1,"Float","Alt/Ref log likelihood ratio of top strand"],"BG":[1,"Float","Alt/Ref log likelihood ratio of bottom strand"],\
    "TC":[4,"Integer","Top strand base count"],"BC":[4,"Float","Bottom strand base count"],\
    "BC":[4,"Integer","Top strand base count"],"BC":[4,"Float","Bottom strand base count"]}
    formatDict = {"AC":[1,"Integer","Count of alt allele"],"RC":[1,"Integer","Count of ref allele"],"DP":[1,"Integer","Depth at the location"]}
    filterDict = {"PASS":"All filter Passed"} 
    vcfLines = createVcfStrings(chromDict,infoDict,formatDict,filterDict,mutsAll)
    with open(args.output + '.vcf','w') as vcf:
        vcf.write(vcfLines)
    
    burden = muts_num/coverage
    efficiency = duplex_num/rec_num

    with open(params["output"]+"_summary.txt",'w') as f:
        f.write(f"Number of Mutations\t{muts_num}\n")
        f.write(f"Effective Coverage\t{coverage}\n")
        f.write(f"Estimated Burden\t{burden}\n")
        f.write(f"Duplex read number\t{duplex_num}\n")
        f.write(f"Efficiency\t{efficiency}\n")
    
    with open(params["output"]+"_duplex_group_stats.txt",'w') as f:
        f.write(f"Duplex Group Strand Composition\tDuplex Group Number\n")
        for read_num in duplex_read_num.keys():
            f.write(f"{read_num}\t{duplex_read_num[read_num]}\n")

    print("..............Completed variant calling " + str((time.time()-startTime)/60) + " minutes..............." )