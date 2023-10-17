from subprocess import check_output
import subprocess
from scipy.stats import binom


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
        if cigar[0] in [3,4]:
            readPos += cigar[1]
        if cigar[0] == 1:
            pos = refPos
            sequence = seq.query_sequence[readPos:readPos + cigar[1]]
            quals = seq.query_qualities[readPos:readPos + cigar[1]]
            indels.update({str(pos)+'I':[sequence,quals,readPos]})
            readPos += cigar[1]
        if cigar[0] == 2:
            pos = refPos
            indels.update({str(pos)+'D':[cigar[1],readPos]})
            refPos += cigar[1]  
    return indels  
    
def extractDepthSnv(bam,chrom,pos,ref,alt):
    has_barcode = True
    for seq in bam.fetch(chrom,pos-1,pos):
        if len(seq.query_name.split('_')[-1].split('+')) != 2:
            has_barcode = False
        break
    if has_barcode:
        refReadBarcodes = {}
        altReadBarcodes = {}
        for pileupcolumn in bam.pileup(chrom,pos-1,pos,min_base_quality=0):
            if pileupcolumn.pos == pos-1:
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del or pileupread.is_refskip or pileupread.alignment.is_secondary or pileupread.alignment.is_supplementary:
                        continue
                    if pileupread.alignment.query_sequence[pileupread.query_position] == alt:
                        barcodes = pileupread.alignment.query_name.split('_')[-1].split('+')
                        if not altReadBarcodes.get(barcodes[0]+'+'+barcodes[1]) or \
                        not altReadBarcodes.get(barcodes[1]+'+'+barcodes[0]):
                            altReadBarcodes.update({barcodes[0]+'+'+barcodes[1]:1})
                    else:
                        barcodes = pileupread.alignment.query_name.split('_')[-1].split('+')
                        if not refReadBarcodes.get(barcodes[0]+'+'+barcodes[1]) or \
                        not refReadBarcodes.get(barcodes[1]+'+'+barcodes[0]):
                            refReadBarcodes.update({barcodes[0]+'+'+barcodes[1]:1})
        altAlleleCount = len(altReadBarcodes.keys())
        refAlleleCount = len(refReadBarcodes.keys())
        depth = altAlleleCount + refAlleleCount
    else:
        altAlleleCount = 0
        refAlleleCount = 0
        for pileupcolumn in bam.pileup(chrom,pos-1,pos,min_base_quality=0):
            if pileupcolumn.pos == pos-1:
                for pileupread in  pileupcolumn.pileups:
                    if pileupread.is_del or pileupread.is_refskip or pileupread.alignment.is_secondary or pileupread.alignment.is_supplementary or pileupread.alignment.is_duplicate:
                        continue
                    if pileupread.alignment.query_sequence[pileupread.query_position] == alt:
                        altAlleleCount += 1
                    else:
                        refAlleleCount += 1
        depth = altAlleleCount + refAlleleCount
    return altAlleleCount,refAlleleCount,depth

def extractDepthIndel(bam,chrom,pos,ref,alt):
    if len(ref) == 1:
        indelPos = str(pos)+'I'
    else:
        indelPos = str(pos)+'D'
    refReadBarcodes = {}
    altReadBarcodes = {}
    for read in bam.fetch(chrom,pos-1,pos):
        matchFlag = 0
        if read.mapping_quality <= 30: continue
        if 'I' in read.cigarstring or 'D' in read.cigarstring:
            indels = findIndels(read)
            if indels.get(indelPos):
                matchFlag = 1
        if matchFlag == 1:
            barcodes = read.query_name.split('_')[-1].split('+')
            if not altReadBarcodes.get(barcodes[0]+'+'+barcodes[1]) or \
            not altReadBarcodes.get(barcodes[1]+'+'+barcodes[0]):
                altReadBarcodes.update({barcodes[0]+'+'+barcodes[1]:1})
        else:
            barcodes = read.query_name.split('_')[-1].split('+')
            if not refReadBarcodes.get(barcodes[0]+'+'+barcodes[1]) or \
            not refReadBarcodes.get(barcodes[1]+'+'+barcodes[0]):
                refReadBarcodes.update({barcodes[0]+'+'+barcodes[1]:1})
    altAlleleCount = len(altReadBarcodes.keys())
    refAlleleCount = len(refReadBarcodes.keys())
    depth = altAlleleCount + refAlleleCount
    return altAlleleCount,refAlleleCount,depth

def createVcfStrings(chromDict,infoDict,formatDict,filterDict,recs):
    lines = ["##fileformat=VCFv4.2"]
    for filterr in filterDict.keys():
        lines.append("##FILTER=<ID={},Description=\"{}\">".format(filterr,filterDict[filterr]))
    for info in infoDict.keys():
        lines.append("##INFO=<ID={},Number={},Type={},Description=\"{}\">".format(info,infoDict[info][0],infoDict[info][1],infoDict[info][2]))
    for form in formatDict.keys():
        lines.append("##FORMAT=<ID={},Number={},Type={},Description=\"{}\">".format(form,formatDict[form][0],formatDict[form][1],formatDict[form][2]))
    for chrom in chromDict.keys():
        lines.append("##contig=<ID={},length={}>".format(chrom,chromDict[chrom]))
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL")
    for rec in recs:
        chrom = rec['chrom']
        pos = rec['pos']
        alt = rec['alt']
        ref = rec['ref']
        infos = rec['infos']
        formats = rec['formats']
        samples = rec["samples"]
        lineEntries = [chrom,str(pos),'.',ref,alt,".","PASS",";".join(["{info}={val}".format(info=info,val=infos[info]) for info in infoDict.keys()]),':'.join(formats),":".join([str(s) for s in samples[0]]),":".join([str(s) for s in samples[1]])]
        lines.append('\t'.join(lineEntries))
    return '\n'.join(lines)+"\n"

def checkCallablity(F1R2,F2R1,F1R2BQ,F2R1BQ,mutRate,perror,pcut):
    mockFBaseProb = genotypeProbSnv(['A']*F1R2,F1R2BQ,[mutRate,mutRate,mutRate,1-3*mutRate],perror)
    mockRBaseProb = genotypeProbSnv(['A']*F2R1,F2R1BQ,[mutRate,mutRate,mutRate,1-3*mutRate],perror)
    if mockFBaseProb[3] <= pcut and mockRBaseProb[3] <= pcut:
        return True
    return False


#def germlineProbilityTest()
#def extractContext(readSet,trim3,trim5,fasta):
    #contextList = []
    #for pos in range(readSet.pivotSeq.reference_start + trim3,readSet.pivotSeq.reference_end - trim5):
        #contextList += [fasta[pos-1:pos+2]]
    #return contextList
