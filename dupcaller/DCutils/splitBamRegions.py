import argparse
from pysam import AlignmentFile as BAM
from multiprocessing import Queue,Process,Lock
import math

def countBamWorker(jobQueue,jobLock,resultQueue):
    while True:
        with jobLock:
            job = jobQueue.get()
            if job == 1: break
            bam = job[0]
            chrom = job[1]
            bamObject = BAM(bam,'rb')
            count = bamObject.count(chrom)
            print("Finished read counting for "+chrom)
            resultQueue.put({chrom:count})

def splitBamRegions(bam,num,contigs):
    #jobQueue = Queue()
    #jobLock = Lock()
    #resultQueue = Queue()
    bamObject = BAM(bam,'rb')
    bamStats = bamObject.get_index_statistics()
    readNumber = 0
    contigCountsDict = {}
    for stat in bamStats:
        if stat.contig in contigs:
            readNumber += stat.total
            contigCountsDict.update({stat.contig:stat.total})
    #readNumber = bamObject.count()
    chunkSize = math.ceil(readNumber / num)
    tidList = [bamObject.get_tid(c) for c in contigs]


    contigCounts = [contigCountsDict[contig] for contig in contigs]
    contigCountsAll = [0 for c in range(bamObject.nreferences)]
    for nn,contig in enumerate(contigs):
        currentTid = bamObject.get_tid(contig)
        contigCountsAll[currentTid] = contigCounts[nn]
    contigCountsTuple = tuple(contigCountsAll) 
    #print(contigCounts)
    contigLengths = [bamObject.get_reference_length(contig) for contig in bamObject.references]
    chunkSize = int(sum(contigCounts)/num)
    #currentCount = 0
    #outBams = [BAM(outPrefix+'_'+str(nn)+".bam",'wb',header=bamObject.header) for nn in range(num)]
    #currentRegionIndex = 0
    cutSite = [(0,0)]
    leftChunkSize = chunkSize
    currentContigCount = 0
    tidPointer = 0
    while len(cutSite) < num:
        if contigCounts[0] == 0:
            contigCounts.pop(0)
            #tidList.pop(0)
            tidPointer += 1
            currentContigCount = 0
        elif leftChunkSize > contigCounts[0]:
            leftChunkSize -= contigCounts.pop(0)
            #tidList.pop(0)
            tidPointer += 1
            currentContigCount = 0
        else:
            currentContigCount += leftChunkSize
            contigCounts[0] -= leftChunkSize
            cutSite.append((tidPointer,int(float(currentContigCount)/float(contigCountsTuple[tidList[0]]) * contigLengths[tidList[0]])))
            leftChunkSize = chunkSize
        #print(cutSite)
    return cutSite,chunkSize

        