import argparse
from gzip import open as gzopen
import os


#from itertools import izip
def trim(readPair,pattern):
    adapterLen = len(pattern)
    read1 = readPair[0]
    read2 = readPair[1]
    name1 = read1[0]
    name2 = read2[0]
    seq1 = read1[1]
    seq2 = read2[1]
    qual1 = read1[2]
    qual2 = read2[2]
    #bc1 = ''.join([a for a,b in zip(seq1[0:adapterLen],pattern) if b == 'N'])
    #bc2 = ''.join([a for a,b in zip(seq2[0:adapterLen],pattern) if b == 'N'])
    bc1 = ''.join([base for nn,base in enumerate(seq1[0:adapterLen]) if pattern[nn] == 'N'])
    bc2 = ''.join([base for nn,base in enumerate(seq2[0:adapterLen]) if pattern[nn] == 'N'])
    if len(bc1) > 0 :
        namenew1 = name1.strip('\n').split(' ')[0]+'_'+bc1+'+'+bc2+' '+ ' '.join(name1.strip('\n').split(' ')[1:])+'\n'
        namenew2 = name2.strip('\n').split(' ')[0]+'_'+bc1+'+'+bc2+' '+ ' '.join(name2.strip('\n').split(' ')[1:])+'\n'
    else:
        namenew1 = name1
        namenew2 = name2
    seqnew1 = seq1[adapterLen:]
    seqnew2 = seq2[adapterLen:]
    qualnew1 = qual1[adapterLen:]
    qualnew2 = qual2[adapterLen:]
    return [namenew1,seqnew1,qualnew1],[namenew2,seqnew2,qualnew2]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim fastq sequences and move barcodes to name')
    parser.add_argument('-i','--fq',type=str,help='fastq file (read 1 if paired)')
    parser.add_argument('-i2','--fq2',type=str,help='read 2 fastq file')
    parser.add_argument('-p','--pattern',type=str,help='pattern of sequence barcode')
    parser.add_argument('-o','--output',type=str,help='prefix of the output fastq files')
    #parser.add_argument('-f','--force',type=bool,default=False,help='overwrite existing files')
    args = parser.parse_args()
    #if os.path.exists(args.output+'_1.fastq') or os.path.exists(args.output+'_2.fastq'):
        #raise ("Output files exist. Please use --f to overwrite")
    if args.fq[-3:] == '.gz':
        fq1 = gzopen(args.fq,'rt')
        fq2 = gzopen(args.fq2,'rt')
    else:
        fq1 = open(args.fq,'r')
        fq2 = open(args.fq2,'r')
    with open(args.output+'_1.fastq','w') as out1:
        out1.write('')
    with open(args.output+'_2.fastq','w') as out2:
        out2.write('')
    fq1Out = open(args.output+'_1.fastq','a')
    fq2Out = open(args.output+'_2.fastq','a')
    lineIndex = 0
    for line1,line2 in zip(fq1,fq2):
        #print(lineIndex)
        #print(line1,line2)
        if lineIndex == 0:
            name1 = line1
            name2 = line2
            lineIndex = 1
        elif lineIndex == 1:
            seq1 = line1
            seq2 = line2
            lineIndex = 2
        elif lineIndex == 2:
            lineIndex = 3
        else:
            qual1 = line1
            qual2 = line2
            readPair = [[name1,seq1,qual1],\
                        [name2,seq2,qual2]]
            #barcodeIndex = [nn for nn,a in enumerate(args.pattern) if a == 'N']
            #print(barcodeIndex)
            read1,read2 = trim(readPair, args.pattern)
            lineIndex = 0     
            #print(read1,read2)
            fq1Out.write(read1[0]+read1[1]+'+\n'+read1[2])
            fq2Out.write(read2[0]+read2[1]+'+\n'+read2[2])
    fq1Out.close()
    fq2Out.close()
    fq1.close()
    fq2.close()




    
