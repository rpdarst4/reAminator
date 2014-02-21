#!/usr/bin/python2.7
# (c) 2013, Russell Darst, University of Florida

"""
* identifies barcodes and records them in sequence IDs
*   as >seqID|barcode
* can name barcodes by sequence (default) or ID
* RPD4 2013
"""

from Bio import SeqIO
from collections import defaultdict
from os import environ,mkdir,path,remove,tempnam
from subprocess import call
import argparse,ConfigParser,warnings

##############################################################################

# dependencies

config=ConfigParser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)),'reAminator.cfg'))

# to supress warnings from use of tempnam
warnings.filterwarnings('ignore')

cmd='{0} -grep -noindex -in {1} -out {2} -max-mismatch {3} -left {4} -seq {5}'

##############################################################################

def make(reads,barcodes,saveAs,name=False,mismatch=0,report=False):

    # fastools does not support pipe
    temp=tempnam(environ['HOME'])
    args=[config.get('paths','fastools'),reads,temp,mismatch]
    calls={seq.id: [] for seq in list(SeqIO.parse(reads,'fasta'))}
    counts=defaultdict(lambda: defaultdict(int))
    for barcode in SeqIO.parse(barcodes,'fasta'):

        # Fastools to check left of each sequence for barcode match
        call(cmd.format(*args+[-len(barcode),barcode.seq]).split())
        for seq in SeqIO.parse(temp,'fasta'): calls[seq.id]+=[barcode]
        remove(temp)
    with open(saveAs,'w') as handle:
        for seq in SeqIO.parse(reads,'fasta'):
            if calls[seq.id]:

                # only use best barcode
                barcode=sorted(calls[seq.id],key=lambda x: sum([
                    seq.seq.find(j,i)-i for i,j in enumerate(x)])-len(x))[0]
                seq.id+='|{0}'.format((barcode.seq,barcode.id)[name])
                counts[barcode][str(seq.seq[:len(barcode)+1])]+=1
            SeqIO.write(seq,handle,'fasta')

    if report:
        with open(path.splitext(saveAs)[0]+'.report.txt','w') as handle:
            handle.write('\n\n'.join(['\n'.join(
                ['{0}: {1}'.format(k.id,k.seq)]+
                ['{0} x {1}'.format(i,j) for i,j in sorted(v.items())])
                                      for k,v in counts.items()]))
                                
##############################################################################

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='')
    parser.add_argument('reads',help='FASTA-format file of 454 reads')
    parser.add_argument('codes',help='FASTA-format file of barcodes')
    parser.add_argument('-name',action='store_true',
                        help='track barcodes by ID (default by sequence)')
    parser.add_argument('-miss',default=0,help='number mismatches allowed')
    parser.add_argument('-report',action='store_true',help='record # aligned')
    parser.add_argument('-save',help='output file')
    args=parser.parse_args()
    if not args.save:
        args.save=path.splitext(args.reads)[0]+'.'+path.basename(args.codes)
    make(args.reads,args.codes,args.save,args.name,args.miss,args.report)
