#!/usr/bin/python2.7
# (c) 2013, Russell Darst, University of Florida

"""
* aligns sequences to references using BLAST
* masks cytosines for unbiased alignment
* returns SQL table
* user can choose to align on either or both strands
* works for cytosine methylation in any context (CG, GC, etc.)
* uses fastools, makeblastdb, blastn
* RPD4 2013
"""

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from os import environ,path,remove,tempnam
from subprocess import call
from time import time
import argparse,ConfigParser
import sqlite3,warnings

##############################################################################

# dependencies

config=ConfigParser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)),'reAminator.cfg'))
bn  =config.get('paths','blastn')
c2bm=config.get('paths','convert2blastmask')
ft  =config.get('paths','fastools')
mbdb=config.get('paths','makeblastdb')

##############################################################################

# to swap extensions
suffix=lambda x,y: path.splitext(x)[0]+'.'+y

# to supress warnings from use of tempnam
warnings.filterwarnings('ignore')

# for ambiguous base calls
ambigIUPAC={'A': 'MRW', 'C': 'MSY', 'G': 'KRS', 'T': 'KWY',
            'K': 'GT','M': 'AC','R': 'AG','S': 'CG','W': 'AT','Y': 'CT'}
ambigIUPAC={i+j: (list(set(list(ambigIUPAC[i]))&set(list(ambigIUPAC[j])))+[
    'N'])[0] for i in 'ACGKMRSTWY' for j in 'ACGKMRSTWY' if i!=j}

##############################################################################

def check_IDs(records):
    names=[record.id for record in SeqIO.parse(records,'fasta')]
    indexed=path.splitext(records)[0]+'_indexed.fa'
    if len(set(names))<len(names):
        N=len(str(len(names)))
        with open(records,'r') as source:
            with open(indexed,'w') as handle:
                for n,rec in enumerate(SeqIO.parse(source,'fasta')):
                    rec.id=str(n).zfill(N)+'|'+rec.id
        return indexed
    else: return records
    
##############################################################################
    
def deaminate(source,strand='ab'):
    seqs,temp=[],tempnam(environ['HOME'])
    for i in strand:

        # fastools to convert strands
        call ([ft,'-in',source,'-out',temp,('-c-to-t','-g-to-a')[i=='b']])

        # Bio.SeqIO to label & orient
        for seq in SeqIO.parse(temp,'fasta'):

            # use A+ or B- strands
            if i=='b': seq.seq=seq.seq.reverse_complement()
            seq.id+='.'+i
            seq.description=''
            seqs.append(seq)
            
        remove(temp)
        remove(temp+'.idx')
        
    with open(suffix(source,strand),'w') as handle:
        SeqIO.write(seqs,handle,'fasta')

############################################################################## 

def parse(whichFile):
    for record in NCBIXML.parse(open(whichFile,'r')):
        EI=[(alignment.hsps[0].expect,i)
            for i,alignment in enumerate(record.alignments)]
        EI.sort()
        for e,i in EI:
            if e/10<=EI[0][0]:
                alignment=record.alignments[i]

                # only take best hsp: no overlapping alignments
                hsp=alignment.hsps[0]
                yield alignment.accession,record.query.split()[0],\
                      hsp.expect,hsp.sbjct_start,hsp.sbjct_end,\
                      hsp.query_start,hsp.query_end,hsp.sbjct,hsp.query

##############################################################################

def put_data(seqs,refs,dest,mask=False,strand='ab',update=False):
    I=tempnam(environ['HOME'])
    
    # ensure all identifiers are unique & count input sequences        
    seqs,refs=check_IDs(seqs),check_IDs(refs)

    # files: (D)atabase, (X)ML blastn output, (Q)uery, (S)ubject
    D,X=suffix(dest,'db'),suffix(seqs,'xml')
    Q,S=suffix(seqs,strand),suffix(refs,'ab')

    # initial setup (skip for second phase of paired-end)
    if not update:
        if path.exists(D): remove(D)

    if not path.exists(X):
        start=time()

        # deaminate query and reference sequences in silico prior to blastn
        print 'Converting bases','...',
        deaminate(seqs,strand)
        deaminate(refs)
        print 'in',time()-start,'\n','Aligning','...',
        
        start=time()
        run_BLAST(Q,S,X,mask)
        print 'in',time()-start
        
    else:
        print 'Warning: XML file already exists'
        print '(delete & run again for fresh alignment).'
        
    print 'Filing data','...',
    start=time()

    # database of initial sequences
    seqs=SeqIO.index_db(I,[seqs,refs],'fasta')

    # set up sql database
    conn=sqlite3.connect(D)
    curs=conn.cursor()
    curs.execute('CREATE TABLE IF NOT EXISTS '
                 'records (id INTEGER PRIMARY KEY NOT NULL, '
                 'read, expt, locus, expect REAL, strand, sequence)')
    curs.execute('CREATE TABLE IF NOT EXISTS '
                 'loci (id INTEGER PRIMARY KEY NOT NULL, locus, sequence)')
    curs.execute('CREATE INDEX IF NOT EXISTS loc_idx ON loci (locus)')
    curs.execute('CREATE INDEX IF NOT EXISTS els_idx ON records '
                 '(expt,locus,strand)')
    curs.execute('CREATE INDEX IF NOT EXISTS rls_idx ON records '
                 '(read,locus,strand)')
    
    for S,Q,e,s5,s3,q5,q3,sBS,qBS in parse(X):
        Q,S=Q.split('.'),S.split('.')
        qStr,sStr=Q[-1],S[-1]
        qID,sID='.'.join(Q[:-1]),'.'.join(S[:-1])

        # discard bogus alignments (i.e. wrong strand)
        if s5>s3: continue
        
        # look up sequences using SeqIO.index_db object
        try: qSeq,sSeq=seqs[qID],seqs[sID]
        except KeyError:
            print qID,sID
            continue

        # had reverse-complemented b strand data
        if qStr=='b': qSeq=qSeq.reverse_complement()
        if sStr=='b': sSeq=sSeq.reverse_complement()
        
        # An A strand query should align to either
        #   an A strand subject (A+), or
        #   the reverse-complement of a B strand subject (B-),
        #   and vice-versa.
        # If both are A, neither is flipped (A+, A+);
        #   if both are B, both are flipped (B-, B-);
        #   if only one is B, only that one is flipped (A+, B-).
        # Thus, we align A+ or B- to A+ or B-.
        # As the strand of the query may not be known, we can try both
        #   A+ (C to T) and B- (G to A, rev. comp.). An A strand query
        #   treated as a B strand, or vice versa, should yield 2-base
        #   nonsense (e.g. AGTA --> AATA); which should not align.
        
        # BLAST coordinates are inclusive
        seq=reaminate(*[list(x.upper()) for x in [
            sSeq[s5-1:s3],qSeq[q5-1:q3],sBS,qBS]])

        # fill out with dashes
        seq='-'*(s5-1)+seq+'-'*len(sSeq[s3:])

        # return to original orientation of reference
        if sStr=='b': seq=seq.reverse_complement()
        
        # store alignments
        qID,xpt=(qID.split('|')+[''])[:2]
        if update:
            done=False

            # paired end reads should have same sequence ID
            for n,old in curs.execute(
                'SELECT id, sequence FROM records WHERE read=? '
                'AND locus=? AND strand=?',[qID,sID,sStr]):
                new=list(old)
                for i,j in enumerate(str(seq).upper()):
                    if j in '-N' or j==new[i]: continue
                    elif new[i] in '-N' : new[i]=j
                    else:

                        # ambiguous base calls
                        new[i]=ambigIUPAC.get(''.join(sorted(new[i]+j)),'N')
                        
                curs.execute('UPDATE records SET sequence = ? '
                             'WHERE id = ?',[''.join(new),n])

                # go to next entry in parse(XML)
                done=True
            if done: continue

        # otherwise, new entry in database
        curs.execute('INSERT INTO records VALUES (NULL,?,?,?,?,?,?)',
                     [qID,xpt,sID,e,sStr,str(seq).upper()])

    print 'in',time()-start,'\n','Adding refs','...',
    start=time()
    for ref, in set(curs.execute(
        'SELECT locus FROM records GROUP BY locus'))\
        -set(curs.execute('SELECT locus FROM loci')):
        curs.execute('INSERT INTO loci VALUES (NULL, ?, ?)',
                     (seqs[ref].id,str(seqs[ref].seq)))
    print 'in',time()-start,'\n','Aligned {0} sequences to {1} loci'.format(
        len(list(curs.execute('SELECT read FROM records GROUP BY read'))),
        len(list(curs.execute('SELECT locus FROM records GROUP BY locus'))))
    print
    
    conn.commit()
    conn.close()
    remove(I)

##############################################################################
            
def reaminate(seq1,seq2,conv1,conv2):
    I=range(len(conv1))

    # align original sequences by BLAST of fully-converted sequences
    for i in I:
        if conv1[i]=='-': seq1[i:i]=['-']
        if conv2[i]=='-': seq2[i:i]=['-']    
    assert len(conv1)==len(seq1)==len(conv2)==len(seq2)

    # clone seq2 as base for final output
    final=list(seq2)

    # resolve indels
    for i in I:
        if seq1[i]=='C'!=seq2[i]:

            # deletion
            if seq2[i]=='-':
                if 'C' in seq2[:i+2][-3:]: final[i]='C'
                elif 'T' in seq2[:i+2][-3:]: final[i]='T'
                else: final[i]='N'
            elif seq2[i]=='T':

                # insertion
                if 'C' in seq2[:i+2][-3:] and '-' in seq1[:i+2][-3:]:
                    final[i]='C'

                # neither
                else: final[i]='t'

    # remove insertions
    return Seq(''.join([final[i] for i in I if seq1[i]!='-']))

##############################################################################

def run_BLAST(query,subject,XML,mask):

    # incorporate lowercase masking
    if mask:
        call([c2bm,'-in',subject,'-out',suffix(subject,'asnb'),
              '-masking_algorithm','repeat',
              '-masking_options','"repeatmasker,default"',
              '-outfmt','maskinfo_asn1_bin','-parse_seqids'])
        call([mbdb,'-dbtype','nucl','-parse_seqids','-in',subject,
              '-mask_data',suffix(subject,'asnb')])

    # makeblastdb
    else: call([mbdb,'-dbtype','nucl','-parse_seqids','-in',subject])

    # don't know why but must use shell=True
    call(' '.join(

        # change complexity filter for query sequences
        [bn,'-task','megablast','-dust','"50 64 1"']

        # include masking information if specified
        +['-db_soft_mask','40']*mask

        # take advantage of supercomputing cluster
        +['-num_threads','8']

        # want XML output
        +['-outfmt','5','-query',query,'-db',subject,'-out',XML]),shell=True)
            
##############################################################################

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='')
    parser.add_argument('seqs',help='FASTA-format file of read sequences')
    parser.add_argument('refs',
                        help='FASTA-format file of reference sequences')
    parser.add_argument('-Save',help='output file name')
    parser.add_argument('-mask',action='store_true',
                        help='enable lowercase masking in references')
    parser.add_argument('-pair',type=str,
                   help='2nd FASTA-format file for paired-end sequencing')
    parser.add_argument('-update',action='store_true')
    parser.add_argument('-1','--a1_or_b1',action='append_const',const='a',
                        dest='strand',help='convert BS-read G to A, i.e.'
                        'seq. primer=a1 or b1')
    parser.add_argument('-2', '--a2_or_b2', action='append_const', const='b',
                   dest='strand',help='convert BS-read C to T, i.e.'
                   'seq. primer=a2 or b2')
    args=parser.parse_args()

    # if strand unspecified, try both
    if not args.strand: args.strand='ab'
    else: args.strand=''.join(sorted(args.strand)).lower()

    # label output
    if not args.Save:
        args.Save=path.join(path.dirname(args.seqs),'.'.join([
            path.split(args.seqs)[-1].split('.')[0][:10],
            path.split(args.refs)[-1].split('.')[0][:10],'']))
    for arg,val in sorted(vars(args).items()): print arg,val
    print

    put_data(args.seqs,args.refs,args.Save,args.mask,args.strand,args.update)
    
    if args.pair:
        args.strand={'a': 'b','b': 'a','ab': 'ab'}[args.strand]
        put_data(args.pair,args.refs,args.Save,args.mask,args.strand,True)
