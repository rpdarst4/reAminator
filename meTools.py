#!/usr/bin/python2.7
# (c) 2013, Russell Darst, University of Florida

from Bio import SeqIO
from numpy import NaN

##############################################################################

def quilt(pattern):

    # no data
    if not pattern.strip():
        return []

    # imaginary white space
    P=pattern+' '

    # get inner bounds B and C, value V
    Q=[{'B': i, 'V': '#*'.find(j),'C': i + P[i:].find(' ') - 1}
       for i,j in enumerate(P) if ' ' == P[i-1] != j]
    Q = [{'A': NaN, 'B': NaN, 'C': NaN, 'D': NaN, 'V': -1}
         ]+Q+[{'A': NaN, 'B': NaN, 'C': NaN, 'D': NaN, 'V': -1}]

    # get outer bounds A and D
    N=range(len(Q))
    for i,j in zip(N[1:],N[:-1]):
        Q[i]['A']=Q[i-1]['C']
        Q[j]['D']=Q[j+1]['B']

    return Q
    
##############################################################################

def read_FASTA(whichFile,**kwargs):
    with open(whichFile,'r') as handle:
        seqs=list(SeqIO.parse(handle,'fasta'))
    return meTable(seqs,**kwargs)

##############################################################################

def seek(obj, col=-1, keep='True', drop='False'):
    #,window=slice(None,None,None)

    """
    yields features for each sequence meeting criterion
    for each patch, produces a DICT:
        V, "value", 1 for methylated, 0 for unmethylated
        AD, distance between outer bounding sites
        BC, distance between inner bounding sites
    """

    # restrict input to allowed values
    keep = eval('lambda v, ad, bc: ' + keep)
    drop = eval('lambda v, ad, bc: ' + drop)
    
    if type(obj) is str: obj = read(obj)
    if type(col) is str: col = obj.__head__.find(col)

    """
    # range in which to seek patches
    if type(obj) is list: span=range(len(obj[0][col]))
    else: span=range(len(obj.ref))[window]
    """

    # feed V, AD, BC values to functions
    test = lambda f, x: f (x['V'], x['D']-x['A']-1, x['C']-x['B']+1)
    for seq in obj:
        patt = quilt(seq[col])

        # avoid errors
        if not len(patt): continue

        # last patch is "q"
        q = patt[0]

        # next patch is "p"
        for p in patt[1:]:

            # merge next patch "p" with last patch "q"
            if test(drop, p) or p['V'] == q['V']:
                q['C'], q['D'] = p['C'], p['D']
                continue

            # cannot merge next patch "p";
            # report last patch "q" if in window of interest
            if test(keep, q):
                # and (q['A']<=span[-1] or q['B']<= span[-1])
                # and (q['C']>=span[ 0] or q['D']>= span[ 0]):

                # to count sites
                ends = slice(q['B'], q['C'] + 1)
                
                yield (

                    # locus, read
                    seq[0],seq[1],

                    # patch ends
                    q['A'],q['B'],q['C'],q['D'],

                    # number sites unmethylated and methylated
                    seq[col][ends].count('#'),seq[col][ends].count('*'))

            # replace last patch "q" with next patch "p"
            q = dict(p)

##############################################################################

class meTable():
    ref,seqs,unambig=None,[],True

    def __init__(self,contig,BS=95,unambig=True,**kwargs):
        self.BS = BS / 100.
        self.__head__    = []
        self.__offsets__ = {}
        self.__sites__   = {}        
        self.ref=contig[0]
        self.seqs={seq.id: seq for seq in contig[1:]}

        # DICT of LAMBDA functions
        self.__rules__={'locus': lambda x: self.ref.id,
                        'seqID': lambda x: x.id}

        # default sites
        if not kwargs: kwargs=dict(CG=1,GC=2)
        for k,v in kwargs.items():

            try: assert type(v) is int and v!=0
            except: IOError('Site parameters must be in format GC=2, etc.')

            # add site to Header
            self.__head__+=[k]

            # get list of instances of site
            self.__sites__[k]=self.__match__(k,v)
            self.__offsets__[k]=v

            # make new rule for methylation site k
            self.__rules__[k]=lambda x,y=k: self.__score__(x,y)
        self.__head__=['locus','seqID']+sorted(self.__head__)
        self.__checkBS__()

        # remove overlapping sites from each list of positions
        if unambig:
            self.__sites__={k: self.__pare__(k,self.__sites__[k])
                            for k in kwargs.keys()}

    def __checkBS__(self):
        B=self.BS/(1-self.BS)
        C=self.__pare__('C',self.__match__('C',1))
        D=lambda x,y: sum([x[c].upper()==y for c in C])        
        self.seqs={i: j for i,j in self.seqs.items() if D(j,'T')>=B*D(j,'C')}

    def __getitem__(self,i):
        if type(i)==str:
            return [self.__rules__[j](self.seqs[i]) for j in self.__head__]
        else:
            return [self.__rules__[j](self.seqs.values()[i])
                    for j in self.__head__]

    def __len__(self):
        return len(self.seqs)

    # get bp positions of all instances of site M (methylated at bp N)
    def __match__ (self,M,N):
        ref=self.ref.seq.upper()
        return [n+N-1 for n in range(len(ref)) if ref.find(M,n)==n]

    # list all bp position from list I that do not overlap with
    # a methylation site other than K
    def __pare__(self,K,I):
        return sorted(set(I)-set(
            [j for k in self.__head__[2:] if k!=K for j in self.__sites__[k]]))

    def __score__(self,read,site):
        check=lambda x: (x=='C')-(x=='T')
        ref_len,seq=len(self.ref.seq),read.seq
        I=self.__sites__[site]
        J=[check(seq[i].upper()) for i in I]
        K=[i for i,j in enumerate(J) if j!=0]
        if len(K)==0: return ' '*ref_len

        # bp before first site left blank
        L=[' '*I[K[0]]]

        # bp between sites colored '-', '+', or ' '
        for i,j in  zip(K[:-1],K[1:]):
            L+=[' *#'[J[i]]]+[' +-'[(J[i]+J[j])/2]]*(I[j]-I[i]-1)

        # bp after last site left blank
        L+=[' *#'[J[-1]],' '*(ref_len-I[K[-1]]-1)]
        return ''.join(L)

    def screen(self,i):
        n,seq=0,self.seqs[i]
        while n<len(self.ref):            
            print 'ref     '+self.ref.seq[n:][:70]
            print 'read    '+seq.seq[n:][:70]
            for site in self.__head__[2:]:
                print site+' '*(8-len(site)) \
                      +self.__rules__[site](seq)[n:][:70]
            n+=78
            print
                
    def write(self,OUT,header=False):
        if type(OUT)==str: OUT=open(OUT,'w')
        OUT.write ('\n'.join(['\t'.join(self.__head__)]*header + [
            '\t'.join([str(i) for i in j]) for j in self[:]]))
