#!/usr/bin/python2.7
# (c) 2013, Russell Darst, University of Florida

from Bio import Cluster
from meTools import *
from os import listdir,path
from PIL import Image
from time import time
import ConfigParser

##############################################################################

# dependencies

config = ConfigParser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)), 'gouache.cfg'))

##############################################################################

def clust(source):
    with open(source,'r') as handle: data=Cluster.read(handle)
    tree=data.treecluster()
    tree.scale()
    data.save(path.splitext(source)[0],tree)

##############################################################################

def color(source):
    code=eval(config.get('KEY','colors'))
    with open(source,'r') as handle: lines=handle.readlines()
    lines=[line.strip('\n').strip('\r').split('\t') for line in lines]
    canvas=Image.new('RGB',(len(lines[0])-4,len(lines)-2),'white')
    for y,line in enumerate(lines[2:]):
        for x,item in enumerate(line[4:]):
            try: canvas.putpixel((x,y),code[float(item)])
            except: continue
    canvas.save(path.splitext(source)[0]+'.png')
    
##############################################################################

def cook(save,*data,**opts):

    # read default values
    cols,incr,least,defs,patch,window=[
        opts.get(x,eval(config.get('args',x)))
        for x in ('cols','incr','least','defs','patch','window')]
    if type(cols)!=list: cols=[cols]
    N=[zip(*defs)[0].index(i)+2 for i in cols]
    save=path.splitext(save)[0]+'.cdt'

    # detect input type
    if not data: print 'No data entered.'
    if not all([type(i)==list for i in data]):
        data=read_from_file(cols,defs,*data)
    
    if not window:
        window=max([len(i.strip()) for row in data for i in row])
    window=range(-window/2,window/2+1)
            
    # header
    results=[['read']+[i for j in cols for i in ['-']*4+[j]+window]]

    for row in data:
        D=[['#- +*'.find(i)-2 for i in row[j]] for j in N]
        for n,q in enumerate([quilt(row[i])[1:-1] for i in N]):
            if not patch: continue
            for p in q:
                
                # methylated: solid yellow
                if p['V']==1:
                    if n==opts.get('mask',-1): C=.5
                    else: C=1

                # unbounded unmethylated: black
                elif not isinstance(p['D']-p['A'],int) or\
                     n==opts.get('mask',-1): C=-1

                # average of AD and BC => purple, blue, or green
                else:
                    m=(p['C']+p['D']-p['A']-p['B'])/2
                    C=(max(0,cmp(incr[0],m))-2-max(0,cmp(m,incr[1])))/4.
                for i in range(p['B'],p['C']+1): D[n][i]=C

        # offset from reference center and mask
        D = [{i - len(d)/2: j for i,j in enumerate(d)} for d in D]
        E = [[k for k,v in d.items() if v != 0] for d in D]
        E = zip(D, [(min(e + [1000]), max(e + [-1])) for e in E])
        D = [{k: v for k,v in d.items() if e[0] <= k <= e[1]} for d,e in E]
        D = [[i.get(j,'') for j in window] for i in D]

        # optional cutoff for data content
        if all([d.count('')*100./(len(d))+least>100 for d in D]): continue

        # add to results
        results+=[['|'.join(row[:2])]+[
            i for d in D for i in ['',10,10,10,'',]+d]]

    # write only if there was data
    if len(results)<3: return
    with open(save,'w') as handle:
        handle.write('\n'.join(['\t'.join(
            [str(i) for i in j]) for j in results]))

    clust(save)
    color(save)

##############################################################################

def read_from_file(cols,defs,*names):
    rows=[]
    for source in names:
        try: MAPit=read_FASTA(source,**dict(defs))
        except: continue
        for i,seq in MAPit.seqs.items():
            row=[MAPit.__rules__[j](seq) for j in ['locus','seqID']+cols]

            # bsDraw reverse strand output
            if path.split(source)[1][:2]=='b-':
                row[2:]=[cell[::-1] for cell in row[2:]]
            rows+=[row]
    return rows
