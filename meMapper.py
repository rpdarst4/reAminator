#!/usr/bin/python2.7
# (c) 2013, Russell Darst, University of Florida

from subprocess import call
import argparse as ap
import ConfigParser
import os.path as path

##############################################################################

config=ConfigParser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)),'reAminator.cfg'))
mm=config.get('paths','methylmapper')

##############################################################################

def plot(source,GCG=False,site=None,TSS=False,weights=None):
    try: TSS=TSS*int(path.splitext(source)[0].split('+')[-1])
    except: TSS=False
    if not site: site='both'
    if not weights: weights='50,50'
    call([mm,'-in',source,'-site',site,'-weight',weights,'-gap','0',
          '-plot',path.splitext(source)[0]+(GCG*'wGCG')+'.'+weights+
          '.png']+bool(TSS)*['-tss',str(TSS)]+GCG*['-gcg','False'])

##############################################################################

if __name__=='__main__':
    p=ap.ArgumentParser(description='')
    p.add_argument('seqs',help='FASTA-format file of sequences')
    p.add_argument('-gcg',action='store_true',help='count GCG sites '
                   '(default is to discard GCG information)')
    p.add_argument(
        '-Site',default='both',
        help='enter site "CG" or "GC" or leave blank for both')
    p.add_argument(
        '-TSS',action='store_true',
        help='read TSS distance after "+" symbol in ref. id (e.g. >YFG+400)')
    p.add_argument(
        '-Weights',default='50,50',
        help='enter CG,GC weighting e.g. "70,30"')
    args=p.parse_args()
    cmdline='{0} -in {1} -plot {2} -site {3} -weight {4} -gap 0'.format(
        mm,args.seqs,path.splitext(args.seqs)[0]+((not args.gcg)*'wGCG')+
        '.png',args.Site,args.Weights)+' -gcg'*(not args.gcg)
    print cmdline
    call(cmdline.split())
