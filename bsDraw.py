#!/usr/bin/python2.7
# (c) 2013, Russell Darst, University of Florida

"""
* extracts sequences from bsBLAST output table
* output is folder of FASTA-formatted contigs and PNG files
* assumes CG and GC methylation
* direct Python commands needed for other methylation patterns
"""

from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from meTools import *
from os import mkdir,path
from time import time
import argparse as ap
import ConfigParser
import sqlite3 as sql

##############################################################################

# dependencies

config = ConfigParser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)), 'reAminator.cfg'))
if config.get('paths','methylmapper'): import meMapper as meMapp
if config.get('exists','PIL').lower() == 'true': import omelet

##############################################################################

def snowflake(seqs, sites):
    U = {''.join([str(seq[c].upper()) for c in sites]): seq for seq in seqs}
    K = []
    for u in sorted(U.keys()):
        for n,k in enumerate(K):

            # k in u: replace k 
            if all([i in ('-', j) for i, j in zip(k, u)]): K[n] = u

            # mismatch: try next
            elif not all([j in ('-', i) for i, j in zip(k, u)]): continue

            # u in k: skip
            break

        # u not found: add to K
        else: K.append(u)
    return [U[k] for k in K]
            
##############################################################################

def unpack(source, dest=None,
           expts=[], loci=[], strands=[],
           methyl={'CG': 1, 'GC': 2},
           min_len='100', min_bs=95, cull_dupes=False):
    
    # choose how to evaluate length
    if '%' in min_len:
        len_test = lambda x, y: len(x) - x.count('-') >= int(
            min_len.strip('%')) / 100. * len(y)
    else:
        len_test = lambda x, y: len(x.replace('-', '')) >= int(min_len)

    # prepare target directory "dest"
    if not dest: dest = path.split(source)[0]
    if not path.exists(dest): mkdir(dest)

    # keep track of outcomes
    report = {}
    
    conn = sql.connect(source)
    curs = conn.cursor()

    # if terms omitted, do all
    T = {'expt': expts, 'locus': loci, 'strand': strands}
    for k in T.keys():
        if T[k]: continue
        curs.execute('SELECT {0} FROM records GROUP BY {0}'.format(k))
        T[k] = [j[0] for j in curs.fetchall()]

    comb = sorted([[x, y, z] for x in T['expt'] for y in T['locus'] \
                   for z in T['strand']], key = lambda x: x[1])

    for X,Y,Z in comb:

        # get reference
        curs.execute('SELECT locus, sequence from loci WHERE locus = ?', (Y,))
        i, j = curs.fetchone()
        ref  = SeqIO.SeqRecord(description = '', id = i, seq = Seq(j.upper()))

        # get sequences
        seqs = [i for i in curs.execute(
            'SELECT read, sequence FROM records WHERE '
            'expt = ? AND locus = ? AND strand = ?', (X,Y,Z))]
        if not seqs: continue
        N = [len(seqs)]

        # test length
        seqs = [SeqIO.SeqRecord(description = '', id = x, seq = Seq(y))
                for x, y in seqs if len_test(y, ref)]
        N.append(len(seqs))

        if Z == 'b':
            for seq in seqs: seq.seq = seq.seq.reverse_complement()
            ref.seq = ref.seq.reverse_complement()
        if not X: X = 'nd'

        # check deamination
        data = meTable([ref]+seqs, BS=min_bs, **methyl)
        seqs = data.seqs.values()
        N.append(len(seqs))

        # remove duplicates
        if cull_dupes:
            seqs = snowflake(seqs, data.__match__('C', 1))
            data = meTable([ref]+seqs, BS=min_bs, **methyl)
            N.append(len(seqs))

        # update report
        report[tuple([X, Y, Z])] = N[:1] + [
            N[n] - N[n+1] for n in range(len(N))[:-1]] + N[-1:]

        # include methylation count
        for j,k in sorted(methyl.items()):
            report[tuple([X, Y, Z])] += [
                sum([l[data.__head__.index(j)].count('*') for l in data])]
            report[tuple([X, Y, Z])] += [
                sum([l[data.__head__.index(j)].count('#') for l in data])]

        # write contig
        if not seqs: continue
        contig=path.join(dest,'-'.join([Z,X,Y]).replace(':','-')+'.fa')
        if path.exists(contig):
            old=list(SeqIO.parse(contig,'fasta'))
            if old:
                combined={seq.id: seq for seq in old[1:]}
                for seq in seqs: combined[seq.id]=seq
                seqs=combined.values()
        with open(contig,'w') as handle:
            SeqIO.write([ref]+seqs,handle,'fasta')
        yield contig

    # write report
    with open(path.join(dest,'report.tsv'),'w') as handle:
        for i in str(datetime.now()).split():
            handle.write('#\t'+i+'\n')
        handle.write('#\tSOURCE\t{}\n#\tDEST.\t{}\n#\tEXPTS.\t{}\n#\tLOCI\t'
                     '{}\n#\tSTRANDS\t{}\n#\tLENGTH\t{}\n#\tDEAM.\t{}%\n#\n'
                     .format(source, dest, expts, loci, strands, min_len,
                             min_bs))
        handle.write('#EXPT\tLOCUS\tSTRAND\tALIGN.\t<{}\t<{}%DEAM.\t{}PASSED'
                     .format(min_len, min_bs, 'duplicates\t' * cull_dupes)
                     + '\t' + '\t'.join([
                         '{}\t{}t{}'.format(k, k[:v][:-1], k[v:])
                         for k,v in sorted(methyl.items())]))
        handle.write('\n'.join([''] + ['\t'.join(list(XYZ) + [
            str(n) for n in report[XYZ]]) for XYZ in sorted(
                report.keys(), key=lambda x: x[1])]))
        
##############################################################################

if __name__ == '__main__':
    p = ap.ArgumentParser (description='')

    # input TSV file from bsBLAST
    p.add_argument ('alignments', default=[], nargs='*',
                    help='bsBlast output table file/s')

    # option to unpack to a chosen directory
    p.add_argument ('-dest', help='output directory')

    # options to extract only specific groups of sequences
    p.add_argument ('-codes', default=[], nargs='*',
                    help='barcode seqs (leave blank for all)')
    p.add_argument ('-refs', default=[], nargs='*',
                    help='ref. IDs (leave blank for all)')
    p.add_argument ('-strand', default='ab',
                    help='strands ("a" = C to T, "b" = G to A)')

    # filters
    p.add_argument ('-bisulfite', default=95, type=int,
                    help='minim. percent deaminated (0-100)')
    p.add_argument ('-length', default='100', 
                    help='minim. bp length (default) or percent')
    p.add_argument ('-uniques', action='store_true',
                    help='use only unique non-deamination patterns')

    # option to change methylation sites
    p.add_argument ('-Sites', type=str, help='enter alternate '
                    'methylation dictionary (e.g. CG=1,CC=1)\n'
                    'NOTE: not compatible with MethylMapper')

    # option to extract TSS position data (3' to 5' end of reference)
    p.add_argument ('-TSS', action='store_true',
                    help='read TSS distance after "+" symbol in ref. id'
                    ' (e.g. >YFG+400)')    

    # option to call MethylMapper
    p.add_argument ('-Weights', default=None,
                    help='ratio of CG to GC weight for MethylMapper'
                    ' (e.g. 50,50)')

    # option to call omelet
    p.add_argument ('-Omelet', action='store_true',
                    help='make omelet .png files')

    # option to adjust window length
    p.add_argument ('-window', default=600, type=int,
                    help='basepair width of omelet .png files '
                    '(or zero for automatic)')
    
    a = p.parse_args()

    if a.Sites: sites = [(k,int(v)) for k,v in [
        x.split('=') for x in a.Sites.split(',')]]
    else: sites = [('CG', 1), ('GC', 2)]
    a.strand = ''.join(sorted(a.strand)).lower()
    
    for arg, val in sorted(vars(a).items()): print arg, val
    print

    for alignment in a.alignments:
        for contig in unpack(
            alignment, a.dest, a.codes, a.refs, a.strand, dict(sites),
            a.length, a.bisulfite, a.uniques):
            if a.Weights or a.TSS:
                if not a.Weights: a.Weights = '50,50'
                meMapp.plot(contig,False,a.Sites,a.TSS,a.Weights)
            if a.Omelet:
                A,B=path.split(contig)
                if config.get('exists','PIL').lower()=='true':
                    omelet.cook(path.join(A,B[2:]),
                                path.join(A,'a'+B[1:]),
                                path.join(A,'b'+B[1:]),
                                methyl=sites, window=a.window)
                else: raise IOError('Called -Omelet but PIL is False')
