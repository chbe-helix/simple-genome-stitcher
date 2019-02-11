#!/usr/bin/env python3

# 
# Copyright 2019, Christopher Bennett
# 
# This is a prototype script to align two linear sequences (genomes)
# 
# fast-genome-aligner is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# fast-genome-aligner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# See <http://www.gnu.org/licenses/> for a copy of the GNU General Public License
# 

import sys, os, subprocess, re
import multiprocessing
from argparse import ArgumentParser

# Set encoder for DNA string
def dna_encoder(seq):
    encode = {'A' : '0',
              'C' : '1',
              'G' : '2',
              'T' : '3'}
    
    table = str.maketrans(encode)
    seq = seq.translate(table)

    seq = [int(nt) for nt in seq]
    return(seq)
    
    
def ncount(plist):
    plist = sorted(plist)
    result = {}

    ncnt, start = 1, plist[0]
    for itr in range(len(plist)):
        if itr+1 == len(plist):
            result.update({ start : ncnt })
            break

        ncurr = plist[itr]
        nnext = plist[itr+1]
        if ncurr + 1 == nnext:
            ncnt += 1
        else:
            result.update({ start : ncnt })
            ncnt, start = 1, plist[itr+1]

    return(result)

# Load DNA string from fasta
class Genome:
    def __init__(self, fasta_fn):
        self.fn = fasta_fn
        self.chr_seq = {}
        self.chr_len = {}
        self.nidx = {}

    def _load(self):
        try:
            rseq = open(self.fn, 'r').read()

        except Exception as ex:
            print(ex)

        rseq = rseq.strip('\n').split('>')[1:]
        while rseq:
            nlpos = rseq[0].find('\n')
            chrom = rseq[0][:nlpos].split()[0]
            if 'chr' in chrom:
                chrom = chrom.replace('chr','')

            self.chr_seq.update({ chrom : rseq[0][nlpos:].replace('\n','') })
            self.chr_len.update({ chrom : len(rseq[0][nlpos:].replace('\n','')) })
            rseq.pop(0)

        # del rseq
        # del splt

    def _compress(self): # Remove N and compress size of sequence
        for chrom, seq in self.chr_seq.items():
            npos = [i.start() for i in re.finditer('N', seq)]

            # npos = [i for i, ktr in enumerate(seq) if ktr == 'N']
            npos = ncount(npos)

            self.nidx.update({ chrom : npos })
            self.chr_seq[chrom] = dna_encoder(seq.replace('N', ''))

    def extract(self):
        self._load()
        # self._compress()


def binary_compare(seq1, seq2, p1, p2, l):
    if l > 1000:
        while True:
            if seq1[p1:p1+l] == seq2[p2:p2+l]:
                p1 += l
                p2 += l
                continue
            else:
                p1, p2 = binary_compare(seq1, seq2, p1, p2, l/2)
    else:
        while seq1[p1] == seq2[p2]:
            p1 += 1
            p2 += 1
        return(p1, p2)
        
def compare(seq1, seq2, l):
    pidx = []
    p1, p2 = 0, 0
    while True:
        if len(seq1[p1:]) == len(seq2[p2:]):
            break
        if p1+l > len(seq1) or p2+l > len(seq2):
            l = int(l/2)

        if seq1[p1:p1+l] == seq2[p2:p2+l]:
            p1 += l
            p2 += l
        else:
            p1, p2 = binary_compare(seq1, seq2, p1, p2, l/2)
            print(p1, p2)
            i = 1
            while i < 15:
                if seq1[p1+i:p1+i+50] == seq2[p2:p2+50]:
                    p1 += i
                    pidx.append([p1, p2])
                    break
                elif seq1[p1:p1+50] == seq2[p2+i:p2+i+50]:
                    p2 += i
                    pidx.append([p1, p2])
                    break
                elif seq1[p1+i:p1+i+50] == seq2[p2+i:p2+i+50]:
                    p1 += i
                    p2 += i
                    break
                i += 1

    return(pidx)

"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='Simple Genome Stitcher')
    parser.add_argument('--g1',
                        dest='genome_1',
                        type = str,
                        required = True,
                        help='First genome of interest')
    parser.add_argument('--g2',
                        dest='genome_2',
                        type = str,
                        required = True,
                        help='Second genome of interest')

    args = parser.parse_args()

    g1 = Genome(args.genome_1)
    g2 = Genome(args.genome_2)

    g1.extract()
    g2.extract()

    for chrom, seq in g1.chr_seq.items():
        if chrom not in g2.chr_seq:
            print('%s not in %s skipping' % (chrom, g2))
            continue
        pidx = compare(seq, g2.chr_seq[chrom], 100)

        print(chrom, pidx)
    
    print(g1.chr_len)
    print(g2.chr_len)

