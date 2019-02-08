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
import multiprocessing, codecs
from argparse import ArgumentParser

# Set encoder for DNA string
class GenomeCodec(codecs.Codec):
    def encode(self, input_, errors='strict'):
        return codecs.charmap_encode(input_, errors, encoding_table)

    def decode(self, input_, errors='strict'):
        return codecs.charmap_decode(input_, errors, decoding_table)

def lookup(name):
    if name != 'dna':
        return None
    return codecs.CodecInfo(
        name = 'dna',
        encode=GenomeCodec().encode,
        decode=GenomeCodec().decode)

decoding_table = ('A' 'C' 'G' 'T')
encoding_table = codecs.charmap_build(decoding_table)

codecs.register(lookup)

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

        for itr in range(len(rseq)):
            splt = rseq[itr].split()
            self.chr_seq.update({ splt[0] : splt[-1] })
            self.chr_len.update({ splt[0] : len(splt[-1]) })

        # del rseq
        # del splt

    def _ncount(plist):
        plist = sorted(plist)
        result = {}

        ncnt, start = 1, plist[0]
        for itr in range(len(plist)):
            ncurr = plist[itr]
            nnext = plits[itr+1]
            if ncurr + 1 == nnext:
                ncnt += 1
            else:
                result.update({ start : ncnt })
                ncnt, start = 1, plist[itr+1]
    
    def _compress(self): # Remove N and compress size of sequence
        for chrom, seq in self.chr_seq.items():
            npos = [i for i, ktr in enumerate(seq) if ktr == 'N']
            npos = _ncount(npos)

            self.nidx.update({ chrom : npos })
            self.chr_seq[chrom] = seq.replace('N', '').encode('dna')

    def extract(self):
        _load(self)

    """
            with open(self.fn, 'r') as fp:
                chrom, Nstr = '', False
                n_pos = 0

                for line in fp:
                    line = line.strip()
                    line = line.split()[0]

                    if line.startswith('>'):
                        chrom, Nstr = line.replace('>', ''), False
                        n_pos = 0
                        if 'chr' in chrom:
                            chrom.replace('chr', '')

                        print('Extracting chr %s' % chrom)
                        self.chr_seq.update({ chrom : '' }) #.encode('dna') })
                        self.chr_len.update({ chrom : 0 })
                        self.n_idx.update({ chrom : {} })
                        continue

                    begin, end = line.find('N'), line.rfind('N') + 1
                    assert (begin == 0 or end == len(line)) or\
                             (begin == -1)

                    if not Nstr:
                        n_pos = self.chr_len[chrom] + len(line[:begin]) + 1

                    if begin == -1:
                        self.chr_seq[chrom] += line # .encode('dna')
                    else:
                        if n_pos not in self.n_idx[chrom]: 
                            self.n_idx[chrom].update({ n_pos : 0 })
                        
                        self.n_idx[chrom][n_pos] += len(line[begin:end])
                        self.chr_seq[chrom] += line[:begin] #.encode('dna')
                        self.chr_seq[chrom] += line[end:] #.encode('dna')

                    Nstr = True if end == len(line) else False
                    self.chr_len[chrom] += len(line)
    """

"""
"""
if __name__ == '__main__':
    # ref = Genome('genome.fa')
    # ref.load()

    ref = open('genome.fa', 'r').read()
    ref = ref.strip('\n').split('>')[1:]

    input('Press any key ....')
    print("Done")

