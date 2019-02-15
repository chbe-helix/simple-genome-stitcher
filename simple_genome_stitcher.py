#!/usr/bin/env python3

# 
# Copyright 2019, Christopher Bennett
# 
# This is a prototype script to align two linear sequences (genomes)
# 
# simple-genome-stitcher is free software: you can redistribute it and/or modify
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

###
### Genome Functions
###
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
        self.fai = {}

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

    def _compress(self): # Remove N and compress size of sequence
        for chrom, seq in self.chr_seq.items():
            npos = [i.start() for i in re.finditer('N', seq)]
            npos = ncount(npos)

            self.nidx.update({ chrom : npos })
            self.chr_seq[chrom] = seq.replace('N', '') # Ideally ACGT would be represented here in 2 bit form

    def load_idx(self):
        gidx = self.fn + '.fai'
        if not os.path.exists(gidx):
            print ("Genome index does not exist")

        try:
            with open(gidx, 'r') as f:
                for line in f:
                    line = line.strip().split()
                    chrom = line[0].replace('chr', '')
                    self.fai.update({ chrom : line [1:] })

        except Exception as ex:
            print(ex)



    def extract(self):
        self._load()
        # self._compress() # Not yet working fast

###
### Stitcher Functions
###
class SW_aligner:
    def __init__(self, seq1, seq2, mw, mmw, gw, anc, ancw, s = -1):
        if len(seq1) < len(seq2): 
            tmp = seq1
            seq1, seq2 = seq2, tmp
            del tmp

        self.seq1       = 'Z' + seq1 # make both seq one base
        self.seq2       = 'Z' + seq2
        self.weights    = [mw, mmw, gw, anc, ancw, s]# match score, mismatch score, gap score, anchor len, anchor score
        self.matrix     = None
        self.pos        = None
        self.alignments = None
        self.align_str  = None

    def _SWmatrix(self): # sequences, match score, mismatch score, gap score
        # Make SW matrix
        seq1, seq2 = self.seq1, self.seq2
        m, mm, g, a, aw, s = self.weights
 
        cols, rows = len(seq1), len(seq2)
        score_matrix = [[0 for col in range(cols)] for row in range(rows)]

        # Fill the scoring matrix.
        max_score = 0
        max_pos   = None    # The row and columbn of the highest score in matrix.
        for i in range(1, rows):
            for j in range(1, cols):
                ms = m if seq2[i] == seq1[j] else mm

                ds = score_matrix[i-1][j-1] + ms
                us = score_matrix[i-1][j] + g
                ls = score_matrix[i][j-1] + g

                score = max(0, ds, us, ls)
                if score > max_score:
                    max_score = score
                    max_pos   = (i, j)

                if (i == j) and (i <= a and j <= a):
                    score_matrix[i][j] = score + aw
                elif (s != -1) and (i == j) and (i > rows-s and j > cols-s):
                    score_matrix[i][j] = score + aw
                else:
                    score_matrix[i][j] = score

        assert max_pos is not None, 'the x, y position with the highest score was not found'
        self.matrix, self.pos = score_matrix, max_pos

    def _traverse(self, score_matrix, i, j):
        diag = score_matrix[i - 1][j - 1]
        up   = score_matrix[i - 1][j]
        left = score_matrix[i][j - 1]
        if diag >= up and diag >= left:     # Tie goes to the DIAG move.
            return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
        elif up > diag and up >= left:      # Tie goes to UP move.
            return 2 if up != 0 else 0      # UP move or end.
        elif left > diag and left > up:
            return 3 if left != 0 else 0    # LEFT move or end.
        else:
            # Execution should not reach here.
            raise ValueError('invalid move during traceback')

    def _traceback(self):
        END, DIAG, UP, LEFT = range(4)
        aligned_seq1 = []
        aligned_seq2 = []
        i, j         = self.pos
        move         = self._traverse(self.matrix, i, j)
        while move != END:
            if move == DIAG:
                aligned_seq1.append(self.seq1[j])
                aligned_seq2.append(self.seq2[i])
                i -= 1
                j -= 1
            elif move == UP:
                aligned_seq1.append('-')
                aligned_seq2.append(self.seq2[i])
                i -= 1
            else:
                aligned_seq1.append(self.seq1[j])
                aligned_seq2.append('-')
                j -= 1

            move = self._traverse(self.matrix, i, j)

        aligned_seq1.append(self.seq1[j])
        aligned_seq2.append(self.seq2[i])

        self.alignments = [''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))]

    def _alignment_string(self):
        idents, gaps, mismatches = 0, 0, 0
        aligned_seq1, aligned_seq2 = self.alignments
        alignment_string = []
        for base1, base2 in zip(aligned_seq1, aligned_seq2):
            if base1 == base2:
                alignment_string.append('|')
                idents += 1
            elif '-' in (base1, base2):
                alignment_string.append('~')
                gaps += 1
            else:
                alignment_string.append(':')
                mismatches += 1

        self.align_str = ''.join(alignment_string), idents, gaps, mismatches

    def dyn_align(self):
        self._SWmatrix()
        self._traceback()
        self._alignment_string()

    def get_align_str(self):
        return self.align_str

    def get_seq_align(self):
        return self.alignments

    def get_matrix(self):
        return self.matrix

def kmer_compare(seq1, seq2, p1, p2, k):
    kmers1, kmers2 = [], []
    p1a, p2a = 0, 0 # Additions to 1 and 2

    i = 0
    while i < len(seq1) - k:
        kmer1, kmer2 = seq1[i:i+k], seq2[i:i+k]
        
        match = False
        if kmer1 == kmer2:
            p1a += i
            p2a += i
            match = True
        elif kmer1 in kmers2:
            p1a += i
            p2a += kmers2.index(kmer1)
            match = True
        elif kmer2 in kmers1:
            p1a += kmers1.index(kmer2)
            p2a += i
            match = True

        if match:
            return(True, p1a, p2a)

        i += 1
        kmers1.append(kmer1)
        kmers2.append(kmer2)

    print('No match in kmer size %s in window %s' % (k, len(seq1)))
    return(False, p1a, p1a)

def binary_compare(seq1, seq2, p1, p2, l):
    p1, p2, l = int(p1), int(p2), int(l)
    # print(p1, p2, l)
    if l > 1000:
        while seq1[p1:p1+l] == seq2[p2:p2+l]:
            p1 += l
            p2 += l

        p1, p2 = binary_compare(seq1, seq2, p1, p2, l/2)

    else:
        while seq1[p1] == seq2[p2]:
            p1 += 1
            p2 += 1
    
    return(p1, p2)

def compare(seq1, seq2, l, s, isize):
    pidx = []
    p1, p2 = 0, 0
    while True:
        if len(seq1[p1:]) == len(seq2[p2:]):
            break
        if p1+l+1 > len(seq1) or p2+l+1 > len(seq2):
            l1 = len(seq1[p1:p1+l])
            l2 = len(seq2[p2:p2+l])
            l =  int(l1/2) if l1 < l2 else int(l2/2)

        # print(p1, p2, l)
        if seq1[p1:p1+l] == seq2[p2:p2+l]:
            p1 += l
            p2 += l
        else:
            p1, p2 = binary_compare(seq1, seq2, p1, p2, (l/2))
 
            i = 1
            while i < isize:
                if seq1[p1+i:p1+i+s] == seq2[p2+i:p2+i+s]:
                    p1 += i
                    p2 += i
                    break                
                elif seq1[p1+i:p1+i+s] == seq2[p2:p2+s]:
                    p1 += i
                    pidx.append([p1, p2])
                    break
                elif seq1[p1:p1+s] == seq2[p2+i:p2+i+s]:
                    p2 += i
                    pidx.append([p1, p2])
                    break
                i += 1

            if i == isize:
                m, mm, g, aw = 3, -1, -1, 10 # SW Weights: match, missmatch, gap, anchor
                b, e = 10, 30               # Bases to take at Begin and End of p1 and p2 site
                k = 50                   # K-mer size
                aseq1, aseq2 = '', ''
                while e < 200:
                    aln = SW_aligner(seq2[p2-b:p2+(e*2)], seq1[p1-b:p1+e], m, mm, g, b, aw)
                    aln.dyn_align()
                    aln_str = aln.get_align_str()[0]
                    if aln_str.endswith('|' * s):
                        aseq2, aseq1 = aln.get_seq_align()
                        break
                    e += s

                if not aseq1 or not aseq2:
                    print('No matches: trying k-mer matching')
                    ws = 5000 # window size for kmer comparison
                    found = False
                    while ws < 100000 and not found: 
                        found, p1a, p2a = kmer_compare(seq1[p1:p1+ws], seq2[p2:p2+ws], p1, p2, k)
                        ws += 5000

                    if p1a == 0 or p2a == 0:
                        print('K-mer matching found insert. Setting Genome_1 to %d and Genome_2 to %d' % (p1+p1a, p2+p2a))
                        p1 += p1a
                        p2 += p2a
                        pidx.append([p1, p2])
                        continue

                    if found:
                        print('K-mer match found at distances of %d in Genome_1 and %d in Genome_2. Resuming indexing' % (p1a, p2a))
                        gg1 = True if p1a > p2a else False
                        p1ckr, p2ckr = p1 + p1a, p2 + p2a
                        p1 -= 1; p2 -= 1
                        p1a += 1; p2a += 1

                        if gg1:
                            ratio = int(p1a/p2a)
                            remain = p1a%p2a
                            while p2a > 0:
                                p2 += 1
                                p2a -= 1
                                if p2a == 0:
                                   p1 += (ratio+remain)
                                else:
                                    p1 += (ratio)
                                    
                                pidx.append([p1, p2])
                                
                        else:
                            ratio = int(p2a/p1a)
                            remain = p2a%p1a
                            while p1a > 0:
                                p1 += 1
                                p1a -= 1

                                if p1a == 0:
                                   p2 += (ratio+remain)
                                else:
                                    p2 += (ratio)
                                pidx.append([p1, p2])
                                        
                        assert p1ckr == p1 and p2ckr == p2, '%s : %d != %d, %d != %d' % (gg1, p1ckr, p1, p2ckr, p2)
                        continue

                    else:
                        print('No K-mer match found. Reevaluate options')
                        exit(1)

                r2 = [p2-b, p2+(e*2)]
                r1 = [p1-b, p1+e]
                print('Aligning Regions:%s and %s in seq2 and seq 1\n%s\n%s\n%s' % (r2, r1, aseq2, aln_str, aseq1))

                if aln_str.startswith('|' * b):
                    p1 -= b
                    p2 -= b
                else:
                    print('Max insert size reached and no definitive match: p1: %s, p2:%s' % (p1, p2))
                    exit(1)                    

                for ntidx in range(len(aln_str)):
                    c = aln_str[ntidx]
                    nt1, nt2 = aseq1[ntidx], aseq2[ntidx]
                    if c != '~':
                        p1 += 1
                        p2 += 1
                    else:
                        if nt1 == '-':
                            p2 += 1
                        if nt2 == '-':
                            p1 += 1
                        pidx.append([p1, p2])

    return(pidx)

def check(g1, g2, chrom_idx):
    g1.load_idx()
    g2.load_idx()

    passed = True
    for chrom, ls in g1.fai.items():
        if chrom not in g2.chr_len:
            print('Skipping %s. Not in Genome_2' % chrom)
            continue
        ls1 = int(ls[0])
        ls2 = int(g2.fai[chrom][0])
        kdiff = abs(ls1-ls2)

        if ls1 != int(g1.chr_len[chrom]):
            print('Chr %s Failed length test in genome_1' % chrom)
            passed = False
        if ls2 != int(g2.chr_len[chrom]):
            print('Chr %s Failed length test in genome_2' % chrom)
            passed = False

        if not chrom_idx[chrom]:
            print('Chr %s NO CHANGE' % chrom)
            continue

        indeces = chrom_idx[chrom][-1]
        diff = abs(indeces[0] - indeces[1])

        if kdiff != diff:
            print('Chr %s FAILED: index conversion of know diff %d and calculated diff %d' % (chrom, kdiff, diff))
            passed = False
        else: 
            print('Chr %s PASSED' % chrom)

    return(passed)


def build_idx(g1, g2, chrom_idx):
    g1name = g1.fn.split('.')[0]
    g2name = g2.fn.split('.')[0]
    with open('%s-%s-coord.idx' % (g1name, g2name), 'w') as fo:
        fo.write('## Coordinant index for %s vs %s\n' % (g1name, g2name))
        fo.write('#ORD\t%s\t%s\n' % (g1name, g2name))
        for chrom, index in chrom_idx.items():
            entry = ''
            entry += chrom + '\t'
            entry += g1.fai[chrom][0] + '\t'
            entry += g2.fai[chrom][0] + '\t'
            for i in index:
                entry += '%d:%d,' % (i[0], i[1])
            entry = entry[:-1] + '\n'
            fo.write(entry)

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
    parser.add_argument('--window-size',
                        dest='window_size',
                        type = int,
                        default = 1000000,
                        help='Size of window to binary search (Default: 1000000)')
    parser.add_argument('--seed-size',
                        dest='seed_size',
                        type = int,
                        default = 10,
                        help='Size of seed to match (Default: 10')
    parser.add_argument('--insert-size',
                        dest='insert_size',
                        type = int,
                        default = 250,
                        help='Base size of insert allowed before alternative methods (Default: 250)')

    args = parser.parse_args()

    g1 = Genome(args.genome_1)
    g2 = Genome(args.genome_2)

    g1.extract()
    g2.extract()

    chrom_idx = {}
    for chrom, seq in g1.chr_seq.items():
        if chrom not in g2.chr_seq:
            print('%s not in %s skipping' % (chrom, g2))
            continue
            
        print('Indexing chrom %s' % chrom)        

        pidx = compare(seq, g2.chr_seq[chrom], args.window_size, args.seed_size, args.insert_size)

        chrom_idx.update({ chrom : pidx })

    passed = check(g1, g2, chrom_idx)

    if not passed:
        print('At least one failure. Not building index')
        exit(1)

    build_idx(g1, g2, chrom_idx)
