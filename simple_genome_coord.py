#!/usr/bin/env python3

# 
# Copyright 2019, Christopher Bennett
# 
# This is a prototype script to align two linear sequences (genomes)
# 
# simple-genome-coord is free software: you can redistribute it and/or modify
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
def read_idx(fn):
    idx = {}
    names = []
    with open(fn, 'r') as fi:
        for line in fi:
            line = line.strip().split('\t')
            if line[0].startswith('##'):
                continue
            if line[0].startswith('#ORD'):
                names = line[-2:]
                continue

            try:
                pos = [x.split(':') for x in line[3].split(',')]
                pos.append([str(int(line[1])+1), str(int(line[2])+1)])
            except:
                pos = []

            idx.update({ line[0] : pos })

    if not names:
        print("No header detected")
        exit(1)

    return(idx, names)

def get_indicies(pidx, i):
    master = [[], []]
    for idx in pidx:
        p1 = int(idx[0])
        p2 = int(idx[1])
        if i == 2:
            master[0].append(p2)
            master[1].append(p1 - p2)
        elif i == 1:
            master[0].append(p1)
            master[1].append(p2 - p1)
    master[1].pop()
    master[1].insert(0, 0)

    return(master)

def convert_coord(nidx, idx_base, fn, form):
    idx, names = read_idx(nidx)

    if '.fa' in idx_base:
        idx_base = idx_base.replace('.fa', '')
    if idx_base not in names:
        print('Index base name not found in index file')
        exit(1)
    if not os.path.exists(fn):
        print('File not found')
        exit(1)

    base = names.index(idx_base) + 1

    for chrom, pos in idx.items():
        if not idx[chrom]:
            continue
        par_pos = get_indicies(pos, base)
        idx[chrom] = par_pos

    fi = open(fn, 'r')
    fo = open('conv_%s' % fn, 'w')

    pos, prev_chr = 0, ''
    for line in fi:
        if line.startswith('#'):
            fo.write(line)
            continue
    
        line = line.split('\t')
        chrom = line[0]
        if 'MT' in chrom: # temp to convert MT to M
            line[0] = 'M'

        if not (chrom in 'MTYX' or chrom.isdigit()):
            continue

        if 'chr' in chrom:
            chrom = chrom.replace('chr', '')
        # else:
        #    line[0] = 'chr' + chrom if chrom.isdigit() else chrom
        if chrom != prev_chr:
            pos = 0

        if not idx[chrom]:
            fo.write('\t'.join(line))
            continue

        if int(line[1]) < idx[chrom][0][pos]:
            if pos == 0: 
                fo.write('\t'.join(line))
                continue
            line[1] = str(int(line[1]) + idx[chrom][1][pos])
        else:
            pos += 1
            line[1] = str(int(line[1]) + idx[chrom][1][pos])

        prev_chr = chrom
        fo.write('\t'.join(line))

    fi.close()
    fo.close()


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Simple Genome Coordinant Converter')
    parser.add_argument('--idx',
                        dest='nidx',
                        type=str,
                        required=True,
                        help='Genome Coordinant index')
    parser.add_argument('--idx-base',
                        dest='idx_base',
                        type=str,
                        help='Name of genome in index used in infile')
    parser.add_argument('-i', '--infile',
                        dest='fn',
                        type=str,
                        required=True,
                        help='File with sorted coordiants to be converted')
    parser.add_argument('--type',
                        dest='form',
                        type=str,
                        help='File format to convert (only VCF supported right now)')

    args = parser.parse_args()

    convert_coord(args.nidx, args.idx_base, args.fn, args.form)
