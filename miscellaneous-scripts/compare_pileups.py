#!/usr/bin/env python3

#import
import argparse
import numpy as np

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-a', '--first')
    parser.add_argument('-b', '--second')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    return args

def read_pileup(path):
    chroms = {k:[] for k in ['W','2L','2R','3L','3R','X']}
    with open(path) as inf:
        for entry in inf:
            spent = entry.strip().split()
            if spent[0] == 'CP001391.1':
                chro = "W" #wolbachia.
            elif spent[0] == 'NT_479533.1':
                chro = '2L'
            elif spent[0] == 'NT_479535.1':
                chro = '3L'
            elif spent[0] == 'NT_479534.1':
                chro = '2R'
            elif spent[0] == 'NT_479536.1':
                chro = '3R'
            elif spent[0] == 'NC_029795.1':
                chro = 'X'
            else:
                continue #skip node and small scaffold 
            chroms[chro].append((int(spent[1]),int(spent[3])))
    return chroms

def combine_pileup(b1, b2):
    bdeps = {}
    for chro, locs in b1.items():
        if chro not in bdeps:
            bdeps[chro] = {}

        for loc, dep in locs:
            bdeps[chro][loc] = [dep, 0]

    for chro, locs in b2.items():
        for loc, dep in locs:
            if loc in bdeps[chro]:
                bdeps[chro][loc][1] = dep
            else:
                bdeps[chro][loc] = [0, dep]
    return bdeps

def main():
    args = argparser()
    #insert code
    b1 = read_pileup(args.first)
    b2 = read_pileup(args.second)
    print("Sum of depth in pileup b1:", sum([sum([v[1] for v in k]) for k in b1.values()]))
    print("Sum of depth in pileup b2:", sum([sum([v[1] for v in k]) for k in b2.values()]))
    bdeps = combine_pileup(b1,b2)
    with open(args.output, 'w+') as outf:
        for chro, locd in bdeps.items():
            for l, pd in sorted(list(locd.items())):
                print('\t'.join([str(v) for v in [chro, l, pd[0], pd[1], pd[1]-pd[0]]]), file = outf)

if __name__ == "__main__":
    main()