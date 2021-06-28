#!/usr/bin/env python3

#import
import argparse
import numpy as np
import skbio.alignment as skaln
#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-a', '--first')
    parser.add_argument('-b', '--second')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    return args

def subtract_quals(q1,q2):
    qd = []
    for i,v in enumerate(q1):
        d = int(q2[i]) - int(v)
        qd.append(d)
    return ','.join([str(v) for v in qd])

def encode_match(s1,s2,order = 0):
    encoding = {('A','T'):'Q',
                ('A','G'):'W',
                ('A','C'):'E',
                ('T','A'):'R',
                ('T','G'):'Y',
                ('T','C'):'U',
                ('C','A'):'I',
                ('C','T'):'O',
                ('C','G'):'P',
                ('G','A'):'S',
                ('G','T'):'D',
                ('G','C'):'F',
                ('N','A'):'H',
                ('N','T'):'J',
                ('N','C'):'K',
                ('N','G'):'L',
                ('A','N'):'Z',
                ('T','N'):'X',
                ('C','N'):'V',
                ('G','N'):'B'}
    #assume s1 and s2 are equal length
    fs = []
    for i,b in enumerate(s1):
        if b == s2[i]:
            fs.append(b)
        elif order == 0: #seq 1 is s1
            fs.append(encoding[(b,s2[i])])
        else:
            fs.append(encoding[(s2[i],b)])
    return ''.join(fs)

def get_alignment(query,target):
    #query = shorter sequence
    #target = longer sequence
    seqaln = skaln.StripedSmithWaterman(query, gap_open_penalty = 5, gap_extend_penalty = 2, match_score = 1, mismatch_score = -1, zero_index = False)
    aln = seqaln(target)
    return aln.query_begin, aln.query_end, aln.target_begin, aln.target_end_optimal

def get_comp(first, second):
    #inputs are split versions of a bam entry.
    #output is a string of read name, map location from the first file (e.g put russ's first), consensus between them both with encoding for differences, and a vector of the differences in their count per base.
    #extract the sequence and qual lines (last 2 entries)
    name = first[0]
    chro = first[2]
    loc = first[3]
    fseq = first[-2]
    fq = first[-1]
    sseq = second[-2]
    sq = second[-1]
    sq1,sq2 = sorted([(fseq,fq,0),(sseq,sq,1)], key = lambda x:len(x[0]))
    qb, qe, tb, te = get_alignment(sq1[0],sq2[0])
    if sq1[2] == 0: #first was the query/shorter
        ncon = encode_match(sq1[0][qb:qe],sq2[0][tb:te])
        nq = subtract_quals(sq1[1][qb:qe],sq2[1][tb:te])
    elif sq1[2] == 1: #second was query
        ncon = encode_match(sq2[0][qb:qe],sq1[0][tb:te])
        nq = subtract_quals(sq2[1][qb:qe],sq1[1][tb:te])
    datastr = '\t'.join([name, chro, loc, ncon, nq])
    return datastr

def main():
    args = argparser()
    with open(args.output, 'w+') as outf:
        with open(args.first) as inf1:
            with open(args.second) as inf2:
                #assume these sams are sorted by read name and contain the same set of read names.
                #samtools sort -n <bam> | samtools view > sorted.sam
                #first iterate both file objects past the headers.
                while True:
                    firent = inf1.readline()
                    if firent[0] != '#':
                        firent = firent.strip().split()
                        break
                while True:
                    secent = inf2.readline()
                    if secent[0] != '#':
                        secent = secent.strip().split()
                        break
                #return to primary iteration, again using a while loop.
                while True:
                    if firent == '' or secent == '': #hit the end of the file.
                        break
                    #the variables saved from the last iteration are my first pair
                    assert firent[0] == secent[0] #names should match
                    matchstr = get_comp(firent, secent)
                    print(matchstr, file = outf)
                    #read the next entry.
                    firent = inf1.readline().strip().split()
                    secent = inf2.readline().strip().split()

if __name__ == "__main__":
    main()