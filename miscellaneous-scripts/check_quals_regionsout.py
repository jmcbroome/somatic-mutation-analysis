#!/usr/bin/env python3
import sys
import numpy as np
def phred_code_qual(sym):
    if isinstance(sym, str):
        return ord(sym) - 33
    elif isinstance(sym, int):
        return chr(sym+33)
    else:
        print("Error: quality score conversion not understood. Check input")
        raise ValueError

quald = {r:{a:[] for a in 'ACGT' if a != r} for r in 'ACGT'} 

for entry in sys.stdin:
    try:
        ref, cigar, con, conq = entry.strip().split('\t')
    except:
        # print('Cant parse entry', entry)
        continue
    if cigar.count('I') > 0 or cigar.count('D') > 0:
        continue
    if len(ref) != len(con):
        #print('Aligned reference and consensus dont match lengths', cigar, len(ref), len(con), len(conq), ref, con, conq)
        continue
    if len(con) != len(conq):
        #print("Target sequence and quality sequence dont match lengths", cigar, len(ref), len(con), len(conq), ref, con, conq)
        continue
    for i,b in enumerate(con):
        if ref[i].upper() != b.upper() and ref[i].upper() != 'N' and b.upper() != "N" and ref[i].upper() != ' ' and b.upper() != ' ':
            qs = phred_code_qual(conq[i])
            quald[b.upper()][ref[i].upper()].append(qs)

for ref, alts in quald.items():
    for a, qs in alts.items():
        binc, bine = np.histogram(qs, bins = 10)
        print(ref, a, np.mean(qs), binc, bine)