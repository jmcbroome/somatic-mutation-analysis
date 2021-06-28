#!/usr/bin/env python

#this script takes the consensus remap alignments and strips out consensus sequences that remapped to places other than the original mapping location, indicating that they're likely some kind of orthology problem.
#feed this in a bam from a bwa mem with the -a parameter activated to pick up secondary alignments
#then this will retain only primary alignments which have secondary alignments which are all either within a short range or have a short cigar on them.

import sys
import re

def check_maxmatch(cig, maxm = 20):
    try:
        match = int(re.split('[A-Z]', cig.partition('M')[0])[-1])
    except:
        return False
    if match < maxm: #I will allow ones where the match is less than the max value to always count as close enough.
        return True
    else:
        return False

def verify_entry(group):
    #extract the original location information from the read name, make sure the primary alignment matches up with it. and that there's only one primary alignment.
    #primary alignment has flag 0, secondaries have 250+, supplementary have 2000+
    #after verifying the primary is where it should be, check each secondary/supplementary that it's either 1. short (<25 bp matches) or 2. within 300 bp of the primary
    #print the primary only if and only if the above is true for all secondary/supplementary alignments
    #name information first, just grabbed from the first entry in group on the assumption they'll all be the same (should be).
    ni = group[0][0].split("_")
    oname = ni[0]
    oloc = ni[-1]
    ochro = '_'.join(ni[1:-1])
    if [d[1] for d in group].count('0') == 1:
        #first, identify and separate the primary
        primary = [d for d in group if d[1] == '0'][0]
        secondaries = [d for d in group if d[1] != '0']
        #now verify that the primary mapping location aligns with the original initial mapping, extracted from the name
        if primary[2] == ochro and  abs((int(oloc) - int(primary[3]))) < 300:
            #now verify that all of the secondaries either have a short cigar or are near the primary mapping location
            if all([check_maxmatch(s[5]) or abs(int(s[3]) - int(primary[3])) < 300 for s in secondaries]):
                #if all of the above is true, I can save the primary alignment!
                print('\t'.join(primary))
    #        else:
    #            print('secondaries do not support primary location')
    #    else:
    #        print('new primary does not match original')
    #else:
    #    print(oname, [d[1] for d in group])
    #    print("wrong primary count")

curname = None
group = []
for sament in sys.stdin:
    #print all headers.
    if sament[0] == '@':
        print(sament.strip())
        continue
    spent = sament.strip().split()
    if spent[0] != curname:
        #print the old group, or don't.
        if group != []: #if there's a group to dump to stdout, verify its integrity, dump it if it passes checks
            verify_entry(group)
        group = [spent]
        curname = spent[0]
    else:
        group.append(spent)