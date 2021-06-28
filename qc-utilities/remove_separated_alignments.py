#!/usr/bin/env python

#this script takes a SAM format file from standard input, collects read groups, and prints read groups where all secondary and supplementary alignments are within 300 bases of one another, while dropping reads otherwise.
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

curname = None
group = []
for sament in sys.stdin:
    #print all headers.
    if sament[0] == '@':
        print(sament.strip())
    spent = sament.strip().split()
    if spent[0] != curname:
        #print the old group, or don't.
        if group != []: #if there's a group to dump to stdout, dump it
            #check to see if all of them have the same chromosome.
            #entry number 2 (name, flag, chro, loc)
            if len(set([g[2] for g in group])) == 1: #exactly 1 unique chromosome string
                #calculate pairwise differences between locations for each value, assert they're all less than 300.
                if all([all([int(g[3]) - int(og[3]) < 300 or check_maxmatch(og[5]) for og in group]) or check_maxmatch(g[5]) for g in group]):
                #if True: #disable this part of the filter.
                    for g in group:
                        print('\t'.join(g))
        #start a new group regardless of the outcome above.
        group = []
        curname = spent[0]
    else:
        group.append(spent)