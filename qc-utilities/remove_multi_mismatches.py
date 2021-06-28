#!/usr/bin/env python3

#this script takes a consensus sam that's been passed through "samtools calmd"


import argparse
import sys
import logging

logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action = 'store_true', help = "Print status updates.")
    parser.add_argument('-m', '--max_mismatches', type = int, default = 1, help = 'Set to a maximum number of mismatches allowed within a single read. Default 1')
    args = parser.parse_args()
    return args

def verify_entry(sament):
    #split up the entry, assert that there are "=" in the line (from samtools calmd)
    sam = sament.strip().split()
    #sequence is entry 9
    seq = sam[9]
    if seq.count("=") == 0:
        if seq.count("N") == len(seq): #some lines are just all N. don't want those either but they're not an error with the setup really.
            return 'J' #j for junk!
        else:
            log.error("Entry does not contain equals signs- did you use samtools calmd before passing the sam in? N count {}".format(seq.count("N")))
            return 0
    #count the number of ACGT that occur in it (don't care about Ns)
    mism = sum([seq.count(b) for b in 'ACGT'])
    if mism < 2: 
        #if it doesn't have multiple mismatches, print it back to standard out
        print(sament.strip())
    return mism

def main():
    args = argparser()
    mctrack = {'J':0}
    for entry in sys.stdin:
        if entry[0] == '@':
            #give header lines back without processing
            print(entry.strip())
            continue
        mc = verify_entry(entry)
        mctrack[mc] = mctrack.get(mc, 0) + 1

    if args.verbose:
        log.info("{} reads with single mismatches, {} reads with grouped mismatches, {} junk reads, {} total reads".format(mctrack[1], sum([v for k,v in {k:v for k,v in mctrack.items() if k != 'J'}.items() if k > 1]), mctrack['J'], sum(mctrack.values())))

if __name__ == "__main__":
    main()