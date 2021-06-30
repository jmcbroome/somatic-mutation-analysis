#!/usr/bin/env python3
import sys
import argparse

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reads', type = str)    
    args = parser.parse_args()
    return args

def read_names(readfile):
    infnames = set()
    with open(readfile) as inf:
        for name in inf:
            infnames.add(name.strip())
    return infnames

def main():
    args = argparser()
    keeps = read_names(args.reads)
    for line in sys.stdin:
        if line[0] == "@":
            print(line.strip())
        else:
            name = line.split()[0]
            if name in keeps:
                print(line.strip())
                
if __name__ == "__main__":
    main()