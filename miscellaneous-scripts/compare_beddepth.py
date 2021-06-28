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
    args = parser.parse_args()
    return args

def read_bedgraph(path):
    b = []
    with open(path) as inf:
        for entry in inf:
            spent = entry.strip().split()
            b.append(spent)
    return b

def main():
    args = argparser()
    #insert code
    b1 = read_bedgraph(args.first)
    b2 = read_bedgraph(args.second)
    print("Sum of depth in bedgraph b1:", sum([int(b[1])*int(b[2]) for b in b1]))
    print("Sum of depth in bedgraph b2:", sum([int(b[1])*int(b[2]) for b in b2]))    
    for i, b1d in enumerate(b1):
        b2d = b2[i]
        diff = int(b2d[2]) - int(b1d[2]) 
        print(i, diff, i*diff)

if __name__ == "__main__":
    main()