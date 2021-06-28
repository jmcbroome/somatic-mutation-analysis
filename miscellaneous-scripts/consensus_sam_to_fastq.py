#!/usr/bin/env python

#this tiny script takes a sam input (get from consensus bam with samtools view)
#extracts the constructed consensus sequences and saves the original mapping location as part of an extended read name
#this file can then be used to remap consensi 
#used as a quality control step to try to handle orthologous sequence
#could try stripping Ns but I'm not sure its necessary/helpful to do so?

import sys
for entry in sys.stdin:
    name, _, chro, loc, _, _, _, _, _, seq, qual = entry.strip().split()
    fqname = '_'.join(['@'+name, chro, loc])
    print(fqname)
    print(seq)
    print('+')
    print(qual)