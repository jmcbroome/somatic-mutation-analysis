#small script which counts the total depth across each reference base in a pileup.
import sys
count = {b:0 for b in 'ACGT'}
for entry in sys.stdin:
    chro, loc, ref, dep, alt, qual = entry.strip().split()
    if ref in count:
        count[ref] += int(dep)
for k,v in count.items():
    print(k + "\t" + str(v))
    