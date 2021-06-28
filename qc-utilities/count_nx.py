#go through a pileup passed to stdin and count the number of sites where an alternative appears once, or twice, and so forth.
#used for qc'ing edge filtering
#note that a pileup may include pcr duplicates, so focus on the 1x and 2x sightings for comparison.
import sys
counts = {}
for entry in sys.stdin:
    chro, loc, ref, dep, alt, qual = entry.strip().split()
    subcount = {b:0 for b in 'ACGT'}
    if int(dep) < 10:
        continue #don't want low-depth germline mutations to count.
    for a in alt:
        if a in subcount:
            subcount[a] += 1
    for b, c in subcount.items():
        if c > 0:
            if c not in counts:
                counts[c] = 0
            counts[c] += 1
print("XSeen", "SiteCount")
for c in sorted(list(counts.keys())):
    print(c, counts[c])