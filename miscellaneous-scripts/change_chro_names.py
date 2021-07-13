#quick script to change the names on the dmel 6.34 from genbank to X, 2L, 3R, etc etc.
import sys
for entry in sys.stdin:
    if entry[0] == ">":
        print(">" + entry.strip().split()[-1])
    else:
        print(entry.strip())