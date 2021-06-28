import sys
seq = ''
title = None
for entry in sys.stdin:
    if entry[0] == '>':
        if seq != '':
            if seq[:3] == 'ATG' and len(seq)%3 == 0: #if it starts with a start and is the right length
                end = seq[-3:]
                if end in ['TAG','TAA','TGA']:
                    print(title)
                    print(seq)
                else:
                    seq = ''
                    continue
            else: 
                seq = ''
                continue
        else:
            seq = ''
            continue
        #reset
        seq = ''
        title = entry.strip()
    else:
        seq += entry.strip()