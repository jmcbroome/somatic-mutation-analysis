import sys
seen = {}
total = 0
for entry in sys.stdin:
    if entry[0] != '#':
        spent = entry.strip().split()
        anns = spent[-3]
        total += int(anns.split(';')[0].strip('DP='))    
        try:
            div = anns.split('|')
            if len(div) > 1:
                atype = div[1]
                seen[atype] = seen.get(atype,0) + 1
        except:
            continue
for k,v in seen.items():
    print(k, v / total)
