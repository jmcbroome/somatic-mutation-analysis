import sys
for entry in sys.stdin:
    spent = entry.strip().split()
    nent = [spent[0], spent[1], spent[2]]
    i = 0
    for read in spent[4]:
        vars = []
        if read != '.' and read in 'ACGTN.':
            i += 1
            if spent[5][i-1] not in '.*/':
                if read != spent[2] and read not in 'N.' and int(spent[5][i-1]) >= 2:
                    vars.append(read)
    if len(vars) > 0:
        nent.append(','.join(vars))
    if len(nent) > 3:
        print('\t'.join(nent + [str(i)]))
    
