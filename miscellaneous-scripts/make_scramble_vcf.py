#this script takes a vcf, strips snpeff annotations, randomizes locations based on hardcoded chromosome lengths, 
import argparse
from Bio import SeqIO as sqio
import random

#global variable dictionary for resolving chromosome name issues.
chro_conv = {
    'NC_004354.4':'X',
    'NT_033779.5':'2L',
    'NT_033778.4':'2R',
    'NT_037436.4':'3L',
    'NT_033777.3':'3R'
}
chro_conv.update({v:k for k,v in chro_conv.items()})

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genome', help = 'path to the genome fasta')
    parser.add_argument('-f', '--vcf', help = 'path to an input vcf to scramble')
    parser.add_argument('-o', '--output', help = 'name of the outputs prefix, default "scramble_"', default = 'scramble_')
    args = parser.parse_args()
    return args

def randomize_location(chro, loc, alts, genome):
    nloc = random.randrange(2, len(genome[chro])-1)
    nref = genome[chro][nloc-1]
    nalts_dict = {a:random.choice([b for b in 'ACGT' if b != nref]) for a in alts}
    return nloc, nref, nalts_dict

def main():
    args = argparser()
    genome = sqio.to_dict(sqio.parse(args.genome, format = 'fasta'))  
    with open(args.vcf) as inf:
        with open(args.output + args.vcf, 'w+') as outf:
            for entry in inf:
                if entry[0] != '#':
                    chro, loc, x, ref, alts, y, pv, info = entry.strip().split()
                    if chro not in chro_conv:
                        continue
                    chro = chro_conv[chro]
                    #remove prior annotation information
                    if ';ANN=' in info:
                        info = info.split(";ANN=")[0]
                    nloc, nref, nalts_dict = randomize_location(chro, loc, alts, genome)
                    nalts = ','.join([nalts_dict[a] for a in alts.split(',')])
                    #update the info line by replacing new values in the third section onwards.
                    ninfo_groups = []
                    for v in info.split(';')[2:]:
                        if v[0] in nalts_dict:
                            ninfo_groups.append(nalts_dict[v[0]] + v[1:])
                    ninfo = ';'.join(info.split(';')[:2]+ninfo_groups)
                    print('\t'.join([str(v) for v in [chro_conv[chro], nloc, x, nref, nalts, y, pv, ninfo]]), file = outf)

if __name__ == '__main__':
    main()
