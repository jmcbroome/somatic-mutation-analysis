import pandas as pd
import numpy as np
from Bio import SeqIO as sqio
import argparse
import random

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mutations', help = "Path to a mutation dataframe (output of make_mutation_frame.py)")
    parser.add_argument('-r', '--reference', help = "Path to a reference fasta")
    parser.add_argument('-c', '--codons', help = "Path to a codon dataframe built from your reference and GTF (output of gtf_to_codons.py)")
    parser.add_argument('-g', '--genelist', help = "Text file containing a list of flybase gene identifiers to calculate the ratio on. Default uses all genes in the codons table", default = "")
    parser.add_argument('-p', '--perms', type = int, help = "Number of permutations to perform for building null expectation of missense/synonymous ratio conditioned on mutation type and number distribution.", default = 100)
    args = parser.parse_args()
    return args

args = argparser()
translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

gdf = pd.read_csv(args.codons, sep = '\t').drop("Unnamed: 0", axis = 1)
gdf = gdf.set_index("GID")
mutdf = pd.read_csv(args.mutations, sep = '\t')
genome = sqio.to_dict(sqio.parse(args.reference, format = 'fasta'))

if args.genelist != "":
    genes = [g.strip() for g in open(args.genelist).readlines()]
    gdf = gdf.loc[genes]


chroconv = {"X":"NC_004354.4", '2L':'NT_033779.5', '2R':'NT_033778.4', '3L':'NT_037436.4', '3R':'NT_033777.3'}
#get the prior and latter bases for each mutation for trinucleotide-based analysis
prior = []
latter = []
for i,d in mutdf.iterrows():
    if d.Chro in chroconv:
        prior.append(genome[chroconv[d.Chro]][d.Loc-2])
        latter.append(genome[chroconv[d.Chro]][d.Loc])
    elif d.Chro in genome.keys():
        prior.append(genome[d.Chro][d.Loc-2])
        latter.append(genome[d.Chro][d.Loc])
    else:
        prior.append("N")
        latter.append("N")
mutdf['MutType'] = mutdf.Ref + '>' + mutdf.Alt
mutdf['Prior'] = prior
mutdf['Latter'] = latter
mutdf['Missense'] = mutdf.Type.apply(lambda x:("missense" in x))
mutdf['Synonymous'] = mutdf.Type.apply(lambda x:("synonymous" in x))

#the following corrects for GC bias.
atvc = mutdf[(~mutdf.Synonymous) & (~mutdf.Missense) & (mutdf.Ref.isin(['A','T']))].MutType.value_counts().apply(lambda x:x/2)
gcvc = mutdf[(~mutdf.Synonymous) & (~mutdf.Missense) & (mutdf.Ref.isin(['C','G']))].MutType.value_counts().apply(lambda x:x/2)
occur = pd.concat([atvc, gcvc])
#print("DEBUG: ", occur)
mutdf['TrinType'] = mutdf.Prior.apply(lambda x:x.upper()) + mutdf.Ref + mutdf.Latter.apply(lambda x:x.upper()) + '>' + mutdf.Alt
triatvc = mutdf[(~mutdf.Synonymous) & (~mutdf.Missense) & (mutdf.Ref.isin(['A','T']))].TrinType.value_counts(normalize = True).apply(lambda x:x/2)
trigcvc = mutdf[(~mutdf.Synonymous) & (~mutdf.Missense) & (mutdf.Ref.isin(['C','G']))].TrinType.value_counts(normalize = True).apply(lambda x:x/2)
trioccur = pd.concat([triatvc, trigcvc])

def build_encoding(cdf, genelist = list(gdf.index)):
    #select the subset of the cdf that I'm going to use (default is the whole thing) based on genelist
    sdf = cdf.loc[genelist]
    #get the total count of each codon in sdf
    sumcount = {c:sum(sdf[c]) for c in translate.keys()}
    #initialize
    bases = {b:[] for b in 'ACGT'}
    #for each reference, for each unique codon, record the codon name, the positions that have the reference in question, and the count of this codon in the query set
    for b in 'ACGT':
        for c in translate.keys():
            #count the number of times the codon appears in the query set
            #don't record if its 0.
            count = sumcount[c]
            if count == 0:
                continue
            #check each location to see if it matches the reference, record if so
            for i,cb in enumerate(c):
                if cb == b:
                    bases[b].append([c, i, count])
    return bases
bencode = build_encoding(gdf)
for b in 'ACGT':
    sumc = sum([sb[2] for sb in bencode[b]]) #this sum value is actually the total number of positions assigned to this base among the query set
    for i,e in enumerate(bencode[b]):
        bencode[b][i].append(e[2]/sumc)
        bencode[b][i] = tuple(bencode[b][i])
#the following code generates a distribution of bootstrapped values for the expected ratio of missense to synonymous mutations under random relocation
#then displays the mean and 95% confidence interval for these ratios. These ratios form the denominator for the estimate of DnDs across the genome or for a specific set of genes.
ratios = []
print("Beginning {} ratio permutations".format(args.perms))
for perm in range(args.perms): #just to get an idea of the variance.
    syn = 0
    non = 0
    for mtype in occur.index:
        c = occur[mtype]
        ref = mtype[0]
        #make a number of choices of location with replacement equal to C and collect it as a value count dictionary
        plocs = random.choices(population = bencode[ref], weights = [sb[-1] for sb in bencode[ref]], k = int(c))
        #boil this vector down into a value counts
        un, co = np.unique(plocs, return_counts = True, axis = 0)
        #mtdict = {tuple(u):c for u,c in zip(list(un), list(co))}
        mtdict = {}
        for u, c in zip(list(un), list(co)):
            #np.unique made entries in the first into strings because its dumb.
            mtdict[(u[0], int(u[1]), int(u[2]), float(u[3]))] = c
        #print(mtype, mtdict)
        for info, count in mtdict.items():
            #check whether this mutation type for this specific codon location is synonymous or nonsynonymous
            #and add the count to the appropriate counter
            ncodon = list(info[0])
            assert ncodon[info[1]] == mtype[0]
            ncodon[info[1]] = mtype[-1]
            if translate[''.join(ncodon)] != translate[info[0]]:
                non += count
            else:
                syn += count
    ratios.append(non/syn)

fifth,mean,ninetyfifth = np.quantile(ratios,[0.025,0.5,0.975])
print("Expected baseline ratio 95 confidence interval is {},{} with a mean {} and standard deviation {}. Use this for DnDs/PiNPiS denominators in this analysis.".format(fifth,ninetyfifth,mean,np.std(ratios)))