import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import math

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mutations', help = "Path to a mutation dataframe (output of make_mutation_frame.py)")
    parser.add_argument('-b', '--basecounts', help = "Path to a text file of base counts to use for rate calculation.")
    parser.add_argument('-o', '--output', help ="Choose a prefix for generated figures (Figure 1 of the manuscript).")
    args = parser.parse_args()
    return args
args = argparser()
mutdf = pd.read_csv(args.mutations, sep = '\t')
mutdf['Type'] = mutdf.Ref + ">" + mutdf.Alt
#print some basic statistics.
ssnvc = mutdf.SSN.value_counts()
for ssn in ssnvc.index:
    print("Individual {} has {} total mutations.".format(ssn, ssnvc[ssn]))
    mtvc = mutdf[mutdf.SSN == ssn].Type.value_counts()
    for t in mtvc.index:
        print("with {} {} mutations".format(mtvc[t],t))

print("Generating graph of mutation rates.")
def read_basecounts(bf):
    bcd = {}
    with open(bf) as inf:
        for entry in inf:
            spent = entry.strip().split()
            if len(spent) == 3:
                ind,base,count = spent
                if ind not in bcd:
                    bcd[ind] = {b:0 for b in "ACGT"}
                bcd[ind][base] = int(count)
            elif len(spent) == 2:
                base,count = spent
                bcd[base] = int(count)
    return bcd

bcd = read_basecounts(args.basecounts)
pairs = []
for a in 'ACGT':
    for b in 'ACGT':
        if a != b:
            pairs.append((a,b))

cdf = {k:[] for k in ['ind', 'Mutation Type', 'Rate (Detections per Site per Depth)']}
for ssn in mutdf.SSN.value_counts().index:
    for ref, alt in pairs:
        subdf = mutdf[(mutdf.SampleFreq < .25) & (mutdf.SSN == ssn) & (mutdf.Ref == ref) & (mutdf.Alt == alt)]
        count = sum([math.floor(v.Depth * v.SampleFreq) for i,v in subdf.iterrows()])
        cdf['ind'].append(ssn)
        cdf['Mutation Type'].append(ref + '>' + alt)
        if ssn in bcd.keys():
            c = bcd[ssn][ref]
        else:
            c = bcd[ref]
        cdf['Rate (Detections per Site per Depth)'].append(count/c)
cdf = pd.DataFrame(cdf)

ax = sns.boxplot(x = 'Mutation Type', y = 'Rate (Detections per Site per Depth)', data = cdf, color = 'grey') #hue = 'ind' if I want to tag it in
#ax.set_yticklabels([0,0, '2e-7', '4e-7', '6e-7', '8e-7', '1e-6', '1.2e-6'])
#print(list(ax.get_yticks()))
ax.set_yticklabels(['{:.2e}'.format(v) for v in list(ax.get_yticks())])
ax.set_xlabel("Mutation Type")
ax.set_ylabel("Rate (Detections per Site per Depth)")
ax.set_title('Somatic Mutation Rates')
plt.savefig(args.output + '_mutations_rates.png',dpi=800)

y,x=np.histogram(np.log10(mutdf.SampleFreq),bins=10,density=True)
x = [round(10**(i),3) for i in x]
ax=sns.barplot(x=x[1:],y=y,color='grey')
#ax.set_xticklabels(['<'+str(xv) for xv in x])
ax.set_title("Binned Log Circle Frequency Spectrum")
ax.set_ylabel("Density")
ax.set_xlabel("Maximum Sample Frequency")
plt.savefig(args.output + "_cfs.png",dpi=800)