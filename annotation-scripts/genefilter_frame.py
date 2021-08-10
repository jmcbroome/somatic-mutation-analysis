import pandas as pd
import numpy as np
import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mutations', help = "Path to a mutation dataframe (output of make_mutation_frame.py)")
    parser.add_argument('-c', '--codons', help = "Path to a codon dataframe built from your reference and GTF (output of gtf_to_codons.py)")
    parser.add_argument('-g', '--bygene', help = "Name of a one-gene-per-row dataframe of mutation information. Default is genedf.tsv", default = 'genedf.tsv')
    parser.add_argument('-o', '--coding', help = "Name of a filtered and reformatted version of the mutation dataframe, used for further analysis of coding mutations downstream. Default coding_mutdf.tsv", default = 'coding_mutdf.tsv')
    parser.add_argument('-t', '--rate', type = float, help = "Maximum number of coding mutations per base allowed for a gene to be included. Default is .2", default = .2)
    parser.add_argument('-p', '--countcap', type = int, help = "Maximum number of coding mutations allowed for a single gene period. Default is 1000", default = 1000)
    parser.add_argument('-i', '--ratio', type = float, help = "Maximum raw ratio of missense to synonymous mutations allowed for an individual gene before removal. Default 20x more missense than synonymous", default = 20)
    parser.add_argument('-d', '--depths', help = "Path to an output table containing depth and normalized count values for spectra generation.", default = "restricted_depth_values.tsv")
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
mutdf = pd.read_csv(args.mutations, sep = '\t')
print("Saving depth table.")
dvc = mutdf.Depth.value_counts(normalize=True)
with open(args.depths, "w+") as outf:
    for i in dvc.index:
        print(str(i) + "\t" + str(dvc[i]), file = outf)
cdf = pd.read_csv(args.codons, sep = '\t')
cdf = cdf.set_index('GID')
print("Constructing by-gene mutation table.")
gtd = {k:{'missense':0, 'synonymous':0, 'other':0, 'mc':0, 'sc':0, 'oc':0} for k in cdf.index}
for i,d in mutdf[mutdf.SampleFreq < .25].iterrows():
    #split the d.GID column into values
    try:
        gvs = d.GID.split("|")
    except AttributeError:
        continue #some gids are float for some unknowable reason. blame the gtf.
    #get the same info for effects
    effects = d.Type.split("|")
    #check which are legit gene IDs I have codon info for
    for i,g in enumerate(gvs):
        if g in gtd:
            #save its info there.
            ef = effects[i]
            if ef == 'missense_variant':
                gtd[g]['missense'] += d.Pi
                gtd[g]['mc'] += 1
            if ef == 'synonymous_variant':
                gtd[g]['synonymous'] += d.Pi
                gtd[g]['sc'] += 1
            else:
                gtd[g]['other'] += d.Pi
                gtd[g]['oc'] += 1
genedf = {k:[] for k in ['GID', 'Strand', 'Ratio', 'Length', 'PiNPiS', 'mis_pi', 'mis_c', 'syn_pi', 'syn_c', 'other_pi', 'other_c']}
for g, md in gtd.items():
    #calculate pin/pis
    if md['sc'] > 0:
        genedf['PiNPiS'].append((md['missense'] / md['synonymous']) / cdf.loc[g]['Ratio'])
        genedf['GID'].append(g)
        genedf['Strand'].append(cdf.loc[g]['Strand'])
        #get the ratio for this gene.
        genedf['Ratio'].append(cdf.loc[g]['Ratio'])
        #calculate the length too
        genedf['Length'].append(sum([cdf.loc[g][c] for c in translate.keys()]))
        #save the total pi
        genedf['mis_pi'].append(md['missense'])
        genedf['mis_c'].append(md['mc'])
        genedf['syn_pi'].append(md['synonymous'])
        genedf['syn_c'].append(md['sc'])
        genedf['other_pi'].append(md['other'])
        genedf['other_c'].append(md['oc'])
genedf = pd.DataFrame(genedf)
genedf['Rate'] = [(d.mis_c + d.syn_c)/d.Length for i,d in genedf.iterrows()]
genedf = genedf[genedf.Rate < args.rate] #filter out the high-mutation ones as being possible orthologous mismapping.
#want to test these filters more.
genedf = genedf[(genedf.mis_c + genedf.syn_c < args.countcap)] #stripping outliers.
genedf = genedf[genedf.PiNPiS < args.ratio] #more outliers out
genedf.to_csv(args.bygene,sep='\t')

print("Generating filtered coding table.")
def generate_smutdf(mutdf):
    #because the way snpeff records annotation is nightmarish, 
    #I have to use a bunch of code to parse out whether any given mutation is synonymous, missense, or whatever
    #going to make a fresh frame that has a unique entry for each gene involved with the mutation
    smutdf = {k:[] for k in ['Chro','Loc','Effect','GID','Pi', 'Overlap', "SSN", "SampleFreq", "Depth"]}
    for i,d in mutdf[mutdf.SampleFreq < .25].iterrows():
        #split the d.GID column into values
        try:
            gvs = d.GID.split("|")
        except AttributeError:
            effect = 'Floater'
            continue #some gids are float for some unknowable reason. blame the gtf.
        #get the same info for effects
        effects = d.Type.split("|")
        #check which are legit gene IDs I have codon info for
        for i,g in enumerate(gvs):
            #only include genes which passed the codon writer filters (e.g. having valid frame, start-stop codons).
            if g in cdf.index:
                #save its info there.
                ef = effects[i]
                if ef == 'missense_variant':
                    smutdf['Chro'].append(d.Chro)
                    smutdf['Loc'].append(d.Loc)
                    smutdf['Effect'].append('missense')
                    smutdf['GID'].append(g)
                    smutdf['Pi'].append(d.Pi)
                    smutdf['Overlap'].append(len(gvs))
                    smutdf["SSN"].append(d.SSN)
                    smutdf["SampleFreq"].append(d.SampleFreq)
                    smutdf["Depth"].append(d.Depth)
                if ef == 'synonymous_variant':
                    smutdf['Chro'].append(d.Chro)
                    smutdf['Loc'].append(d.Loc)
                    smutdf['Effect'].append('synonymous')
                    smutdf['GID'].append(g)
                    smutdf['Pi'].append(d.Pi)
                    smutdf['Overlap'].append(len(gvs))
                    smutdf["SSN"].append(d.SSN)
                    smutdf["SampleFreq"].append(d.SampleFreq)
                    smutdf["Depth"].append(d.Depth)
                else:
                    smutdf['Chro'].append(d.Chro)
                    smutdf['Loc'].append(d.Loc)
                    smutdf['Effect'].append('other')
                    smutdf['GID'].append(g)
                    smutdf['Pi'].append(d.Pi)
                    smutdf['Overlap'].append(len(gvs))
                    smutdf["SSN"].append(d.SSN)
                    smutdf["SampleFreq"].append(d.SampleFreq)
                    smutdf["Depth"].append(d.Depth)
    smutdf = pd.DataFrame(smutdf)
    #also strip duplicate entries.
    smutdf = smutdf.drop_duplicates(subset = ['Loc', 'GID', 'SSN', 'Pi', 'Effect'])
    return smutdf

smutdf = generate_smutdf(mutdf)
#filter smutdf based on the content of genedf.
#the resulting dataframe will contain only coding mutations for genes which pass the above quality checks.
smutdf = smutdf[(smutdf.Effect != 'other') & (smutdf.GID.isin(genedf.GID))]
smutdf.to_csv(args.coding,sep='\t')
