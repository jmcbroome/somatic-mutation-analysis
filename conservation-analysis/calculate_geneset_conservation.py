#!/usr/bin/env python3

#this large script replicates code largely hashed out in and copied from a series of analysis notebooks
#it takes a mutation dataframe and some additional information about genes and produces one or two dataframes of output representing the overall conservation of those genes.

import argparse
import numpy as np
import pandas as pd
import random
import seaborn as sns
import matplotlib.pyplot as plt
import math
from tqdm import tqdm
from scipy.stats import binom, mannwhitneyu, percentileofscore
import gffutils
import sys

translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def construct_genedf(mutdf, cdf):
    '''
    Code to construct a by-gene mutation frame for bootstrapping.
    '''
    gtd = {k:{'missense':0, 'synonymous':0, 'other':0, 'mc':0, 'sc':0, 'oc':0} for k in cdf.index}
    print("QC: Mutdf shape", mutdf[mutdf.SampleFreq < .25].shape, file = sys.stderr)
    c = 0
    nc = 0
    for i,d in mutdf[mutdf.SampleFreq < .25].iterrows():
        #split the d.GID column into values

        try:
            gvs = d.GID.split("|")
            #c += 1
        except AttributeError:
            #nc += 1
            continue #some gids are float for some unknowable reason. blame the gtf.
        #get the same info for effects
        effects = d.Type.split("|")
        #check which are legit gene IDs I have codon info for
        for i,g in enumerate(gvs):
            if g in gtd:
                c += 1
                #save its info there.
                ef = effects[i]
                if ef == 'missense_variant':
                    gtd[g]['missense'] += d.Pi
                    gtd[g]['mc'] += 1
                elif ef == 'synonymous_variant':
                    gtd[g]['synonymous'] += d.Pi
                    gtd[g]['sc'] += 1
                else:
                    gtd[g]['other'] += d.Pi
                    gtd[g]['oc'] += 1
            else:
                nc += 1
    print("QC: gtd", list(gtd.keys())[:10], file = sys.stderr)
    print("QC: successful gid filtering", c, nc, file = sys.stderr)
    #once I have a dictionary of information, I can calculate a frame
    genedf = {k:[] for k in ['GID', 'Strand', 'Ratio', 'Length', 'PiNPiS', 'mis_pi', 'mis_c', 'syn_pi', 'syn_c', 'other_pi', 'other_c', 'Rate']}
    c = 0
    nc = 0
    for g, md in gtd.items():
        #calculate pin/pis
        if md['sc'] > 0:
            c += 1
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
            genedf['Rate'].append((md['mc'] + md['sc'] + md['oc'])/sum([cdf.loc[g][c] for c in translate.keys()]))
        else:
            nc += 1
    print("QC: filter passed entries", c, nc, file = sys.stderr)
    genedf = pd.DataFrame(genedf)
    print('QC: genedf shape pre filtering', genedf.shape)
    genedf = genedf[genedf.Rate < .1] #filter out the high-mutation ones as being possible orthologous mismapping.
    genedf = genedf[(genedf.mis_c + genedf.syn_c < 1000)] #stripping outliers.
    genedf = genedf[genedf.PiNPiS < 20] #more outliers out    
    return genedf

def parse_genes(txt):
    '''
    Read in a simple gene text file.
    '''
    gn = []
    with open(txt) as inf:
        for entry in inf:
            gn.append(entry.strip())
    return gn

def point_bootstrap(group, genedf, ratio = 2.75, size = 1000):
    '''
    This is a block of code to create a bootstrap distribution of ratio estimates for a group of genes by sampling genes with replacement to compose the final set
    This is intended to identify potential problems with individual gene variance, as genes vary widely in size and number of mutations
    And to get accurate depictions of the conservation of an overall group, not have the signal be driven by a single large gene that may be conserved for other reasons not related to this term
    '''
    sdf = genedf[genedf.GID.isin(group)]
    perms = []
    for i in range(size):
        sample = sdf.sample(frac=1, replace = True)
        ratio = (sum(sample.mis_pi)/sum(sample.syn_pi))/ratio
        perms.append(ratio)
    return perms

def split_go_term(mutdf, gdb):
    #split the mutdf into subframes stored in a dictionary.
    god = {}
    for gid in mutdf.GID.value_counts().index:
        try:
            terms = gdb[gid].attributes['Ontology_term']
        except:
            continue
        for t in terms:
            if t not in god:
                god[t] = set()
            god[t].add(gid)
    god = {k:v for k,v in god.items() if 20 < len(v) < 1000} #hard-coded limits to remove outlier terms or terms without enough information.
    tdfs = {}
    for term, gv in god.items():
        subdf = mutdf[mutdf.GID.isin(gv)]
        tdfs[term] = subdf
    return tdfs

def build_encoding(cdf, genelist):
    #select the subset of the cdf that I'm going to use (default is the whole thing) based on genelist
    sdf = cdf.loc[genelist].dropna() #this introduces a ton of blank entries for some reason.
    #print(sdf)
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

#here's a good method for calculating ratio as the joint likelihood of a nonsynonymous divided by the joint likelihood of a synonymous mutation
def get_ratio(occur, cdf, genelist):
    #get the bencode object for this set of genes
    #this script is a wrapper for the above + a likelihood calculation
    bencode = build_encoding(cdf, genelist)
    for b in 'ACGT':
        sumc = sum([sb[2] for sb in bencode[b]]) #this sum value is actually the total number of positions assigned to this base among the query set
        for i,e in enumerate(bencode[b]):
            bencode[b][i].append(e[2]/sumc)
            bencode[b][i] = tuple(bencode[b][i])
    non_like = 0
    syn_like = 0
    for ref, codonvec in bencode.items():
        for cdata in codonvec:
            for alt in 'ACGT':
                if ref != alt:
                    #check whether this is a synonymous change or not.
                    ncodon = list(cdata[0])
                    ncodon[cdata[1]] = alt
                    like = occur[ref + ">" + alt] * cdata[-1]
                    if translate[''.join(ncodon)] != translate[cdata[0]]:
                        #this is missense. its likelihood will go there
                        #now calculate that likelihood- the chance of getting this mutation type (occur normalized) * the chance of it being in this particular scenario
                        non_like += like
                    else:
                        syn_like += like
    return non_like / syn_like #the ratio.

def do_genebag_analysis(mutdf, cdf, genedf, genepaths, verbose = False, outn = 'geneset_conservation.tsv'):
    '''
    Perform analysis by gene group and save the results to an output dataframe.
    '''
    mutdf['MutType'] = mutdf.Ref + '>' + mutdf.Alt
    if verbose:
        print("Calculating overall occurrence values...", file = sys.stderr)
    atvc = mutdf[(~mutdf.Synonymous) & (~mutdf.Missense) & (mutdf.Ref.isin(['A','T']))].MutType.value_counts().apply(lambda x:x/2)
    gcvc = mutdf[(~mutdf.Synonymous) & (~mutdf.Missense) & (mutdf.Ref.isin(['C','G']))].MutType.value_counts().apply(lambda x:x/2)
    occur = pd.concat([atvc, gcvc])    
    #prepare a dataframe similar to the one in the ontology analysis
    odf = {k:[] for k in ['Term', 'GeneCount', 'MeanGeneLen', 'GeneLenStd', 'mis_pi', 'syn_pi', 'Ratio', 'BaseRatio', 'PermRatio','PermRatioStd']}
    for genepath in genepaths:
        #read in the gene set information
        tid = genepath.split("_")[1]
        if verbose:
            print("Evaluating term {}".format(tid), file = sys.stderr)
        geneset = parse_genes(genepath)
        #get the ratio for this gene set
        brm, brstd = get_ratio(occur, cdf, geneset)
        #get a bootstrap distribution
        perms = point_bootstrap(geneset, genedf, ratio = brm)
        pmean = np.mean(perms)
        pstd = np.mean(perms)
        #assign the values to the dataframe output
        #in this case "term" isn't a classic ontology term but a generic reference to any group of genes id
        #one may notice that this script overall is laid out a bit weird and that these and the ontology group are really similar structures
        #this is true but also this code is copied from more than jupyter notebook and this method requires less reformatting + allows easy differentiation by input
        #subset genedf to the relevant terms.
        sgdf = genedf[genedf.GID.isin(geneset)]
        odf['Term'].append(tid)
        odf['GeneCount'].append(len(geneset))
        odf['MeanGeneLen'].append(np.mean(sgdf.Length))
        odf['GeneLenStd'].append(np.std(sgdf.Length))
        odf['mis_pi'].append(sum(sgdf.mis_pi))
        odf['syn_pi'].append(sum(sgdf.syn_pi))
        odf['Ratio'].append((sum(sgdf.mis_pi)/sum(sgdf.syn_pi)) / brm) #closest to "reality"
        odf['PermRatio'].append(pmean)
        odf['PermRatioStd'].append(psd)
        odf['BaseRatio'].append(brm)
    odf.to_csv(outn, sep = '\t')
    if verbose:
        print("Gene group analyses complete.")

def do_ontology_analysis(mutdf, cdf, genedf, gffdb, verbose = False, outn = 'ontology_conservation.tsv'):
    '''
    Perform ontology analysis and save the results to an output dataframe. 
    '''
    gdb = gffutils.FeatureDB(gffdb)
    mutdf['MutType'] = mutdf.Ref + '>' + mutdf.Alt
    if verbose:
        print("Calculating overall occurrence values...", file = sys.stderr)
    atvc = mutdf[(~mutdf.Synonymous) & (~mutdf.Missense) & (mutdf.Ref.isin(['A','T']))].MutType.value_counts().apply(lambda x:x/2)
    gcvc = mutdf[(~mutdf.Synonymous) & (~mutdf.Missense) & (mutdf.Ref.isin(['C','G']))].MutType.value_counts().apply(lambda x:x/2)
    occur = pd.concat([atvc, gcvc])
    if verbose:
        print("Subsetting by ontology terms...")
    ontd = split_go_term(genedf, gdb)
    print("QC: OntDict Size", len(ontd), file = sys.stdout)
    print("QC: GeneDF Size", genedf.shape, file = sys.stdout)
    odf = {k:[] for k in ['Term', 'GeneCount', 'MeanGeneLen', 'GeneLenStd', 'mis_pi', 'syn_pi', 'Ratio', 'BaseRatio']}
    for ont, sdf in ontd.items():
        odf['Term'].append(ont)
        odf['GeneCount'].append(sdf.shape[0])
        odf['MeanGeneLen'].append(np.mean(sdf.Length))
        odf['GeneLenStd'].append(np.std(sdf.Length))
        odf['mis_pi'].append(sum(sdf.mis_pi))
        odf['syn_pi'].append(sum(sdf.syn_pi))
        brm, brstd = get_ratio(occur, cdf, list(sdf.GID))
        odf['BaseRatio'].append(brm)
        odf['Ratio'].append((sum(sdf.mis_pi)/sum(sdf.syn_pi))/brm)
    odf = pd.DataFrame(odf)
    odf.set_index('Term')
    #get the bootstrap distribution
    #specifically, this is intended to bypass problems with single large genes, containing many mutations, from individually driving effects in specific ontology groups
    if verbose:
        print("Bootstrapping ontology distributions...")
    data = {}
    for term, sdf in tqdm(ontd.items()):
        #make bootstraps
        perms = []
        for i in range(1000):
            sample = sdf.sample(frac=1, replace = True)
            ratio = (sum(sample.mis_pi)/sum(sample.syn_pi))/odf.loc[term]['BaseRatio']
            perms.append(ratio)
        mean = np.mean(perms)
        stdev = np.std(perms)
        perc = percentileofscore(perms, 1)
        data[term] = [mean, stdev, perc]
    odf['BootRatio'] = odf.index.to_series().apply(lambda x:data[x][0])
    odf['BootStd'] = odf.index.to_series().apply(lambda x:data[x][1])
    odf['Perc'] = odf.index.to_series().apply(lambda x:data[x][2])
    print("QC: Ontology Frame Shape", odf.shape, file = sys.stderr)
    odf.to_csv(outn, sep = '\t')
    if verbose:
        print("Saved ontology table.")

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-m', '--mutations', help = 'Path to the mutation dataframe (with snpeff annotations)')
    parser.add_argument('-c', '--codons', help = 'Path to a codon dataframe file (created by gtf_to_lookup)')
    parser.add_argument('-g', '--gffdb', help = 'Path to a gffutils database representing the gene set for this species. If included, performs ontology analysis.', default = None)
    parser.add_argument('-oo', '--ont_out', help = "Name of the file to save the ontology dataframe to, if created. Default is ontology_conservation.tsv")
    parser.add_argument('-go', '--genes_out', help = 'Name of the file to save the gene bag dataframe to, if created. Default is geneset_conservation.tsv')
    parser.add_argument('genebags', nargs = '+', help = 'Path to any number of text files containing gene IDs. Read naming is expected to be species_type_genes.txt (e.g. mel_ribosome_genes.txt)')
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    #get the requisite frame information
    mutdf = pd.read_csv(args.mutations, sep = '\t')
    cdf = pd.read_csv(args.codons, sep = '\t')
    cdf = cdf.set_index("GID")
    if args.verbose:
        print("Building by-gene table", file = sys.stderr)
    genedf = construct_genedf(mutdf, cdf)
    print("QC: gene frame shape", genedf.shape, file = sys.stderr)
    #two steps, can take either or both
    #one is searching through each of the gene bag sets
    #the other is doing ontology analysis with the gffdb
    if args.gffdb != None:
        if args.verbose:
            print("Doing by-ontology analysis", file = sys.stderr)
        do_ontology_analysis(mutdf, cdf, genedf, args.gffdb, args.verbose, args.ont_out)
    if len(args.genebags) > 0:
        if args.verbose:
            print("Analyzing {} gene groups".format(len(args.genebags)), file = sys.stderr)
        do_genebag_analysis(mutdf, cdf, genedf, args.genebags, args.verbose, args.genes_out)
    if args.verbose:
        print("Analysis complete.", file = sys.stderr)

if __name__ == "__main__":
    main()