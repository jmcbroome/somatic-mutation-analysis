#!/usr/bin/env python3
#this script implements part of the code from "Exploration Round 2", calculating a DnDs (PiN/PiS) for sites included in the input frame for each annotation column.
import numpy as np
import pandas as pd
import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', help = 'set to True to print status updates', default = True)
    parser.add_argument('-d', '--maxdepth', type = int, help = 'Maximum depth to include in the analysis. Default 2000', default = 2000)
    parser.add_argument('-m', '--mindepth', type = int, help = 'Minimum depth to include in the analysis. Default 0', default = 0)
    parser.add_argument('-r', '--ratio', type = float, help = 'Ratio of S to NS thetas for the target species. Default 3.2', default = 3.2)
    parser.add_argument('-f', '--frame', help = 'Path to the tab-separated mutation dataframe produced by the earlier steps in the pipeline. Required.')
    #parser.add_argument('-g', '--go', help = 'Set a path to a gffdb file created by gffutils to divide mutations by GO term for analysis. Optional', default = None)
    parser.add_argument('-a', '--annotations', nargs = '+', help = 'Any number of columns containing annotation information for calculation. Default is "MyAnn" only.', default = ['MyAnn'])
    parser.add_argument('-p', '--plot_prefix', help = 'Prefix to use for saving KDE Frequency plots. Default no prefix', default = "")
    #parser.add_argument('-s', '--site_count', type = int, nargs = '2', help = 'Two integer values representing the number of synonymous sites and nonsynonymous sites in the genome. Countable from a lookup file. Default is number for simulans reference genome.')
    parser.add_argument('-l', '--lookup', help = 'Path to a lookup file to count sites from for PiS/PiN calculation. Used for all annotation columns.')
    parser.add_argument('-o', '--output', help = 'Name of a file to save the output values to.', default = 'dnds_predicted.tsv')
    args = parser.parse_args()
    return args

def get_counts(path):
    pairs = []
    for a in 'ACGT':
        for b in 'ACGT':
            if a != b:
                pairs.append((a,b))
    def count_sites(lookup):
        sd = {b:0 for b in pairs}
        nd = {b:0 for b in pairs}
        with open(lookup) as inf:
            for entry in inf:
                spent = entry.strip().split('\t')
                if len(spent) != 9:
                    continue
                if spent[0] != 'Chro':
                    if spent[3] != '':
                        for v in spent[3].split(','):
                            if v != '':
                                key = (spent[2],v)
                                if key in sd:
                                    sd[key] += 1
                    if spent[4] != '\n':
                        for v in spent[4].split(','):
                            if v != '':
                                key = (spent[2],v)
                                if key in nd:
                                    nd[key] += 1
        return sd, nd
    sd,nd = count_sites(path)
    sync = sum(sd.values())
    nonc = sum(nd.values())
    return sync, nonc

def calculate_dnds(mutdf, sync = 16672385, nonc = 53473118, col = 'MyAnn'):
    #first for the whole group, then for each individual, save the results in a dictionary of lists.
    data = {k:[] for k in ['SSN', 'Som' + col, 'Germ' + col]}
    sample_germ = mutdf[~mutdf.Somatic][col].value_counts()
    sample_som = mutdf[mutdf.Somatic][col].value_counts()    
    data['SSN'].append('all')
    data['Som' + col].append((sample_som['non']/nonc)/(sample_som['syn']/sync))
    data['Germ' + col].append((sample_germ['non']/nonc)/(sample_germ['syn']/sync))
    #now individuals.
    for ssn in mutdf.SSN.value_counts().index:
        sample_germ = mutdf[(mutdf.SSN == ssn) & (~mutdf.Somatic)][col].value_counts()
        sample_som = mutdf[(mutdf.SSN == ssn) & (mutdf.Somatic)][col].value_counts()
        data['SSN'].append(ssn)
        try:
            gdnds = (sample_germ['non']/nonc)/(sample_germ['syn']/sync)
            data['Germ' + col].append(gdnds)
            sdnds = (sample_som['non']/nonc)/(sample_som['syn']/sync)
            data['Som' + col].append(sdnds)
        except:
            data['Germ' + col].append('fail')
            data['Som' + col].append('fail')
    return data

def get_dnds(frame, lookup, cols = ['MyAnn'], verbose = False):
    outdf = {k:[] for k in ['AnnColumn', 'SSN', 'Som', 'Germ']}
    mutdf = pd.read_csv(frame, sep = '\t')
    sync, nonc = get_counts(lookup)
    for c in cols:
        if verbose:
            print('Calculating with annotation column {}'.format(c))
        cd = calculate_dnds(mutdf, sync, nonc, c)
        #rearrange inputs.
        outdf['AnnColumn'].extend([c for i in range(len(cd['SSN']))])
        outdf['SSN'].extend(cd['SSN'])
        outdf['Som'].extend(cd['Som' + c])
        outdf['Germ'].extend(cd['Germ' + c])
    return pd.DataFrame(outdf)

def main():
    args = argparser()
    df = get_dnds(args.frame, args.lookup, args.annotations, args.verbose)
    df.to_csv(args.output, sep = '\t')

if __name__ == '__main__':
    main()


