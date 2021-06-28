#!/usr/bin/env python3
#this script takes a mark-only dataframe, creates various combinations of filters, and calculates results for DnDs and DFE for each filter combination.
#used for QC purposes.

import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import dadi
from dadi import Godambe
import math
import Optimize_Functions
import argparse
import Selection
from scipy.stats import binom, gamma, percentileofscore
import gffutils
from calculate_dnds_from_frame import get_counts
from calculate_dfe_from_frame import *
from tqdm import tqdm
import itertools

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', help = 'set to True to print status updates', default = True)
    parser.add_argument('-s', '--samples', type = int, help = 'Number of samples to create a virtual SFS for. Default 1000', default = 1000)
    parser.add_argument('-d', '--maxdepth', type = int, help = 'Maximum depth to include in the analysis. Default 2000', default = 2000)
    parser.add_argument('-m', '--mindepth', type = int, help = 'Minimum depth to include in the analysis. Default 0', default = 0)
    parser.add_argument('-c', '--cutoff', type = float, help = 'Exclude sites with this sample frequency and higher. Default 1', default = 1)
    parser.add_argument('-r', '--ratio', type = float, help = 'Ratio of S to NS thetas for the target species. Default 3.2', default = 3.2)
    parser.add_argument('-b', '--bootstraps', type = int, help = 'Number of bootstraps to use for Godambe uncertainty analysis. Default 25', default = 25)
    parser.add_argument('-bs', '--bootstrap_size', type = float, help = 'Size of bootstraps to use. Default .1', default = .1)
    parser.add_argument('-n', '--model_name', help = 'Name to save the optimized demographic information under. Default "demog_opt"', default = "demog_opt")
    parser.add_argument('-f', '--frame', help = 'Path to the mutation dataframe produced by the earlier steps in the pipeline. Must be MARK-ONLY. Required.')
    parser.add_argument('-p', '--prefix', help = 'Prefix to use when saving residual and SFS comparison plots. Default "plot"', default = 'plot')
    parser.add_argument('-t', '--type', help = 'Choose a feature to split the mutation dataframe into for multiple runs. Optional', default = None)
    parser.add_argument('-g', '--go', help = 'Set a path to a gffdb file created by gffutils to divide mutations by GO term for analysis. Optional', default = None)
    parser.add_argument('-e', '--germline', type = bool, help = 'Set to True to infer DFE for germline mutations instead of somatic mutations (cutoff between types set by argument). Default False', default = False)
    parser.add_argument('-l', '--lookup', help = 'Path to a lookup file to count sites from for PiS/PiN calculation. Used for all annotation columns.')
    parser.add_argument('-o', '--output', help = 'Name of a file to save the output values to.', default = 'filter_combinations.tsv')
    parser.add_argument('-a', '--annotations', nargs = '+', help = 'Any number of columns containing annotation information for calculation. Default is "MyAnn" only.', default = ['MyAnn'])
    parser.add_argument('-i', '--filters', nargs = '+', help = 'Names of any number of filter columns to iterate with.')
    parser.add_argument('-u', '--skip_uncert', type = bool, help = 'Set to True to skip calculating Godambe uncertainty. Used for debugging. Default False', default = False)
    args = parser.parse_args()
    return args

def add_dfe_shape(n, a, s):
    values = {}
    example = gamma.rvs(a = a, scale = s, size = 1000)
    for t in range(0,4):
        values[str(10**t)] = percentileofscore(example, 10**t)
        #slightly weird because lazy adaptation from two separate notebook cells.
    values['0'] = n
    raw = [float(values[v]) for v in ['0', '1', '10', '100', '1000']]
    combined = [(raw[0] + value) * (100-raw[0])/100 for value in [0] + raw[1:]]
    each = [combined[0]] + [combined[i]-combined[i-1] for i in list(range(len(combined)))[1:]]
    each.append(100 - combined[-1])
    return {c:each[i] for i,c in enumerate(['0','0-1','1-10','10-100','100-1000','1000+'])}

def compute_sfs(mutdf, coln, args):
    '''
    The first step in the procedure. Uses binomial draws with the predicted sample frequency to estimate the number of inviduals in a virtual population that would have each mutation, then computes the SFS from that distribution.
    '''
    non_sfs = sfs_from_binomial(mutdf, 'non', cutoff = args.cutoff, maxd = args.maxdepth, mind = args.mindepth, samples = args.samples, mode = coln, germ = args.germline)
    syn_sfs = sfs_from_binomial(mutdf, 'syn', cutoff = args.cutoff, maxd = args.maxdepth, mind = args.mindepth, samples = args.samples, mode = coln, germ = args.germline)
    return non_sfs, syn_sfs

def perform_analysis_return(mutdf, prefix, col, args, demog = False):
    non_sfs, syn_sfs = compute_sfs(mutdf, col, args)
    
    if args.verbose:
        print("Fitting demographic parameters...")
    if not demog: #if there isn't already a fit demography available.
        fit_demography(syn_sfs, args)
    #otherwise, it already exists, pull its parameters.
    like, theta, demog_params = get_best_demog('.'.join([str(args.samples), args.model_name, 'optimized', 'txt']))
    if args.verbose:
        print("Bootstrapping for Godambe uncertainty")
    bootstraps = make_bootstrap_sfs_binom(mutdf, 'syn', samples = args.samples) #for Godambe uncertainty of demographic parameters.
    if not args.skip_uncert:
        uncert = Godambe.GIM_uncert(dadi.Numerics.make_extrap_func(exponential_development), [args.samples + 10, args.samples + 20, args.samples + 30], bootstraps, demog_params, syn_sfs)
    else:
        uncert = 'disabled'
    if args.verbose:
        print("Calculating spectra")
    spectra, theta_ns = create_spectra(demog_params, theta, args)
    #not using simple model. complex only
    if args.verbose:
        print("Fitting DFE model")
    complex_popt = fit_neugamma_dfe(spectra, theta_ns, non_sfs)
    #not currently using Godambe uncertainty for dfe parameters, though hypothetically possible, haven't worked out issues with implementation
    if args.verbose:
        print("Collecting final results")
    #save_results(non_sfs, syn_sfs, spectra, simple_popt, complex_popt, like, theta, demog_params, theta_ns, uncert, prefix, args)
    #not plotting residuals within this loop
    #add the dfe shape results based on the content of complex_popt
    n, a, s = complex_popt[1]
    shape_cols = add_dfe_shape(n,a,s)
    resv = [theta, demog_params, uncert, a, s, n] + [shape_cols[i] for i in ['0','0-1','1-10','10-100','100-1000','1000+']]
    return {k:resv[i] for i,k in enumerate(['Theta', 'DemogParams', 'DemogUncert', 'DFE_a', 'DFE_s', 'CompNeuProp', '0','0-1','1-10','10-100','100-1000','1000+'])}

def calculate_dnds(mutdf, sync = 16672385, nonc = 53473118, col = 'MyAnn', skip_ind = False):
    #first for the whole group, then for each individual, save the results in a dictionary of lists.
    data = {k:[] for k in ['SSN/GO', 'SomPiN/PiS', 'GermPiN/PiS']}
    sample_germ = mutdf[~mutdf.Somatic][col].value_counts()
    sample_som = mutdf[mutdf.Somatic][col].value_counts()    
    data['SSN/GO'].append('all')
    data['SomPiN/PiS'].append((sample_som['non']/nonc)/(sample_som['syn']/sync))
    data['GermPiN/PiS'].append((sample_germ['non']/nonc)/(sample_germ['syn']/sync))
    if skip_ind:
        return data
    else:
        #now individuals.
        for ssn in mutdf.SSN.value_counts().index:
            sample_germ = mutdf[(mutdf.SSN == ssn) & (~mutdf.Somatic)][col].value_counts()
            sample_som = mutdf[(mutdf.SSN == ssn) & (mutdf.Somatic)][col].value_counts()
            data['SSN/GO'].append(ssn)
            try:
                gdnds = (sample_germ['non']/nonc)/(sample_germ['syn']/sync)
                data['GermPiN/PiS'].append(gdnds)
                sdnds = (sample_som['non']/nonc)/(sample_som['syn']/sync)
                data['SomPiN/PiS'].append(sdnds)
            except:
                data['GermPiN/PiS'].append('fail')
                data['SomPiN/PiS'].append('fail')
        return data

def do_analysis(subdf, args, sync, nonc):
    '''
    Performs DnDs and DFE analysis on the sub dataframe passed in.
    '''
    result_vectors = []
    for ac in args.annotations:
        current_results = {k:None for k in ['AnnColumn', 'SSN/GO', "SomPiN/PiS", "GermPiN/PiS", 'Theta', 'DemogParams', 'DemogUncert', 'DFE_a', 'DFE_s', 'CompNeuProp', '0','0-1','1-10','10-100','100-1000','1000+']}
        dnds_datad = calculate_dnds(subdf, col = ac) #covers columns 'SSN/GO' through "GermPiN...".
        for i,ssn in enumerate(dnds_datad['SSN/GO']):
            if args.verbose:
                print("Analyzing individual {} with {} mutations".format(ssn, subdf[subdf.SSN == ssn].shape)) #this will say 0 mutations are assigned to all atm, fyi. just a display issue
            dfe_datad = perform_analysis_return(subdf, args.prefix, ac, args)
            #set values in the result vector.
            current_results['AnnColumn'] = ac
            current_results.update({k:v[i] for k,v in dnds_datad.items()})
            current_results.update(dfe_datad) #DFE is a single row, while the DNDs is a NxM set of rows.
            result_vectors.append(current_results)
        #now repeat the above steps for all valid go term subframes. Yeesh.
        gd = split_go_term(subdf, args)
        for term, sdf in gd.items():
            if args.verbose:
                print("Evaluating {} mutations assigned to term {}".format(sdf.shape[0], term))
            current_results = {k:None for k in ['AnnColumn', 'SSN/GO', "SomPiN/PiS", "GermPiN/PiS", 'Theta', 'DemogParams', 'DemogUncert', 'DFE_a', 'DFE_s', 'CompNeuProp', '0','0-1','1-10','10-100','100-1000','1000+']}
            dnds_datad = calculate_dnds(sdf, col = ac, skip_ind = True) #covers columns 'SSN/GO' through "GermPiN...". Skip individuals for go analysis, just care about the whole set.
            #replace the "all" value with the current go term.
            dnds_datad['SSN/GO'] = [term]
            for i,go in enumerate(dnds_datad['SSN/GO']): #same format as above, though these are single entry lists.
                dfe_datad = perform_analysis_return(sdf, args.prefix, ac, args, demog = True)
                #set values in the result vector.
                current_results['AnnColumn'] = ac
                current_results.update({k:v[i] for k,v in dnds_datad.items()})
                current_results.update(dfe_datad) #DFE is a single row, while the DNDs is a NxM set of rows.
                result_vectors.append(current_results)

    return result_vectors #one per row, including annotation columns differnces.

def run_combinations(mutdf, args):
    '''
    Create various combinations of subdataframes using all filter columns denoted in args.
    '''
    results_df = {k:[] for k in args.filters + ['AnnColumn', 'SSN/GO', "SomPiN/PiS", "GermPiN/PiS", 'Theta', 'DemogParams', 'DemogUncert', 'DFE_a', 'DFE_s', 'CompNeuProp', '0','0-1','1-10','10-100','100-1000','1000+']}
    #count the sites for doing dnds out of the lookup outside of the loop.
    sync, nonc = get_counts(args.lookup)
    #first, enumerate combinations of filters, including no filters and all filters.
    filter_combos = [[]]
    for l in range(len(args.filters)):
        filter_combos.extend(itertools.combinations(args.filters, l))
    if any([len(c) == 0 and c != [] for c in filter_combos]):
        print("QC: Removing bad filters detected in ", filter_combos)
        filter_combos = [c for c in filter_combos if len(c) > 0 or c == []]
    print("QC: Will run over filters: ", filter_combos)
    for fc in tqdm(filter_combos):
        if args.verbose:
            print("Analyzing filter set {}".format(fc))
        #get the subframe.
        if args.verbose: 
            print("QC: fc var = {}, mutdf shape {}".format(fc, mutdf.shape))
        subdf = mutdf[mutdf[fc].apply(all, axis = 1)]
        results_vecs = do_analysis(subdf, args, sync, nonc)
        #create the entries and update the results df.
        for rv in results_vecs:
            #for each annotation column used
            for f in args.filters:
                if f in fc:
                    results_df[f].append(True)
                else:
                    results_df[f].append(False)
            for k,v in rv.items():
                results_df[k].append(v) #added results items.
    #finally finished!
    return pd.DataFrame(results_df)

def main():
    args = argparser()
    #read in the mark-only dataframe.
    mutdf = pd.read_csv(args.frame, sep = '\t')
    results = run_combinations(mutdf, args)
    results.to_csv(args.output, sep = '\t')

if __name__ == '__main__':
    main()