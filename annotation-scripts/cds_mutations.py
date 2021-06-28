#!/usr/bin/env python3

#import
import argparse
import gffutils
import numpy as np
from sqlite3 import OperationalError
import pandas as pd
from multiprocessing import Pool
import os
#define functions/classes
translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
            'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
            'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G', 
            None:'0'} #if a returned AA has a 0, means the AA at that index is unknown

def parse_variants(path):
    variants = {}
    with open(path) as inf:
        for entry in inf:
            cds, loc, ref, rdep, alts, quals = entry.strip().split()
            if cds not in variants:
                variants[cds] = {"tb":0,"data":[]}
            variants[cds]['tb'] += len([b for b in alts if b in 'ACGT.'])
            for a in alts:
                if a in 'ACGT': #optional additional filtering by quality possible here
                    variants[cds]['data'].append((int(loc), ref + '->' + a))
    return variants

def search_cds(path, variants, gdb, outf = 'ontologies.csv', verbose = False):
    collect = False
    sequence = ''
    ontd = {}
    #print("QC:InitiatingSearch")
    with open(path, 'r') as fastain:
        for entry in fastain:
            if entry[0] == '>':
                spent = entry.strip().split()

                if collect == True: #done collecting from the last one; spent previously defined on earlier loop
                    if sequence[:3] != 'ATG' or translate[sequence[-3:]] != 'Stop':
                        print("Sequence is not proper coding sequence", sequence[:3], sequence[-3:])
                        sequence = ''
                    elif len(sequence) > 3 and len(sequence)%3 == 0:
                        parent = data[7].strip('parent=').split(',')[0]
                        onts = gdb[parent]['Ontology_term']
                        total, sil, non, mis = assign_muts(sequence, variants[data[0][1:]])
                        for o in onts:
                            if o not in ontd:
                                ontd[o] = {"Total Bases":0,'Silent':0,'Nonsense':0,'Missense':0}
                            ontd[o]['Total Bases'] += total
                            ontd[o]['Silent'] += sil
                            ontd[o]['Nonsense'] += non
                            ontd[o]['Missense'] += mis
                        #if verbose:
                            #print("Processed entry with ontology {}".format(onts))
                    else:
                        print("Invalid sequence {} length of {}".format(data[0][1:],len(sequence)))
                    sequence = '' #reset sequence before moving on.

                if spent[0][1:] not in variants:
                    #print('QC:Skipping Entry',spent[0][1:])
                    sequence = ''
                    collect = False
                    continue #no mutations here, keep iterating.
                else: #start collecting this one
                    #print("QC:Recording Entry:")
                    collect = True
                    data = spent

            elif collect:
                sequence += entry.strip()
    #reformat to a dataframe friendly option now that everything is added up
    #convert and save
    #print(variants.keys())
    reform = {t:[] for t in ['Term', 'Total Bases', 'Silent', 'Nonsense', 'Missense']}
    for term,catd in ontd.items():
        reform['Term'].append(term)
        for c,count in catd.items():
            reform[c].append(count)
    reform = pd.DataFrame(reform)
    reform.to_csv(outf)

def get_codons(sequence):
    codons = [sequence[i:i+3] for i in range(0,len(sequence), 3)]
    try:
        assert codons[0] == 'ATG' and translate[codons[-1]] == 'Stop'
    except:
        print("Can't translate codons from", len(sequence), sequence)
        raise AssertionError
    return codons

def assign_muts(sequence, locs):
    #locs is a subset of the Variants dictionary, which has an outer layer of CDS key
    #the locs here is one layer in, with keys 'tb' and 'data' holding a count of all bases see in the circle dataset and the location/transition data in the standard format in 'data'
    #sequence is a raw string, first I have to divide it into exons.
    codons = get_codons(sequence)
    total = locs['tb']
    sil = 0
    non = 0
    mis = 0
    #now process the mutations.
    for l, mtype in locs['data']:
        l = l-1 #possible index issue?
        codind = int((l-l%3)/3) #index of the codon this belongs to
        cind = int(l%3) #index within said codon
        refcodon = codons[codind]
        if 'N' in refcodon:
            continue #can't identify whatever mutation would go here.
        try:
            assert refcodon[cind] == mtype[0]
        except:
            print("Error placing mutation in", refcodon, cind, mtype)
            continue
        newcodon = list(refcodon)
        newcodon[cind] = mtype[-1]
        newcodon = ''.join(newcodon)
        if translate[newcodon] == translate[refcodon]:
            sil += 1
        elif translate[newcodon] == 'Stop' or translate[refcodon] == 'Stop':
            non += 1
        else:
            mis += 1

    return total, sil, non, mis

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-c', '--cds', help = 'Path to a fasta file of CDS sequences, to which the variants are mapped.')
    parser.add_argument('-a', '--annotation', help = 'Path to a genome annotation file in gff3 format.')
    parser.add_argument('-m', '--mutation', help = 'Path to a mpileup file in text format.')
    parser.add_argument('-d', '--database', help = 'Name of the database file to create out of the annotation. If it already exists, will use it instead of overwriting. Default is annotations.db', default = 'annotations.db')
    parser.add_argument('-o', '--csv', help = "Name of the output data table file for visualization and downstream processing. Default is rates.csv", default = 'rates.csv')
    # parser.add_argument('-t', '--threads', type = int, help = 'Number of threads allowed. Default 1', default = 1)
    # parser.add_argument('-n', '--ontologies', nargs = '+', help = 'Name specific ontologies to evaluate. Default evaluates all ontologies', default = None)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    try:
        db = gffutils.create_db(args.annotation, dbfn=args.database, force=False, keep_order=False, merge_strategy='merge', sort_attribute_values=False)
    except OperationalError:
        db = gffutils.FeatureDB(args.database)
    if args.verbose:
        print("Parsing variant file")
    variants = parse_variants(args.mutation)
    if args.verbose:
        print("{} CDSs included in variants".format(len(variants)))
        print("Iterating through CDS")
    search_cds(args.cds, variants, db, outf = args.csv, verbose = args.verbose)

if __name__ == "__main__":
    main()
