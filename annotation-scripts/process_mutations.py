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

def pool_make_exon(pinput):
# for f in db.features_of_type('gene'):
    f, fasta, dbname = pinput
    db = gffutils.FeatureDB(dbname)
    gid = f['ID'][0]
    onts = f['Ontology_term']
    objs = []
    for c in db.children(gid, featuretype = ('CDS')): #these are actually CDS. Could be future issues with overlap of the same mutated exons in many CDSs.
        if c.end - c.start > 3: #can't believe I actually need this conditional.
            seq = c.sequence(fasta)
            obj = Exon(gid, c.chrom, [c.start, c.end], seq, onts, c.frame, c.strand)
            objs.append(obj)
    return objs
    
def create_exons(annotation, fasta, dbname = 'annotations.db', threads = 40, ontologies = None):
    '''
    This function parses a GFF annotation file with the goal of collecting the following attributes for exon spans:
    Gene ID
    Chromosome
    Coordinates
    Ontology
    Sequence 
    
    Then it creates an Exon class object for that entry. This is done for all entries and a list of Exon objects is returned

    Retrieving sequence requires a reference FASTA genome file.
    Relies on GFFUtils and SQLLite.
    '''
    try:
        db = gffutils.create_db(annotation, dbfn=dbname, force = False, keep_order = False, merge_strategy = 'merge', sort_attribute_values = False)
    except OperationalError: #it already exists
        db = gffutils.FeatureDB(dbname)
    p = Pool(threads) #multithread this process because its slow and exons don't need to be in order.
    if ontologies != None:
        if type(ontologies) != list:
            onts = [ontologies]
        else:
            onts = ontologies
        exon_lvec = p.map(pool_make_exon, [(f, fasta, dbname) for f in list(db.features_of_type('gene')) if any([o in f['Ontology_term'] for o in onts])]) #only create exons for entries I care about then.
    else:
        exon_lvec = p.map(pool_make_exon, [(f, fasta, dbname) for f in list(db.features_of_type('gene'))])        
    p.close()
    p.join()
    exonvec = []
    for group in exon_lvec:
        exonvec.extend(group)
    return exonvec

def parse_pileup(path):
    mutations = {}
    with open(path, 'r') as inf:
        for entry in inf:
            spent = entry.strip().split()
            ref = spent[2].upper()
            cir = ''.join([c for c in spent[4].upper() if c in 'ACGTN.'])
            if ref != 'N' and any([c != 'N' for c in cir]): #don't even parse these. Don't want to record a site with all N because its a waste of memory to have a location with 0 depth and no mutations. They can't produce a mutation so they shouldn't contribute to depth
                chro = spent[0].split('_')[-1] #remove the Scf_ and Scf_NODE_ part.
                if chro not in mutations:
                    mutations[chro] = {}
                loc = int(spent[1])
                if loc not in mutations[chro]:
                    mutations[chro][loc] = []
                assert int(spent[3]) == len(cir)
                depth = len([c for c in cir if c != 'N']) #ignoring N entries for depth.
                mutations[chro][loc] = [depth, ref] #depth and reference. 
                for mut in cir:
                    if mut != "N" and mut != '.': #don't both recording these.
                        change = ref + "->" + mut #for now, ignoring lower-case data.
                        mutations[chro][loc].append(change) #resulting list is [depth, reference base, any mutations present..]
    return mutations

class Exon:
    def __init__(self, gene_id, chrom, coords, sequence, ontology = None, frame = 0, strand = '+'):
        self.gene_id = gene_id
        self.chrom = chrom
        self.coords = coords
        self.strand = strand
        self.frame = frame
        self.translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                            'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
                            'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                            'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                            'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                            'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                            'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G', None:'0'} #if a returned AA has a 0, means the AA at that index is unknown
        try:
            # assert abs(coords[1]-coords[0]+1) == len(sequence) #doesn't make any sense if this isn't true. +1 for being in a 0 coordinate system.
            assert len(sequence) > 3
        except:
            # print("Discordance in length and coordinates", coords[0], coords[1], len(sequence))
            print("Impossibly short sequence", coords, len(sequence))
            raise AssertionError
        self.codons = self._process_sequence(sequence, frame)
        self.ontology = ontology
    #two small functions for qc and such.
    def get_aa(self):
        return ''.join([self.translate[codon] for codon in self.codons])
    def get_seq(self):
        return ''.join(self.codons)

    def _process_sequence(self, sequence, frame):
        #split the sequence string into a list of codons.
        #could make this more complex in the future with reading frames.
        codons = [sequence[i:i+3] for i in range(int(frame), len(sequence), 3)]
        #some codons have N's in them, making them untranslatable; remove these entries and replace them with Nones
        filtered = []
        for entry in codons:
            if 'N' in entry:
                filtered.append(None)
            else:
                filtered.append(entry)
        codons = filtered
        if codons[-1] != None and len(codons[-1]) != 3 and len(codons[-1]) != 0:
            print("Uneven number of codons- check reading frame")
            print(self.coords, self.strand, len(codons), len(sequence), len(codons[-1]), int(frame)) #I always have this issue and I don't know why. Revisit this if positions and coding sequences don't behave as expected.
            print("Trimming {}".format(codons[-1]))
            codons = codons[:-1]
            if self.strand == '+': #sequence runs from start to end, update end coordinate by moving it downwards
                self.coords[1] = self.coords[1] - len(codons[-1])
            elif self.strand == '-': #sequence runs from end to start, update start coordinate by moving it upwards
                self.coords[0] = self.coords[0] + len(codons[-1])
        if self.translate[codons[-1]] != 'Stop':
            print("Does not end with stop codon; check sequence")
            if self.translate[codons[-2]] != "Stop":
                print("Neither does one further back. Check sequence for validity")
                self.codons = codons
                print(self.get_aa()[-5:])
                print(self.get_seq()[-15:])
            # if self.translate[codons[0]] == 'Stop':
                # print("Codons may be right direction but wrong order")
            # if self.translate[codons[0][::-1]] == "Stop":
                # print("The sequence appears to be backward!")
        else:
            print("Ends with a stop codon, as is appropriate")
        return codons
    
    def count_sites(self):
        '''
        This important function evaluates self.codons and counts the number of sites which could potentially produce mutations of each type and what those mutations are.
        Agnostic to genome coordinate but not to the type of mutation or where it is in its codon.
        '''
        totals = {'silent':[],'missense':[],'nonsense':[]}
        for codon in self.codons:
            if codon == None: #there was an N when parsing this codon.
                continue
            for cindex in [0,1,2]:
                for altbase in [b for b in 'ACGT' if b != codon[cindex]]:
                    ncodon = list(codon)
                    ncodon[cindex] = altbase
                    ncodon = ''.join(ncodon)
                    mutation = (cindex, codon[cindex] + '->' + altbase)
                    if self.translate[codon] == self.translate[ncodon]:
                        totals['silent'].append(mutation)
                    elif self.translate[codon] == "Stop" or self.translate[ncodon] == 'Stop':
                        totals['nonsense'].append(mutation)
                    else:
                        totals['missense'].append(mutation)
        return totals

    #functions for evaluating a set of mutations
    def process_mutations(self, rawmuts):
        '''
        Input: a list of mutations in [(coord,type)] format. Assumption is that coord will lie within the span, otherwise its ignored with a printed warning message.
        Output: a processed list of mutations in the format [(coord, type, category)] where category is whether it is silent, missense, or nonsense with a reading frame left to right starting at index 0.
        '''
        pmuts = []
        for coord, etype in rawmuts:
            cat = self._check_mut(coord, etype)
            if cat != 'ambig': #if it doesn't land in a codon containing an N and flushed, or the 'alternate allele' isn't an N
                pmuts.append((coord, etype, cat, (coord-min(self.coords))%3))
        return pmuts

    def _check_mut(self, coord, etype):
        '''
        Does the legwork of finding the spot in the sequence and checking what category it is.
        '''
        if coord < min(self.coords) or coord > max(self.coords):
            print("Mutation not in this span; check input")
            return None
        #first need to identify the index of this mutation in the sequence.
        if self.coords[1] > self.coords[0]: #span runs forward.
            index = coord - self.coords[0] #- 1
        else: #span runs backward.
            index = coord - self.coords[1] #- 1
        try:
            assert index >= 0
        except:
            print(coord, index, self.coords)
            raise AssertionError
        #find the index where the codon this mutation belongs to starts.
        refcodon = self.codons[int((index-index%3)/3)] #for example, index 7 will return codon index 2, which is the third in the list of self.codons (0 indexed)
        if refcodon == None:
            return 'ambig'
        refaa = self.translate[refcodon]
        try:
            assert refcodon[index%3] == etype[0] #should be represented properly in etype
        except:
            print("Error index issue")
            print(self.coords, refcodon, index, index%3, etype)
            print(self.codons[int((index-1-(index-1)%3)/3)], self.codons[int((index+1-(index+1)%3)/3)], self.get_seq()[index], self.frame, self.strand)
            # raise AssertionError
            return 'ambig' #ignore this mutation.
        ncodon = list(refcodon)
        ncodon[index%3] = etype[-1]
        ncodon = ''.join(ncodon)
        if 'N' in ncodon: #don't want to introduce these.
            return 'ambig' #ambig gets tossed out by the master function
        naa = self.translate[ncodon]
        if refaa == naa: #its silent
            return 'silent'
        elif naa == 'Stop' or refaa == 'Stop': #truncation or extension
            return 'nonsense'
        else:
            return 'missense' #changed one AA

def pool_ont_rater(pinput): #structure is weird because adding pooling later and lazy about restructuring. Also pool has weird demands.
    term, exonvec, mutations, verbose, output = pinput
    if verbose:
        print("Evaluating term {}".format(term))
    data_d = {k:[] for k in ['Term', '# of CDSs', 'Avg Gene Length', 'Overall Rate', 'silent', 'nonsense', 'missense', 'A','C','G','T',0,1,2]}
    exons = [ex for ex in exonvec if term in ex.ontology]
    if verbose:
        print("Term contains {} coding sequences".format(len(exons)))
    rates = calculate_mutations(mutations, exons, verbose = False)
    #rates is a dictionary with keys 'overall', 'silent...', 'A...', '0...' and values that are mutation rate floats for this set of genes conditioned on the key
    data_d['Term'] = term
    data_d['# of CDSs'] = len(exons)
    data_d['Avg Gene Length'] = np.mean([ex.coords[1] - ex.coords[0] for ex in exonvec])
    data_d['Overall Rate'] = rates['overall']
    for k,v in rates.items():
        if k != 'overall':
            data_d[k] = v
    # return data_d #too much memory load.
    if verbose:
        print("Term {} completed, saving data to temp.txt".format(term))
    with open(output, 'a+') as outf:
        for k,v in data_d.items():
            print(str(k) + ',' + str(v), end = '\t', file = outf) #tab delineation between entries.
        print('', file = outf) #inserts a newline

def examine_ontologies(mutations, exonvec, output = 'rates.csv', threads = 40, verbose = False, ontologies = None):
    #sort the exon vector into ontologies and run calculate_mutations separately for each ontology
    #collect these results and load them into the output dataframe, which will be saved to a csv for importing in a notebook for visualization and statistics.
    if ontologies == None:
        ontset = set()
        for exon in exonvec:
            for ont in exon.ontology:
                ontset.add(ont)
    else:
        ontset = set(ontologies)
    if verbose:
        print('{} ontologies collected, processing separately'.format(len(ontset)))
    alldata_d = {k:[] for k in ['Term', '# of CDSs', 'Avg Gene Length', 'Overall Rate', 'silent', 'nonsense', 'missense', 'A','C','G','T','0','1','2']}
    p = Pool(threads)
    p.map(pool_ont_rater, [(term, exonvec, mutations, verbose, 'temp.txt') for term in list(ontset)]) #unsorted.
    p.close()
    p.join()
    if verbose:
        print("Collecting ontology mutation data")
    with open('temp.txt') as inf:
        for entry in inf:
            spent = entry.strip().split('\t') #should be a randomly sorted series of entries of column,value 
            for cold in spent:
                title, data = cold.split(',')
                alldata_d[title].append(data)
    # for od in ontds:
        # for k,v in od.items():
            # alldata_d[k].extend(v)
    if verbose:
        print("Saving output")
    df = pd.DataFrame(alldata_d)
    df.to_csv(output)
    os.remove('temp.txt') #remove temporary file.

def calculate_mutations(mutations, exonvec, verbose = False):
    #assign the mutations by location to exons.
    #mutations is a dictionary of chromosome with values of dictionaries of locations with values of lists of mutations if any.
    #exonvec is an unsorted list of Exon objects
    #first, sort exonvec by chromosome and then by location in that chromosome, in the order X / 2L / 2R / 3L / 3R
    #and make sure that the chromosome key is appropriately formatted as the above, removing the Scf_ that precedes these in sim for example.
    #this function is built around calculating rates of these values for a single ontology.
    possibles = {p:{b:0 for b in 'ACGT'} for p in [0,1,2]} #category is not a sequence feature and should not be included here
    actuals = {c:{p:{b:0 for b in 'ACGT'} for p in [0,1,2]} for c in ['silent','nonsense','missense']} #layering is somewhat arbitrary
    depths = []
    for chro in ['X','2L','2R','3L','3R']:
        if verbose:
            print("Examining exons in chromosome {}".format(chro))
        #note that this will also remove exons which don't map to the primary arms.
        current_exons = [ex for ex in exonvec if ex.chrom.split('_')[-1] == chro]
        #now iterate through the exons of this chromosome and check whether each base in them has an entry in mutations.
        for exobj in current_exons:
            mset = []
            for coord in range(exobj.coords[0], exobj.coords[1]):
                if coord in mutations[chro]:
                    depth = mutations[chro][coord][0]
                    ref = mutations[chro][coord][1]
                    cindex = (coord - exobj.coords[0])%3 #this may be a problem if the db likes to put the start after the end for reverse strands.
                    possibles[cindex][ref] += depth
                    depths.append(depth)
                    #store mutations
                    for mtype in mutations[chro][coord][2:]:
                        mset.append((coord, mtype))

            #process collected mutations and record                    
            catlist = exobj.process_mutations(mset)
            for coordinate, mtype, cat, cindex in catlist: #coordinate doesn't matter in this case
                actuals[cat][cindex][mtype[0]] += 1
            # #structure of sited is category key, value list of (codon index, mutation type)
            # sited = exobj.count_sites()
            # for cat, locs in sited.items():
            #     for cindex, mtype in locs:
            #         possibles[cindex][mtype[0]] += 1 #agnostic to category, which is not a denominator feature.

    mean_depth = np.median(depths) #okay, its a median and not a mean, whatever. Use it for scaling rate estimation
    propd = {k:[] for k in ['category','position','base','prop']}
    for cat, bp in actuals.items():
        for pos, based in bp.items():
            for base in based.keys():
                # try:
                if possibles[pos][base] != 0:
                    prop = actuals[cat][pos][base]/possibles[pos][base] #additional divide by 3 to scale the first value to per sequence base
                else:
                    continue
#                     prop = actuals[cat][pos][base]/(sum([cd[pos][base] for k,cd in possibles.items()])*depth) #category shouldn't be in the denominator
#                     prop = sum([cd[pos][base] for k,cd in actuals.items()])/(sum([cd[pos][base] for k,cd in possibles.items()])*depth)
                # except:
                    # continue
                propd['category'].append(cat)
                propd['position'].append(pos)
                propd['base'].append(base)
                propd['prop'].append(prop)
#     dfs = {c:pd.DataFrame(d) for c,d in propd.items()} 
    propd = pd.DataFrame(propd)
    #do some summary statistics
    rated = {}
    rated['overall'] = sum([sum(propd.loc[propd['category']==cat]['prop'].tolist())/3/mean_depth for cat in ['silent','missense','nonsense']])
    #return a data vector that will represent the ontology presented here.
    #mutation rates conditional on category
    for cat in ['silent','missense','nonsense']:
#         propd = dfs[cat]
        rated[cat] = sum(propd.loc[propd['category']==cat]['prop'].tolist()) /3 /mean_depth #divide by 3 because 3 categories, gets rid of scaling issue.
    for base in 'ACGT':
        rated[base] = sum(propd.loc[propd['base']==base]['prop'].tolist()) /mean_depth
    for pos in [0,1,2]:
        rated[pos] = sum(propd.loc[propd['position']==pos]['prop'].tolist()) /mean_depth
    return rated

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-r', '--reference', help = 'Path to a genome reference fasta file.')
    parser.add_argument('-a', '--annotation', help = 'Path to a genome annotation file in gff3 format.')
    parser.add_argument('-m', '--mutation', help = 'Path to a simplified pileup file in text format.')
    parser.add_argument('-d', '--database', help = 'Name of the database file to create out of the annotation. If it already exists, will use it instead of overwriting. Default is annotations.db', default = 'annotations.db')
    parser.add_argument('-c', '--csv', help = "Name of the output data table file for visualization and downstream processing. Default is rates.csv", default = 'rates.csv')
    parser.add_argument('-t', '--threads', type = int, help = 'Number of threads allowed. Default 1', default = 1)
    parser.add_argument('-o', '--ontologies', nargs = '+', help = 'Name specific ontologies to evaluate. Default evaluates all ontologies', default = None)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    if args.verbose:
        print("Using {} threads".format(args.threads))
        print("Creating/parsing exon database")
    exonv = create_exons(args.annotation, args.reference, dbname = args.database, threads = args.threads, ontologies = args.ontologies)
    if args.verbose:
        print("Parsing mutation file")
    mutations = parse_pileup(args.mutation)
    if args.verbose:
        print("Examining ontologies")
    examine_ontologies(mutations, exonv, output = args.csv, verbose = args.verbose, threads = args.threads, ontologies = args.ontologies)
    
if __name__ == "__main__":
    main()