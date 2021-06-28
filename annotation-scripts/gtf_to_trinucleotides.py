#!/usr/bin/env python
from Bio import SeqIO as sqio
import pandas as pd
import argparse
import itertools

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--annotation', help = 'path to the gtf')
    parser.add_argument('-g', '--genome', help = 'path to the genome fasta')
    parser.add_argument('-i', '--ids', help = 'file with gene ids to use. Optional', default = None)
    parser.add_argument('-o', '--output', help = 'name of the output file.', default = 'codon_tracker.tsv')
    parser.add_argument('-d', '--detail', type = bool, help = 'toggle to further break down nonsynonymous mutation into conservative and radical changes. Default False', default = False)
    args = parser.parse_args()
    return args

class Coding:
    '''
    This class contains the methods and information required to go from a bundle of split GTF entries to a translateable full CDS
    Used for synonymity analysis and annotation of mutations.
    Dependent on Biopython.
    '''
    translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    aagroups = {'aliphatic':['G','A','V','L','I'], 'hydroxyl':['S','C','U','T','M'], 'cyclic':['P'], 'aromatic':['F','Y','W'], 'basic':['H','K','R'], 'acidic':['D','E','N','Q'], 'N/A':['Stop']}

    def __init__(self, gtfe):
        '''
        Save the location, gene, and strand information for these entries and assert that they match appropriately.
        '''
        self.gid = gtfe[0][-1].split()[1].strip("';\"") #they should all have the same gene on being passed in.
        self.strand = gtfe[0][6] #they should also all have the same strand.
        self.cds = []
        for ge in gtfe:
            #double check their gene and strand attributes match.
            assert gtfe[0][-1].split()[1].strip("';\"") == self.gid
            assert gtfe[0][6] == self.strand
            #proceed to parse it
            if ge[2] == 'start_codon':
                self.start = (ge[0], int(ge[3]), int(ge[4]))
            elif ge[2] == 'stop_codon':
                self.stop = (ge[0], int(ge[3]), int(ge[4]))
            else:
                #need to retain both the location and the frame.
                self.cds.append([ge[0], int(ge[3]), int(ge[4]), int(ge[7])])
    
    def locate(self):
        '''
        Return a tuple representing the overall start and stop location of the sequence
        Uses the original coordinates from the GFF.
        '''
        if self.strand == '+':
            return (self.start[0], self.start[1], self.stop[2])
        elif self.strand == '-':
            return (self.stop[0], self.stop[1], self.start[2])
        
    def extract_sequence(self, fasta):
        '''
        Extract the stitched coding sequence for each attribute in cds, plus the stop codon, from the target fasta.
        Reverse complements as necessary.
        Asserts that it starts with a start and stops with a stop, otherwise will raise an error.
        Returns the stitched coding sequence as a Biopython.SeqRecord object.
        '''
        genome = sqio.to_dict(sqio.parse(fasta, format = 'fasta'))
        if any([e[0] not in genome for e in self.cds]):
            print("Sequence {} is on a chromosome {} not ID'd in the reference".format(self.gid, self.cds[0][0]))
            return "NN" #should be skipped automatically.
        if self.strand == '+':
            seq = self._read_positive(genome)
        elif self.strand == '-':
            seq = self._read_negative(genome)
        #sanity check the sequence
        if seq.seq[:3].upper() != 'ATG':
            print("Sequence {} has an incorrect start {}".format(self.gid, seq.seq[:3]))
        if seq.seq[-3:].upper() not in ['TAA', 'TAG', 'TGA']:
            print('Sequence {} has an incorrect stop {}'.format(self.gid, seq.seq[-3:]))
        if len(seq.seq) %3 != 0:
            print("Sequence {} is not divisible by 3".format(self.gid))
        seq.id = self.gid
        seq.description = 'CDS extracted by Coding custom class'
        seq.name = fasta
        return seq
    
    def _read_positive(self, genome):
        #read the fasta on the positive strand.
        #remember the start codon is included in the cds, but the end codon is not.
        #seqs = [genome[chro][start-1+frame:stop] for chro, start, stop, frame in self.cds]
        seqs = [genome[chro][start-1:stop] for chro, start, stop, frame in self.cds]
        fseq = seqs[0].seq
        for sr in seqs[1:]:
            fseq += sr.seq.upper()
        fseq += genome[self.stop[0]][self.stop[1]-1:self.stop[2]]
        return fseq 
    
    def _read_negative(self, genome):
        rcds = self.cds[::-1]
        #seqs = [genome[chro][start-1:stop-frame] for chro, start, stop, frame in rcds]
        seqs = [genome[chro][start-1:stop] for chro, start, stop, frame in rcds]        
        #seqs = [genome[chro][start-1:stop] for chro, start, stop, frame in self.cds]
        fseq = genome[self.stop[0]][self.stop[1]-1:self.stop[2]]
        for sr in seqs: 
            fseq += sr.seq.upper()
        return fseq.reverse_complement()
    
    def get_index(self):
        '''
        This method returns the vector of all indeces for bases that would go into a stitched sequence, in order that they're read.
        Doesn't need a fasta file. 1-based coordinates.
        '''
        idx = []
        if self.strand == '+':
            #start is included in the first cds, so we can pass it.
            for chro, start, stop, frame in self.cds:
                #for i in range(start-1+frame, stop):
                for i in range(start-1, stop):
                    idx.append(i)
            idx.extend([i for i in range(self.stop[1]-1, self.stop[2])])
        elif self.strand == '-':
            idx.extend([i for i in range(self.stop[1]-1, self.stop[2])])
            rcds = self.cds[::-1]
            for chro, start, stop, frame in rcds:
                for i in range(start-1,stop): #-frame):
                    idx.append(i)
        return idx
    
    def get_index_codons(self, genome):
        fullseq = self.extract_sequence(genome)
        index = self.get_index()
        #zip these vectors together, then divide them into groups of 3.
        combined = list(zip(index, str(fullseq.seq)))
        codons = [combined[i:i+3] for i in range(0, len(index), 3)]
        return codons

    def annotate(self, genome, detail = False):
        '''
        This method uses the GFF to create an annotation structure representing the effects of all possible changes to this coding sequence.
        Returns a dictionary keyed with (index, ref, alt) that returns one of several possible effect codes.
        '''
        codons = self.get_index_codons(genome)
        complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
        effects = {}
        for c in codons:
            intent = self.translate.get(''.join([cv[1] for cv in c]), None)
            if intent == None:
                print("Can't translate original", c)
                continue
            for i,ib in enumerate(c):
                #i represents the location within the codon
                #gi represents the location in the genome
                #b is the actual base.
                gi, b = ib
                for ab in 'ACGT':
                    if ab != b.upper():
                        nc = [cv[1] for cv in c]
                        nc[i] = ab
                        change = self.translate.get(''.join(nc), None)
                        if change == None:
                            print("Can't translate alternative", nc)
                            continue
                        #finally, compare the intent and the change.
                        #have to complement the bases if this gene is negative-strand before saving the data.
                        if self.strand == '+':
                            key = (int(gi)+1, b, ab) #+1 brings it into line with the 1-based indexing of the pileup
                        elif self.strand == '-':
                            key = (int(gi)+1, complement[b], complement[ab]) 

                        if intent == change:
                            effects[key] = 'synonymous' #no change, definitely okay
                        elif intent != 'Stop' and change == 'Stop': #stop gained mid-sequence, very bad
                            effects[key] = 'stop_gained'
                        elif intent == 'Stop' and change != 'Stop': #stop lost at the end, also bad
                            effects[key] = 'stop_lost'
                        else:
                            if detail:
                                #check whether the intent and the change are in the same group.
                                intg = [k for k,v in self.aagroups.items() if intent in v][0] #should be length 1 anyways.
                                if change in self.aagroups[intg]:
                                    #its conservative.
                                    effects[key] = 'conservative'
                                else:
                                    #its radical
                                    effects[key] = 'radical'
                            else:
                                effects[key] = 'nonsynonymous' #different amino acid, maybe bad
        return effects

def group_data(gtf_path, ids = None):
    curgene = None
    curgene_ents = []
    with open(gtf_path) as inf:
        for entry in inf:
            spent = entry.strip().split('\t')
            gid = spent[-1].split()[1].strip("';\"")
            if ids != None:
                if gid not in ids:
                    continue
            #print(gid)
            if curgene == None:
                curgene = gid
            if gid == curgene and spent[2] in ['start_codon', 'CDS', 'stop_codon']:
                curgene_ents.append(spent)
            elif gid != curgene:
                curgene = gid
                if spent[0] in ['NT_033777.3','NT_037436.4','NT_033778.4','NT_033779.5','NC_004354.4'] and len(curgene_ents) > 0:
                    yield curgene_ents #skip scaffolds/chr4, for now at least
                curgene_ents = []
        #yield also at the end of iteration.
        if spent[0] in ['NT_033777.3','NT_037436.4','NT_033778.4','NT_033779.5','NC_004354.4'] and len(curgene_ents) > 0:
            yield curgene_ents #skip scaffolds/chr4, for now at least 
                                
def check_effects(codon, cind):
    #calculate and return a dictionary with key of codon, value of another dictionary with key of effect and value of count of mutations
    #used to add up the number of changes for codon usage.  
    effects = {}
    b = codon[cind]
    intent = Coding.translate.get(codon, None)
    for ab in 'ACGT':
        if ab != b.upper():
            nc = list(codon)
            nc[cind] = ab
            change = Coding.translate.get(''.join(nc), None)
            #finally, compare the intent and the change.
            key = b + ">" + ab #keying the counter on mutation type for this implementation
            effects[key] = {k:0 for k in ['synonymous', 'stop_gained', 'stop_lost', 'nonsynonymous']}
            if intent == change:
                effects[key] = 'synonymous' #no change, definitely okay
            elif intent != 'Stop' and change == 'Stop': #stop gained mid-sequence, very bad
                effects[key] = 'stop_gained'
            elif intent == 'Stop' and change != 'Stop': #stop lost at the end, also bad
                effects[key] = 'stop_lost'
            else:
                effects[key] = 'nonsynonymous'
    return effects

def create_trin_frame(gtf_path, genome, outf = 'codons_tracker.tsv', ids = None, detail = False):
    tv = []
    for trip in [''.join(v) for v in itertools.product('ACGT', repeat = 3)]:
        for ab in 'ACGT':
            if ab != trip[1]:
                tv.append(trip + '>' + ab)

    cdf = {k:[] for k in ["GID", 'Strand', 'Type', 'Count', 'synonymous','stop_gained','stop_lost','nonsynonymous']}
    try:
        for cd in group_data(gtf_path, ids):
            cobj = Coding(cd)
            try:
                if len(cobj.extract_sequence(genome).seq)%3 != 0:
                    continue #skip ones that come out weird for whatever reason.
            except:
                continue
            print("Examining gene {}".format(cobj.gid))

            trinseen = {k:0 for k in [''.join(v) for v in itertools.product('ACGT', repeat = 3)]}
            trintrack = {k:{'nonsynonymous':0, 'synonymous':0, 'stop_gained':0, 'stop_lost':0} for k in tv}
            cseq = cobj.get_index_codons(genome)
            codons = [''.join([sc[1] for sc in c]) for c in cseq]
            #count the number of times this trinucleotide appears in this codon set
            #whenever it does, check the effect of each mutation type on the current codon.
            for ci, cod in enumerate(codons):
                for subi, b in enumerate(cod):
                    b = b.upper()
                    if subi == 0:
                        if ci > 0:
                            prior = codons[ci-1][2].upper()
                            latter = cod[1].upper()
                            thistrin = prior + b + latter
                        else:
                            continue #don't record trinucleotides I don't have context for.
                    elif subi == 1:
                        thistrin = cod.upper() #at the center, the codon is the whole context
                    elif subi == 2:
                        prior = cod[1].upper()
                        if ci < len(codons)-1:
                            latter = codons[ci+1][0].upper()
                        else:
                            continue
                        thistrin = prior + b + latter
                    trinseen[thistrin] += 1
                    #now check mutations at this site
                    ed = check_effects(cod.upper(), subi)
                    for mtype, eff in ed.items():
                        assert mtype[0] == thistrin[1] #middle should be the reference.
                        trintrack[thistrin + '>' + mtype[-1]][eff] += 1

            for trintype, effects in trintrack.items():
                cdf['GID'].append(cobj.gid)
                cdf['Strand'].append(cobj.strand)
                cdf['Type'].append(trintype)
                cdf['Count'].append(trinseen[trintype[:-2]])
                for eff, c in effects.items():
                    cdf[eff].append(c)

    except KeyboardInterrupt:  
        #if I decide to cancel early, I can cntrl+c once and it will immediately save all processed entries.
        #useful if I decide to jump the gun and get an at least partial idea of codon usage bias, if not the complete set.             
        pass
    cdf = pd.DataFrame(cdf)
    #add some more columns of information real quick.
    cdf["Ratio"] = cdf.nonsynonymous / cdf.synonymous #this can vary quite a bit between different genes- absolutely a source of bias for go terms and such.
    cdf.to_csv(outf, sep = '\t')

def read_geneid(path):
    gids = set()
    with open(path) as inf:
        for entry in inf:
            gids.add(entry.strip())
    return gids

def get_paralogs(gids, species = 'sim', path_to_fb_convert = 'dmel_orthologs_in_drosophila_species_fb_2020_02.tsv'):
    #flybase has a particular file they provide that lists all paralogs.
    #this function returns the set of IDs passed in through GIDs.
    pcv = {}
    with open(path_to_fb_convert) as inf:
        for entry in inf:
            if entry[0] != '#':
                spent = entry.strip().split()
                melid = spent[0]
                spec = spent[6].split('\\')[0]
                if spec == 'Dsim' and spent[0] in gids:
                    pcv[melid] = spent[5]
    return pcv

def main():
    args = argparser()
    if args.detail:
        print("using detailed annotation")
    if args.ids != None:
        gids = read_geneid(args.ids)
        #print(len(gids))
        pcv = get_paralogs(gids) #not going to be reusable as is in the future, lazy for now.
        print("Using {} genes".format(len(pcv)))
        create_trin_frame(args.annotation, args.genome, args.output, pcv.values(), detail = args.detail)
    else:
        create_trin_frame(args.annotation, args.genome, args.output, detail = args.detail)

if __name__ == '__main__':
    main()


