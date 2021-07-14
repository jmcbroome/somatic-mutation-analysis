# somatic-mutation-analysis
This repository contains scripts, source code, and descriptions of the pipeline for the analysis of mutations identified via the [Circleseq Pipeline](https://github.com/jmcbroome/circleseq). This includes primary results and subdirectories for supplementary analysis. All scripts are written in Python 3 or bash.

## Dependencies
SNPEff 
samtools
Biopython
gffutils

## Full Analysis Procedure

### Primary Table Construction

The primary analysis pipeline begins with the output x_sorted.bam files from the [Circleseq Pipeline](https://github.com/jmcbroome/circleseq), which can be easily ran in a batch with run_all.py. Bams are merged into single files representing individuals with samtools merge as needed. The general principle of this procedure is to be as conservative as possible, in order to ensure that all detections are of high-quality and are confidently true somatic mutations for downstream analysis.  

Next, we follow conservative filtering procedures while constructing a pileup. Circleseq consensus bams are flagged as QC pass/fail based on whether any single subsection fails to remap to the target region; reads which fail this are removed. We additionally apply a qc script to conservatively remove reads which contain more than 1 mismatch, though these reads may have been already removed upstream.. 

samtools calmd -ue merged_sorted.bam | python3 qc-utilities/remove_multi_mismatches.py | samtools view -b > nomm.bam

While the calmd bam contains no multiple mismatches at this point, calmd replaces normal sequence strings with N and =. To match the full sequence alignments to the same filters, 
we extract the set of read names we want to retain then refilter the original file accordingly.

samtools view nomm.bam | awk '{print $1}' > nomm_reads.txt

samtools view -h merged_sorted.bam | python3 qc-utilities/remove_reads.py -r nomm_reads.txt | samtools view -b > merged_sorted.nomm.bam

We can proceed to pileup at this point.

samtools mpileup -Q 17 -B -ff QCFAIL -f your_reference.fa merged_sorted.nomm.bam > merged_sorted.nomm.pileup

We remove potential PCR duplicates based on the fact that they have overtly similar mapping positions.

python3 qc-utilities/remove_bad_entries.py -i merged_sorted.nomm.pileup -o merged_sorted.nomm.nodup.pileup

We convert the pileup to a vcf format. There's some old code here related to annotating a betabinomial inference of the true body frequency of the data, but this information is not included in the finalized primary analysis. Obtain a vcf header for your genome with samtools or bcftools.

python3 annotation-scripts/pileup_to_annotated_vcf.py -i merged_sorted.nomm.nodup.pileup -a reference_vcf_header.txt -o merged_sorted.nomm.nodup.vcf

The next step is to annotate the vcf with mutation effects, so they can be sorted into coding and noncoding. There are multiple options for this, including applying some custom code included in this repository, but for the purposes of the primary paper analysis I will rely on snpEff. You may need to replace chromosome ID strings to ensure that the snpEff database names match the chromosome names in your reference genome file. 

java -Xmx4g -jar ~/snpEff/snpEff.jar -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -no INTRAGENIC -v your_genome merged_sorted.nomm.nodup.vcf > merged_sorted.nomm.nodup.annotated.vcf

Now we construct a dataframe table for primary analysis, which will include all mutations coding and noncoding. This script expects a specific format of file name for the vcf when parsing multiple files; specifically the second '\_' delimited field is parsed as single letters for strain, stage, and sample number (.e.g something_ra6_something.vcf is parsed as r strain, a stage, number 6 in the output table). 

mv merged_sorted.nomm.nodup.annotated.vcf merged_tt1_sorted.nomm.nodup.annotated.vcf

python3 make_mutation_frame.py -s -o all_mutations.tsv merged_tt1_sorted.nomm.nodup.annotated.vcf

This completes the primary data construction step. The all_mutations.tsv is the base material for all further downstream analysis.

NOTE: Many scripts after this point are not written as flexible software with argument parsing but blocks of code with hard-coded variable and file paths to perform specific analytical tasks. To replicate, you may need to edit the script or your filenames as appropriate.

### Mutation Rates and Base Ratio Calculation

We can immediately apply this table to generate our first figure and calculate the numbers reflected in the first part of our results section.

To do mutation rate calculation, we need to know the overall coverage of each base in each original pileup.

python3 miscellaneous-scripts/count_bases.py < merged_sorted.nomm.nodup.pileup > basecounts.tsv

python3 graph-scripts/noncoding_graphs_statistics.py -m redux_allmuts.tsv -b basecount.tsv -o fig1

The next major step in pursuit of somatic conservation is to establish an estimate for the expected ratio of missense to synonymous mutations across the genome, accounting for both codon usage bias and for the basic rate of different types of mutation. This additionally serves as a branching-off point for the analysis of mutational signatures among somatic mutations. 

To calculate this, first we must calculate genome-wide codon usage. You will need a gtf for your reference genome. 

python3 annotation_scripts/gtf_to_codons.py -g ref.fa -a ref.gtf -o codons.tsv

This table can also be used to establish by-gene and by-groups-of-genes codon usage counts. It also supports gene-level quality filtering of mutations- removing mutations belonging to genes like Myosin Heavy Chain or other outliers in terms of repetitive structure, mutation, or other issues.

Finally we estimate the base ratio value we should use for conservation analysis going forward. This is done by permuting a distribution of ratio values conditioned on real codon usage across the genome and the real set of mutations. This is calculated via permutation of the location of mutations while holding their types constant and checking the set of coding mutations which result.

python3 annotation-scripts/calculate_base_ratio.py -m all_mutations.tsv -r ref.fa -c codons.tsv 

For Drosophila melanogaster under the somatic distribution of mutations, this should be very close to 2.75 with a 95% confidence interval of (2.72,2.77).

#### Supplementary: Somatic Signatures

One research area supported by this type of data is that of mutational signatures, or how specific types of mutations appear in specific trinucleotide contexts. Part of the supplement includes basic mutational signature analysis, with the primary goal of further elucidating the distribution of mutation types over different frequencies of somatic mutation (and therefore over the course of early development). 

To accomplish this, we used the somatic signature analysis tools [Helmsman](https://github.com/carjed/helmsman) and the R package [MutationalPatterns](https://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html).

python3 helmsman_txt_from_frame.py -f all_mutations.tsv > all_mutations_helmsman.txt

python3 helmsman.py -M txt -i all_mutations_helmsman.txt -f ref.fa -p helmsman/

Then load mainscript.R in RStudio and run all to inspect the resulting plots.

### Conservation Analysis

python3 annotation-scripts/genefilter_frame.py -m all_mutations.tsv -o coding_mutations.tsv -g bygene_mutations.tsv -c codons.tsv

The coding_mutations table is used for downstream conservation analysis, as it is contains only missense and synonymous mutations for genes which pass quality filtering.

