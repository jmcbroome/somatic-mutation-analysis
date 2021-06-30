# somatic-mutation-analysis
This repository contains scripts, source code, and descriptions of the pipeline for the analysis of mutations identified via the [Circleseq Pipeline](https://github.com/jmcbroome/circleseq). This includes primary results and subdirectories for supplementary analysis. All scripts are written in Python 3 or bash.

## Dependencies
SNPEff 
samtools
Biopython
gffutils

## Full Analysis Procedure

The primary analysis pipeline begins with the output x_sorted.bam files from the [Circleseq Pipeline](https://github.com/jmcbroome/circleseq), which can be easily ran in a batch with run_all.py. Bams are merged into single files representing individuals with samtools merge as needed. The general principle of this procedure is to be as conservative as possible, in order to ensure that all detections are of high-quality and are confidently true somatic mutations for downstream analysis.  

Next, we follow conservative filtering procedures while constructing a pileup. Circleseq consensus bams are flagged as QC pass/fail based on whether any single subsection fails to remap to the target region; reads which fail this are removed. We additionally apply a qc script to conservatively remove reads which contain more than 1 mismatch. 

samtools calmd -ue merged_sorted.bam | python3 qc-utilities/remove_multi_mismatches.py | samtools view -b > nomm.bam

While the calmd bam contains no multiple mismatches at this point, calmd replaces normal sequence strings with N and =. To match the full sequence alignments to the same filters, 
we extract the set of read names we want to retain then refilter the original file accordingly.

samtools view nomm.bam | awk '{print $1}' > nomm_reads.txt

samtools view -h merged_sorted.bam | python3 qc-utilities/remove_reads.py -r nomm_reads.txt | samtools view -b > merged_sorted.nomm.bam

We can proceed to pileup at this point.

samtools mpileup -Q 17 -B -ff QCFAIL -f your_reference.fa merged_sorted.nomm.bam > merged_sorted.nomm.pileup
