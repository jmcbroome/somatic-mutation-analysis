# somatic-mutation-analysis
This repository contains scripts, source code, and descriptions of the pipeline for the analysis of mutations identified via the [Circleseq Pipeline](https://github.com/jmcbroome/circleseq). This includes primary results and subdirectories for supplementary analysis. All scripts are written in Python 3 or bash.
Dependent on SNPEff and a number of python packages.

## Pipeline

The primary analysis pipeline begins with the output x_sorted.bam files from the [Circleseq Pipeline](https://github.com/jmcbroome/circleseq), which can be easily ran in a batch with run_all.py.
