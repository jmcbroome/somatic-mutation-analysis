
#from Rachel Mendelson!
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import seaborn as sns
import scipy.stats as stats
import argparse
import subprocess
import sys

# overall template stuff
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output_file', default='read_depths.png')
parser.add_argument('-n', '--num_bins', default=20, help='optional adjustment for number of bins in hist, default=20')
parser.add_argument('-i', '--input_files', default=None, nargs='*', help='reads files formatted with only one depth'
                                                                         'per line if you only want to awk once, '
                                                                         'otherwise reads depths piped in from awk')
parser.add_argument('-f', '--pileup_files', default=None, nargs='*', help='pileup files to automatically read depths from')
args = parser.parse_args()
#
inputFiles = args.input_files
output_file = args.output_file
num_bins = args.num_bins
pileup_files= args.pileup_files

def fileReader(file_name, pileup=False):
    """
    return a list of depths
    """
    depths = []

    if pileup:

        with open(file_name, 'r') as fin:
            for lin in fin:
                if len(lin) <= 3:
                    continue
                depths.append(int(lin.split('\t')[3].strip()))

    else:

        with open(file_name, 'r') as fin:
            for lin in fin:
                depths.append(int(lin.strip()))

    return depths

def main():
    figureWidth = 6
    figureHeight = 11
    plt.clf()
    plt.figure(figsize=(figureWidth, figureHeight))

    depthDict = {}

    # read from files if given
    if inputFiles is not None:
        for i in inputFiles:
            depthDict[i.rsplit('.')[0]] = fileReader(i)

    # otherwise read piped in depths
    elif pileup_files is not None:

        for f in pileup_files:
            depthDict[f.rsplit('.')[0]] = fileReader(f, pileup=True)


    print(depthDict.keys())
    fig, axs = plt.subplots(len(depthDict)+1, figsize=(figureWidth, figureHeight))
    fig.tight_layout()

    # plot each histogram individually
    dict_keys = list(depthDict.keys())
    for i in range(len(depthDict)):

        #set title
        axs[i].set_title(dict_keys[i], fontweight='bold', fontsize=10)
        sns.distplot(depthDict[dict_keys[i]], bins=num_bins, ax=axs[i], fit=stats.gamma)

        axs[i].set_xlim((0, 5000))
        axs[i].set_ylim(0, .0015)

    axs[len(depthDict)].boxplot(depthDict.values(), labels=dict_keys, showfliers=False)

    #plt.show()
    plt.savefig(output_file)


if __name__ == "__main__":
    main()
