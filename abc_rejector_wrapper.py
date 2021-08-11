import pandas as pd
import numpy as np
import math
from scipy.stats import gamma, expon
import multiprocessing as mp
import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", '--coding_mutations',help = 'Path to a coding mutation table output by -o from genefilter_frame.py', default = 'coding_mutations.tsv')
    # parser.add_argument("-g", '--genmap', help = 'Path to a processed genmap table from a call to genmap with -e 2 and -k 100, intersected with genes and processed downstream.', default = 'gene_genmap.tsv')
    # parser.add_argument("-b", '--badsites', help = 'Path ')
    parser.add_argument('-f', '--frequencies', help = 'Path to the output table from precomputed_spectra.', default = "selinvert_pis.tsv")
    args = parser.parse_args()
    return args

args = argparser()

#https://stackoverflow.com/questions/18622781/why-is-numpy-random-choice-so-slow
def fast_weighted_choice(options,probs):
    x = np.random.rand()
    cum = 0
    for i,p in enumerate(probs):
        cum += p
        if x < cum:
            break
    return options[i]

def nCr(n,r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)

smutdf = pd.read_csv(args.coding_mutations,sep='\t')
#filtering steps.
smutdf = smutdf[smutdf.Effect != 'other'].drop_duplicates(['SSN','Loc',"Effect",'Chro'])
# sumdf = pd.read_csv(args.genmap,sep = '\t')
# goodgids = sumdf[sumdf.MeanScore == 1].GID
# gchro = smutdf.Chro.value_counts().index
# skips = {g:set() for g in gchro}
# with open('badsites.bed') as inf:
#     for entry in inf:
#         chro, loc, locl = entry.strip().split()
#         skips[chro].add(int(loc))
# mask = []
# for i,d in smutdf.iterrows():
#     if d.Loc in skips[d.Chro]:
#         mask.append(True)
#     else:
#         mask.append(False)
# smutdf['Mask'] = mask
# smutdf = smutdf[smutdf.GID.isin(goodgids) & (~smutdf.Mask)]
obs = (sum(smutdf[smutdf.Effect == 'missense'].Pi),sum(smutdf[smutdf.Effect == 'synonymous'].Pi))

def simulator_function_v5(prop_neutral, alpha, mutation_decline = 0.01, size = smutdf.shape[0], pi_mode = True):
    #the simulation process involves a whole lot of statistical drawing
    #for each sample.
    #three summary statistics. Total PiN. Total PiS. Sum of sample frequencies.  
    #leave frequencies out for now since I have no fucking idea what to do with it.
    tpin = 0
    tpis = 0
    detected = 0
    #gens = list(range(min(cdf.generation),max(cdf.generation)+1))
    gens = list(range(2,25))
    declined_gen = [(g**(2*mutation_decline)) for g in gens]
    declined_gen_like = [g/sum(declined_gen) for g in declined_gen]
    gamma_model = gamma(alpha)
    attempt = 0
    while (detected < size):
        attempt += 1
    #for m in [1]: #for testing, fit it on just one mutation to speed things up.
        #first we pick an s based on the input parameters.
        if np.random.rand() < 0.265:
            #this mutation is synonymous
            s = 0
            is_syn = True
        else:
            is_syn = False
            if np.random.rand() < prop_neutral:
                s = 0
            else:
                gs = gamma_model.rvs()
                gs = round(gs,2)
                if gs >= 1:
                    s = .99
                elif gs < 0.0:
                    s = 0
                else:
                    s = gs
        #round rvs to the nearest two decimal places.
        #pick a g based on the relative likelihood of time of mutation.
        #in this case, we're assuming a flat likelihood of mutation across the model, unlike SLiM.
        #g = np.random.choice(list(range(2,25)), p = gen_like)
        g = fast_weighted_choice(gens,declined_gen_like)
        #print(g)
        #with the third iteration, most of the downstream logic is preloaded into the precomputed "spectra"
        #including Pi, which is calculated for each combination generated ahead of simulation time.
        #first, we check whether we're detected at all.
        #if we are, we select a pi value.
        if (g,s) not in pid:
            #if its not in PID, it must never have been detected over 10k generations. Assume this one goes undetected.
            continue
        piv = pid[(g,s)]
        if len(piv) == 0:
            #if its in PID but doesn't have any values saved somehow, skip that also.
            continue
        elif np.random.rand() < detd[(g,s)]:
            #summarize over all the circle 0, which contribute 0 pi. Chance that this combination is not detected.
            continue
        else:
            detected += 1
            #print(attempt,detected)
            if pi_mode:
                lpv = len(piv)
                if (lpv == 1):
                    #print(pid[(g,s)])
                    pi = piv[0]
                else:
                    ind = np.random.randint(0,lpv-1)
                    pi = piv[ind]
                if is_syn:
                    tpis += pi
                else:
                    tpin += pi
            else:
                if is_syn:
                    tpis += 1
                else:
                    tpin += 1
    #convert sfv into relevant summary statistics. 
    #total Pi (sum of circle freq) is a good one to try first.
    return tpin, tpis

cdf = pd.read_csv(args.frequencies,sep='\t')
assert (0 in cdf.selection)
cdf = cdf[cdf.circle > 0]
comtrack = {}
piv = []
for i,d in cdf.iterrows():
    if d.depth in comtrack:
        cd = comtrack[d.depth]
    else:
        cd = nCr(d.depth,2)
        comtrack[d.depth] = cd
    pi = (d.circle * (d.depth - d.circle))/cd
    piv.append(pi)
cdf['Pi'] = piv

cdf['SampleFreq'] = cdf.circle/cdf.depth
detd = {}
pid = {}
for gs,subdf in cdf[cdf.SampleFreq<.1].groupby(['generation','selection']):
    prob_det = subdf.shape[0] / 10000
    detd[gs] = prob_det
    pid[gs] = subdf.Pi.to_numpy()

print('Beginning rejection sampling, obs is ', obs)

def mthread_simulator(params):
    simtpin, simtpis = simulator_function_v5(params[0],params[1],mutation_decline=params[2])
    return simtpin, simtpis, (params[0], params[1], params[2])

def rejection_simulator_pinpis(unimin = 0, unimax = 1, thresh = 1000, target = 100, threads = 22):
    #proportion neutral draws from a uniform (0,1]
    #gamma alpha draws from an exponential
    #decline can be fixed for now. 
    acc_pneu = []
    acc_galpha = []
    acc_mdec = []
    attempts = 0
    try:
        while len(acc_pneu) < target:
            attempts += threads
            #generate threads total parameters to use in this pool run.
            pneu = np.random.uniform(unimin, unimax, size = threads)
            #pneu = .75
            galpha = expon.rvs(0, size = threads)
            mdec = -expon.rvs(0, size = threads)#norm.rvs(0)
            paramset = [(pneu[i],galpha[i],mdec[i]) for i in range(threads)]
            with mp.Pool(threads) as pool:
                for simtpin, simtpis, params in pool.imap_unordered(mthread_simulator, paramset):
                    #print(attempts, pneu, galpha, mdec, simtpin, simtpis)
                    if abs(obs[0]-simtpin) < thresh and abs(obs[1]-simtpis) < thresh:
                        #accepted
                        print('accepted',len(acc_pneu))
                        acc_pneu.append(params[0])
                        acc_galpha.append(params[1])
                        acc_mdec.append(params[2])
    except KeyboardInterrupt:
        pass
    return acc_pneu, acc_galpha, acc_mdec

apn, aga, amd = rejection_simulator_pinpis(unimin = 0.5, thresh = 50, target = 500)
outdf = pd.DataFrame({"PropNeu":apn, "GAlpha":aga, "MDec":amd})
outdf.to_csv("redux_rejection_pests_nmv3.tsv",sep='\t')
print("Completed; found {} acceptable parameter combinations.".format(len(apn)))
