import subprocess
import glob
for f in glob.glob('input/*fq.gz'):
    info = f.split('/')[1].split('_')
    name = info[0]
    if info[-1] == '_R1.fq.gz':
        continue #only need it once per pair.
    print('Running sample {}'.format(name))
    subprocess.check_call(['snakemake', '-j', '12', name + '_dmel6_sorted.bam'])
print("All runs complete.")
