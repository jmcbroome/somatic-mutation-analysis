import pysam
import sys

def reverse_complement(seq):
    bm = {"A":"T","C":"G","G":"C","T":"A"}
    ret = []
    for b in seq:
        ret.append(bm.get(b.upper(),"N"))
    ret.reverse()
    return "".join(ret)

def identify_variants_duplex(read):
    variants = [] 
    seq = read.get_forward_sequence()
    #for some reason, the duplex consensus output bam isn't correctly tagged, so the get_forward_sequence() yields the reverse sequence and it comes out as nonsense.
    #so we have to manually reverse the damn thing.
    if read.is_reverse:
        seq = reverse_complement(seq)
    pairs = read.get_aligned_pairs(matches_only=True,with_seq=True)
    quals = read.get_forward_qualities()
    for ql, rl, refbase in pairs:
        if ql == None or rl == None or refbase== None:
            continue
        b = seq[ql]
        refbase = refbase.upper()
        if b != 'N' and b != refbase:
            q = quals[ql]
            variants.append((ql, refbase, b, q))
    return variants

bamf = pysam.AlignmentFile(sys.argv[1],"rb")
print(bamf.header,end="")
for r in bamf.fetch():
    toss = False
    try:
        variants = identify_variants_duplex(r)
    except:
        toss = True
        variants = []
    for v in variants:
        readpos = v[0]
        if readpos < 4 or readpos > r.query_length - 4:
            toss = True
    if not toss:
        print(r.to_string())