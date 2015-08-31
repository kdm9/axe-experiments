from __future__ import print_function
from sys import stderr
import random
import itertools
from itertools import izip, cycle
from collections import Counter

barcodes = [
    "TATTTTT",
    "TAATA",
    "GTAA",
    "CCGGATAT",
    "ATGCCT",
    "TTCTG",
    "CTAGG",
    "CCTAG",
    "AGGAT",
    "CGCGGAGA",
]

def mut(seq, dist):
    idx = list(range(len(seq)))
    random.shuffle(idx)
    nts = list("AGTC")
    seq = list(seq)
    for i in range(dist):
        rep = random.choice(nts)
        while seq[idx[i]] == rep:
            rep = random.choice(nts)
        seq[idx[i]] = rep
    return "".join(seq)

def barcode_chain():
    for barcode in itertools.cycle(barcodes):
        if random.random() < 0.0005:
            yield (barcode,  mut(barcode, 2))
        elif random.random() < 0.01:
            yield ("None",  "")
        elif random.random() < 0.05:
            yield (barcode, mut(barcode, 1))
        else:
            yield (barcode, barcode)

def fqitr(fp):
    """Fastq iterator. Whoop for python golf!"""
    for h, s, _, q in izip(fp, fp, fp, fp):
        yield (h.strip(), s.strip(), q.strip())

def add_barcodes(fname):
    ifh = open(fname)
    bcd_ctr = Counter()
    bcd_chain = barcode_chain()
    for hdr, seq, qual in fqitr(ifh):
        bcd, mutbcd = next(bcd_chain)
        bcd_ctr[bcd] += 1
        print(hdr)
        print(mutbcd + seq)
        print("+")
        print(mutbcd + qual)
    for bcd, cnt in bcd_ctr.items():
        print(bcd, "\t", cnt, file=stderr)

if __name__ == "__main__":
    import sys
    add_barcodes(sys.argv[1])
