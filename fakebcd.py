from __future__ import print_function, division

from collections import Counter, defaultdict
import itertools
from itertools import izip, cycle
import json
from math import log
import random
from sys import stderr

import screed

mutrate = 1/1000
barcodes = [
    "AAAA",
    "CCCC",
    "GGGG",
    "TTTT",
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
        mm = int(log(random.random()) / log(mutrate))
        yield (barcode, mut(barcode, mm), mm)


def add_barcodes(fname):
    barcode_mm = defaultdict(Counter)
    bcd_ctr = Counter()
    bcd_chain = barcode_chain()
    with screed.open(fname) as reads:
        for read in reads:
            bcd, mutbcd, mismatch = next(bcd_chain)
            bcd_ctr[bcd] += 1
            barcode_mm[bcd][mismatch] += 1
            print("@", read.name, " ", mismatch)
            print(mutbcd + read.sequence)
            print("+")
            print("I" * len(mutbcd) + read.quality)
    for bcd, cnt in bcd_ctr.items():
        print(bcd, "\t", cnt, file=stderr)
    with open('barcode_stats.json', 'w') as fh:
        json.dump(barcode_mm, fh)


if __name__ == "__main__":
    import sys
    add_barcodes(sys.argv[1])
