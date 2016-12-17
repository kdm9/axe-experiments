#!/usr/bin/env python3
from __future__ import print_function
from collections import Counter, defaultdict
import screed
import json
import sys
from os import path
import re
from sys import stderr, stdout

try:
    from itertools import zip
except ImportError:
    pass


def assess_file(filename):
    fq = screed.open(filename)

    ctr = Counter()
    for read in fq:
        rname = read.name
        rname, data = rname.split('\t')
        if '/2' in rname:
            continue
        data = json.loads(data)
        ctr[data['id']] += 1
    return ctr


def calc_sample_demux(table, stats):
    demuxed_to = defaultdict(dict)

    # Re-code to be by source/orig sample then demuxed/assigned sample
    for demux_sample, counts in table.items():
        for orig_sample, count in counts.items():
            demuxed_to[orig_sample][demux_sample] = count

    # Calculate the number of unassigned for each sample
    for sample, counts in demuxed_to.items():
        sum_demuxed = sum(counts.values())
        unassigned = stats[sample] - sum_demuxed
        demuxed_to[sample]['unassigned'] = unassigned

    return demuxed_to


if __name__ == "__main__":
    with open(sys.argv[1]) as fh:
        stats = json.load(fh)
    table = {}
    for fn in sys.argv[2:]:
        samp = re.search(r'([\w-]+).fastq', path.basename(fn))
        if samp is None:
            print("Invalid sample name:", fn, file=stderr)
            sys.exit(1)
        samp = samp.group(1)
        if samp == "unkown":
            continue
        ctr = assess_file(fn)
        table[samp] = ctr

    demux_stats = calc_sample_demux(table, stats)
    samples = list(sorted(stats.keys())) + ['unassigned']

    print('\t' + '\t'.join(samples))
    for orig_sample, counts in sorted(demux_stats.items()):
        print(orig_sample, end='')
        for c in samples:
            print('\t{}'.format(counts.get(c, 0)), end='')
        print()
