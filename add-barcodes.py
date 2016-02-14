#!/usr/bin/env python
from __future__ import print_function
import screed
import docopt
import numpy as np
import random
from sys import stderr, stdout
import sys

CLI = '''
USAGE:
    add-barcodes.py [options] KEYFILE FASTQ [FASTQ2]

OPTIONS:
    -o OUTPUT       Output file [default: stdout]
    -G GAMMA_SHAPE  Shape of gamma distribution, source of sample frequencies.
                    [default: 2]
'''


def read_axe_key(path):
    samples = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.lower().startswith('barcode') or \
                    line.startswith('#') or \
                    line.startswith(';'):
                continue
            fields = line.split('\t')
            barcodes = fields[:-1]
            if len(barcodes) != 2:
                barcodes.append('')
            samples[fields[-1]] = barcodes

    return [{'id': x, 'bcd': samples[x]} for x in sorted(samples)]


def samplefreq(samples, gamma_shape=2):
    # Sample frequencies are taken from a gamma distribution, and normalised to
    # sum to 1.
    freqs = np.random.gamma(len(samples), gamma_shape)
    freqs /= np.sum(freqs)

    return freqs.cumsum()


def randexp(n):
    '''Random number from the exponential distrbuiton (base 10) <= n.'''
    x = n + 1
    while x > n:
        x = int(np.random.exponental(np.log(10)))
    return x


def mutate(seq, dist, alphabet='ACGT'):
    dist = min(dist, len(seq))

    # distance of 0 is unchanged
    if dist < 1:
        return seq

    idx = list(range(len(seq)))
    random.shuffle(idx)
    seq = list(seq)
    for i in range(dist):
        replacement = random.choice(alphabet)
        while seq[idx[i]] == replacement:
            replacement = random.choice(nts)
        seq[idx[i]] = replacement
    return "".join(seq)


def add_barcode_to_read(read_pair, samples, cumsum_prob, max_mismatch=0.5):
    idx = np.searchsorted(cumsum_prob, np.random.uniform())

    sample = samples[idx]

    def read_str(read, barcode):
        mismatch = randexp(int(max_mismatch * len(barcode)))

        sample_tag = json.dumps({'id': sample['id'],
                                 'barcodes': sample['bcd'],
                                 'mismatches': mismatch})

        return '@{}\t{}\n{}{}\n+\n{}{}'.format(
                read.name, sample_tag,
                bcd_seq,  read.sequence,
                fake_qual, read.quality
        )

    r1, r2 = read_pair
    b1, b2 = sample['bcd']

    return read_str(r1, b1) + read_str(r2, b2)


def read_interleaved_or_paired(fq1, fq2=None):
    if fq2:
        r1s = screed.open(fq1)
        r2s = screed.open(fq2)
        for r1, r2 in zip(r1s, r2s):
            yield (r1, r2)
    else:
        reads = screed.open(fq1)
        for r1, r2 in zip(reads, reads):
            yield (r1, r2)


def add_barcodes(fq1, fq2, axe_key, outfile='stdout', gamma_shape=2):

    read_pairs = read_interleaved_or_paired(fq1, fq2)
    samples = read_axe_key(axe_key)

    if outfile == 'stdout':
        outfile = stdout
    else:
        outfile = open(outfile, 'w')

    cumsum_sample_freqs = samplefreq(samples)

    print("CLI is: '{}'".format(' '.join(sys.argv)), file=stderr)
    for i, rp in enumerate(read_pairs):
        if i % 10000 == 0:
            print('\rProcessed {} reads'.format(i), file=stderr)
        print(add_barcode_to_read(rp, sample_tag, cumsum_sample_freqs),
              file=outfile)
    outfile.close()


if __name__ == '__main__':
    opts = docopt.docopt(CLI)
    add_barcodes(opts['FASTQ'],
                 opts['FASTQ1'],
                 opts['KEYFILE'],
                 opts['-o'],
                 float(opts['-G']))
