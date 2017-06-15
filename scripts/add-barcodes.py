#!/usr/bin/env python3
from __future__ import print_function, division
from collections import Counter
import json
import sys
from sys import stderr, stdout

import docopt
import numpy as np
import screed

try:
    from itertools import zip
except ImportError:
    # Python 3's zip is python2's itertools.zip
    pass


# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


CLI = '''
USAGE:
    add-barcodes.py [options] KEYFILE FASTQ [FASTQ2]

OPTIONS:
    -s SEED         Output file [default: 1234]
    -S              Single end mode (only output R1)
    -r RE_SITE      Restriction enzyme site (inserted between barcode and read)
                    [default: ]
    -o OUTPUT       Output file [default: stdout]
    -y STATS        Output file for statistics
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
    freqs = np.random.gamma(gamma_shape, size=len(samples))
    freqs /= np.sum(freqs)

    return freqs.cumsum()


def mutrate(n, phredscore=25):
    '''Random number from the exponential distrbuiton <= n.'''
    error_rate = 10**(phredscore / 10.)
    x = n + 1
    while x > n:
        x = int(np.random.exponential(1./np.log(error_rate)))
    return x


def mutate(seq, dist, alphabet='ACGT'):
    dist = min(dist, len(seq))

    # distance of 0 is unchanged
    if dist < 1:
        return seq

    idx = np.arange(len(seq))
    np.random.shuffle(idx)
    seq = list(seq)
    for i in range(dist):
        replacement = seq[idx[i]]
        while seq[idx[i]] == replacement:
            replacement = alphabet[np.random.choice(len(alphabet))]
        seq[idx[i]] = replacement
    return "".join(seq)


def add_barcode_to_read(read_pair, samples, cumsum_prob, max_mismatch=0.5,
                        re_site='', gibberish_prob=0., single_only=False):
    r = np.random.uniform()
    idx = np.searchsorted(cumsum_prob, r)

    sample = samples[idx]

    def read_str(read, barcode):
        addded_length = len(barcode) + len(re_site)
        # Prevent bugs when len(bcd) > len(read): Expand read length to be
        # longer than required
        fake_qual = read.quality * (addded_length//len(read.quality) + 1)
        # Then shrink back to what's required
        fake_qual = fake_qual[:addded_length]

        if fake_qual:
            avg_qual = sum(ord(x)-33 for x in fake_qual) / float(len(fake_qual))
        else:
            avg_qual = ord(read.quality[0])-33

        if barcode:
            mismatch = mutrate(int(max_mismatch * len(barcode)),
                               phredscore=avg_qual)
            if np.random.uniform() < gibberish_prob:
                # Mutate barcode fully to create random gibberish as a barcode
                mismatch = len(barcode)
        else:
            mismatch = 0

        if re_site:
            mismatch = mutrate(len(re_site), phredscore=avg_qual)
            re_seq = mutate(re_site, mismatch)
        else:
            re_seq = ''

        bcd_seq = mutate(barcode, mismatch)

        sample_tag = json.dumps({'id': sample['id'],
                                 'barcodes': sample['bcd'],
                                 'mismatches': mismatch})

        seq = '@{}\t{}\n{}{}{}\n+\n{}{}'.format(
                read.name, sample_tag,
                bcd_seq, re_seq, read.sequence,
                fake_qual, read.quality
        )
        return {'seq': seq,
                'id': sample['id'],
                'mm': mismatch}

    r1, r2 = read_pair
    b1, b2 = sample['bcd']

    if single_only:
        return read_str(r1, b1)
    else:
        r1 = read_str(r1, b1)
        r2 = read_str(r2, b2)

        seq = r1['seq'] + '\n' + r2['seq']
        id = r1['id']
        mm = (r1['mm'], r2['mm'])
        return {'seq': seq,
                'id': id,
                'mm': mm}


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


def add_barcodes(fq1, fq2, axe_key, outfile='stdout', gamma_shape=2,
                 re_site='', single_only=False, stats_file=None):
    read_pairs = read_interleaved_or_paired(fq1, fq2)
    samples = read_axe_key(axe_key)

    if outfile == 'stdout':
        outfile = stdout
    else:
        outfile = open(outfile, 'w')

    ids = Counter()

    cumsum_sample_freqs = samplefreq(samples)

    print("CLI is:", sys.argv, file=stderr)
    for i, rp in enumerate(read_pairs):
        if i % 10000 == 0:
            print('\rProcessed {} reads'.format(i), file=stderr, end='')
            stderr.flush()
        reads = add_barcode_to_read(rp, samples, cumsum_sample_freqs,
                                    re_site=re_site, single_only=single_only)
        print(reads['seq'], file=outfile)
        ids[reads['id']] += 1
    print('\rProcessed {} reads. Done!'.format(i + 1), file=stderr)
    outfile.close()

    if stats_file is not None:
        with open(stats_file, 'w') as fh:
            # Ensure samples with no reads are included
            stats = {s['id']: ids[s['id']] for s in samples}
            json.dump(stats, fh)


if __name__ == '__main__':
    opts = docopt.docopt(CLI)
    np.random.seed(int(opts['-s']))
    add_barcodes(opts['FASTQ'],
                 opts['FASTQ2'],
                 opts['KEYFILE'],
                 outfile=opts['-o'],
                 gamma_shape=float(opts['-G']),
                 re_site=opts['-r'],
                 single_only=opts['-S'],
                 stats_file=opts['-y'])
