from __future__ import print_function, division

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


def lex_keyfile(filename):
    lines = []
    with open(filename) as fh:
        for line in fh:
            for ignored in ["#", ";"]:
                if line.startswith(ignored):
                    break
            else:
                lines.append(line.rstrip().split('\t'))
    return lines


def count_keyfile_lines(filename):
    """Count number of barcodes in keyfile"""
    return len(lex_keyfile(filename))


def keyfile_names(filename):
    """Parse names only from keyfile"""
    return [bcd[-1] for bcd in lex_keyfile(filename)]
