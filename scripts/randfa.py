#!/usr/bin/env python3
from __future__ import print_function, division
import random
import sys

size = int(sys.argv[1])
seed = int(sys.argv[2])

random.seed(seed)

print(">random")
for i in range(size):
    print(random.choice("ACGT"), end='')
    if i % 80 == 79:
        print()
print()
