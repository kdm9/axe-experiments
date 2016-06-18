#!/usr/bin/env python3

import numpy as np
import pandas as pd
from sys import argv

d = pd.read_table(argv[1], index_col=0)
d = np.asarray(d.as_matrix())
print(100*(d.diagonal().sum()/d.sum()))
