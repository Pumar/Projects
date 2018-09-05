import os
import sys
import configparser

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from root_pandas import read_root

### Insert path to most recent directory with slow data
_,mostRecentDir = sys.argv

rootf=[]
for file in os.listdir(mostRecentDir):
	if file.endswith('.root'):
		rootf += [file]

framesroot = [read_root(mostRecentDir+'/'+rf) for rf in rootf]

fastdata = framesroot[0].concat(framesroot[1:],axis=1, ignore_index=True)

print(fastdata.head())
