import os
import sys
import configparser

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from root_pandas import read_root

### Insert path to most recent directory with fast data
_,mostRecentDir = sys.argv

rootf=[]
for file in os.listdir(mostRecentDir):
	if file.endswith('.root'):
		rootf.append(file)

df = read_root('/home/caio/Documents/processing_data/mx_b_20180716_1742/mx_b_20180716_1742_000000.root', columns=['channel','integral','time','istestpulse','error','baseline','rms','ratio','height'])
df = df[(df['channel']==3) & (df.error==0)]

plot = df.plot.scatter(x='time',y='integral')
ax = plt.gca()
fig = plot.get_figure()
fig.savefig('testfig.png')
"""
if len(rootf)>1:
	framesroot = [read_root(mostRecentDir+'/'+rf, columns=[]) for rf in rootf]
	fastdata = pd.concat(framesroot,axis=1)

else:
	fastdata = read_root(mostRecentDir+'/'+rootf[0])
fastdata['time'] = fastdata['time'] - 2208988800 # Convert from Mac to UNIX time
fastdata['time'] = pd.to_datetime(fastdata['time'], unit = 's')


print(fastdata.columns)
"""
