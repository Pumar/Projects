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

#df = read_root('/home/caio/Documents/processing_data/mx_b_20180716_1742/mx_b_20180716_1742_000000.root', columns=['channel','integral','time','istestpulse','error','baseline','rms','ratio','height'])

if len(rootf)>1:
	framesroot = [read_root(mostRecentDir+'/'+rf, columns=['channel','time','error','integral']) for rf in rootf]
	df = pd.concat(framesroot,axis=0)

else:
	df = read_root(mostRecentDir+'/'+rootf[0], columns=['channel','time','error','integral'])


df = df[(df['channel']==3) & (df['error']==0)]

df = df[['time','integral']]
df['time'] = df['time'] - 2208988800
df['time'] = pd.to_datetime(df['time'], unit = 's')

df = df.set_index('time').resample('10T').count().dropna().reset_index()

plot = df.plot(x='time',y='integral')
ax = plt.gca()
fig = plot.get_figure()
fig.savefig('testfig.png')

