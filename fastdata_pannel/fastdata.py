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

config = configparser.ConfigParser()
config.read('config.ini')

plotoutDir = config['Control']['plotoutdir']

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

df = df[df.error==0]

df_spec = df.groupby(pd.cut(df['integral'],bins=500)).agg('count').rename(columns={'integral' : 'emission count'}).reset_index()
df_spec = df_spec.rename(columns={'integral' : 'energy'})
df_spec['energy'] = df_spec['energy'].astype('str')
df_spec['energy'] = df_spec['energy'].str.split(',').str[0].str.split('.').str[0].str[1:]

print (df_spec['energy'].head(10))

"""
for i in range (0,8):
    df_channel = df[(df['channel']==i)]
    df_channel = df_channel[['time','integral']]
    df_channel['time'] = df_channel['time'] - 2208988800
    df_channel['time'] = pd.to_datetime(df_channel['time'], unit = 's')
    df_channel = df_channel.set_index('time').resample('10T').count().dropna().reset_index()
    df_channel = df_channel.iloc[1:-1]
    plot = df_channel.plot(x='time',y='integral')
    ax = plt.gca()
    fig = plot.get_figure()
    fig.savefig(plotoutDir+ '/' + 'emission_count_channel'+str(i)+'.png')
"""
