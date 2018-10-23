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


def append_dfs(rootf):
    if len(rootf)>1:
	    framesroot = [read_root(mostRecentDir+'/'+rf, columns=['channel','time','error','integral']) for rf in rootf]
	    df = pd.concat(framesroot,axis=0)
    else:
	    df = read_root(mostRecentDir+'/'+rootf[0], columns=['channel','time','error','integral'])
    return df

dataframe = append_dfs()

df = df[(df.error==0) & (df.integral>0)]

channel_label = ['Cs-137','Cs-137','Co-60','Co-60','Background','Background','Ti-44','Ti-44']

for label,i in zip(channel_label,range (0,8)):
    df_channel = df[(df['channel']==i)][['time','integral']]
    df_spec = df_channel.groupby(pd.cut(df_channel['integral'],bins=500)).agg('count').rename(columns={'integral' : '#emission'}).reset_index()
    df_spec = df_spec.rename(columns={'integral' : 'energy'})
    df_spec['energy'] = df_spec['energy'].astype('str')
    df_spec['energy'] = df_spec['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df_spec = df_spec[['energy','#emission']]
    df_channel['time'] = df_channel['time'] - 2208988800
    df_channel['time'] = pd.to_datetime(df_channel['time'], unit = 's')
    df_channel = df_channel.set_index('time').resample('10T').count().dropna().reset_index().rename(columns={'integral' : '#emission'})
    df_channel = df_channel.iloc[1:-1]
    plot1 = df_channel.plot(x='time',y='#emission', title=label+' - channel '+str(i))
    plt.ylabel('Emissions Count')
    fig1 = plot1.get_figure()
    fig1.savefig(plotoutDir+ '/emissions_' + label+'_channel'+str(i)+'.png')
    plot2 = df_spec.plot(x='energy',y='#emission', title = label+' - channel '+str(i),logy=True)
    plt.ylabel('Emission Count')
    fig2 = plot2.get_figure()
    fig2.savefig(plotoutDir+'/spectrum_'+label+'_channel'+str(i)+'.png')

