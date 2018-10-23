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
"""
_,mostRecentDir = sys.argv
"""
config = configparser.ConfigParser()
config.read('config.ini')

plotoutDir = config['Control']['plotoutdir']
"""
rootf=[]
for file in os.listdir(mostRecentDir):
	if file.endswith('.root'):
		rootf.append(file)

def append_dfs(dirpath):
    if len(rootf)>1:
	    framesroot = [read_root(mostRecentDir+'/'+rf, columns=['channel','time','error','integral']) for rf in dirpath]
	    df = pd.concat(framesroot,axis=0)
    else:
	    df = read_root(mostRecentDir+'/'+rootf[0], columns=['channel','time','error','integral'])
    return df
"""
def get_mc_sims(channel):
    if channel == 0 or channel == 1:
        return read_root(os.getcwd()+'/monte_carlo/MC_cs137_modulation.root')
    elif channel == 2 or channel == 3:
        return read_root(os.getcwd()+'/monte_carlo/MC_co60_modulation.root')
    elif channel == 6 or channel == 7:
        return read_root(os.getcwd()+'/monte_carlo/MC_ti44_modulation.root')
    

#dataframe = append_dfs()

#df = df[(df.error==0) & (df.integral>0)]

def timedelta(df):
    df['time'] = df['time'] - 2208988800
    df['time'] = pd.to_datetime(df['time'], unit = 's')
    return (df['time'].min()-df['time'].max()).seconds


def specdf(df,channel):
    df = df[(df['channel']==channel)][['time','integral']]
    td = timedelta(df)
    df = df.groupby(pd.cut(df['integral'],bins=500)).agg('count').rename(columns={'integral' : '#emission'}).reset_index()
    df = df.rename(columns={'integral' : 'energy'})
    df['energy'] = df['energy'].astype('str')
    df['energy'] = df['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df = df[['energy','#emission']]
    df['#emission'] = df['#emission']/(td)
    return df

def plot_spectra(df,channel):
    df_spec = specdf(df[df.error==0],channel)
    df_spec_e = specdf(df[df.error!=0],channel)
    plt.figure()
    ax = df_spec.plot(x='energy',y='#emission',logy=True)
    df_spec_e.plot(x='energy',y='#emission',logy=True, ax = ax)
    ax.legend(['Error free emissions','With error'])
    plt.ylabel('Rate (Hz) / ' +str( ((df_spec['energy'].max()-df_spec['energy'].min())/500) )+' KeV' )
    plt.savefig(plotoutDir+'/test_spectra_cs.png')

df  = get_mc_sims(0)
plt.figure()
df.plot()
plt.savefig(pltoutDir+'/mc_plot_cs.png')



