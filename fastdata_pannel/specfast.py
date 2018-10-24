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


def append_dfs():
    if len(rootf)>1:
	    framesroot = [read_root(mostRecentDir+'/'+rf, columns=['channel','time','error','integral']) for rf in rootf]
	    df = pd.concat(framesroot,axis=0)
    else:
	    df = read_root(mostRecentDir+'/'+rootf[0], columns=['channel','time','error','integral'])
    return df

dataframe = append_dfs()


dataframe = dataframe[dataframe.integral>0]

channel_label = ['Cs-137','Cs-137','Co-60','Co-60','Background','Background','Ti-44','Ti-44']

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

def plot_log_spectra(df,channel):
    df_spec = specdf(df[df.error==0],channel)
    df_spec_e = specdf(df[df.error!=0],channel)
    plt.figure()
    ax = df_spec.plot(x='energy',y='#emission',logy=True)
    df_spec_e.plot(x='energy',y='#emission',logy=True, ax = ax)
    ax.legend(['Error free emissions','With error'])
    plt.ylabel('Rate (Hz) / ' +str( ((df_spec['energy'].max()-df_spec['energy'].min())/500) )+' KeV' )
    plt.xlabel('Energy (KeV)')
    plt.savefig(plotoutDir+'/test_spectra_cs.png')

plot_log_spectra(dataframe,0)

"""
for label,i in zip(channel_label,range (0,1)):
    df_channel = df[(df['channel']==i)][['time','integral']]
    df_spec = df_channel.groupby(pd.cut(df_channel['integral'],bins=500)).agg('count').rename(columns={'integral' : '#emission'}).reset_index()
    df_spec = df_spec.rename(columns={'integral' : 'energy'})
    df_spec['energy'] = df_spec['energy'].astype('str')
    df_spec['energy'] = df_spec['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df_spec = df_spec[['energy','#emission']]
    for er in df_spec['error'].unique():
        if er==0:
            dfe = df_spec[df_spec['error']==0]
            plot2 = dfe.plot(x='energy',y='#emission', title = label+' - channel '+str(i),logy=True)
        else:
            dfe = df_spec[df_spec['error']==er]
            dfe.plot(x='energy',y='#emission',logy='True')
    plt.ylabel('Emission Count')
    fig2 = plot2.get_figure()
    fig2.savefig(plotoutDir+'/spectrum_'+label+'_channel'+str(i)+'.png')
"""
