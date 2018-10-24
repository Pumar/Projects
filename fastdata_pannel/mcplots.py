import os
import sys
import configparser

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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

def get_mc_arr(channel):
    if channel == 0 or channel == 1:
        return pd.read_csv(os.getcwd()+'/monte_carlo/cs137.csv', header = None, names = ['energy','Count']),(580,750)
    elif channel == 2 or channel == 3:
        return pd.read_csv(os.getcwd()+'/monte_carlo/co60.csv', header = None, names = ['energy','Count']),(1050,1450)
    elif channel == 6 or channel == 7:
        return pd.read_csv(os.getcwd()+'/monte_carlo/ti44.csv', header = None, names = ['energy','Count']),(430,1720)

def timedelta(df):
    df['time'] = df['time'] - 2208988800
    df['time'] = pd.to_datetime(df['time'], unit = 's')
    return (df['time'].min()-df['time'].max()).seconds

def specdf(df,channel):
    df = df[(df['channel']==channel)][['time','integral']]
    td = timedelta(df)
    df = df.groupby(pd.cut(df['integral'],bins=500)).agg('count').rename(columns={'integral' : 'Count'}).reset_index()
    df = df.rename(columns={'integral' : 'energy'})
    df['energy'] = df['energy'].astype('str')
    df['energy'] = df['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df = df[['energy','Count']]
    df['Count'] = df['Count']/(td)
    return df

def plot_log_spectra(df,channel):
    df_spec = specdf(df[df.error==0],channel)
    df_spec_e = specdf(df[df.error!=0],channel)
    bins = str((df_spec['energy'].max()-df_spec['energy'].min())/500)
    plt.figure()
    ax = df_spec.plot(x='energy',y='Count',logy=True)
    df_spec_e.plot(x='energy',y='Count',logy=True, ax = ax)
    #inset_ax = inset_axes(ax,  width="50%",height=1.0, loc=1)
    #plot = df_spec.plot(x='energy',y='Count')
    ax.legend(['Error free emissions','With error'])
    plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
    plt.xlabel('Energy (KeV)')
    plt.savefig(plotoutDir+'/spectrum'+str(channel)+'.png')
    plt.close()

def plot_spectra(df,channel):
    if channel == 4 or channel == 5:
        return
    df_spec = specdf(df[df.error==0],channel)
    df_mc,lims = get_mc_arr(channel)
    bins = str((df_spec['energy'].max()-df_spec['energy'].min())/500)
    df_spec = df_spec[(df_spec['energy'] >= lims[0]) & (df_spec['energy'] <= lims[1])]
    df_mc = df_mc[(df_mc['energy'] >= lims[0]) & (df_mc['energy'] <= lims[1])]
    df_mc['Count'] = (df_spec['Count'].max()/df_mc['Count'].max())*df_mc['Count']
    plt.figure()
    ax = df_spec.plot(x='energy',y='Count')
    df_mc.plot(x='energy',y='Count', ax = ax)
    ax.legend(['Experimental Peaks','Simulated Peaks'])
    plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
    plt.xlabel('Energy (KeV)')
    plt.savefig(plotoutDir+'/mc'+str(channel)+'.png')
    plt.close()

dataframe = append_dfs()
dataframe = dataframe[dataframe.integral > 0]

dataframe['time'] = dataframe['time'] - 2208988800
dataframe['time'] = pd.to_datetime(dataframe['time'], unit = 's')

print(dataframe[dataframe.channel==2].head(20))
"""
for i in range(0,8):
    plot_spectra(dataframe,i)
    plot_log_spectra(dataframe,i)
"""



