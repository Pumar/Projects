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

def time_conversion(df):
    df['time'] = df['time'] - 2208988800
    df['time'] = pd.to_datetime(df['time'], unit = 's')
    return df

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
    ax.legend(['Error free emissions','With error'])
    plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
    plt.xlabel('Energy (KeV)')
    plt.savefig(plotoutDir+'/log_spectrum_channel'+str(channel)+'.png')
    plt.close()
    return

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
    plt.savefig(plotoutDir+'/spectrum_mc_channel'+str(channel)+'.png')
    plt.close()
    return

def plot_time_series(df,channel):
    #if chn == 4 or chn == 5:
        #return
    df = df[df.error == 0]
    df = time_conversion(df[df.channel==channel])
    df_ana = read_root(mostRecentDir.split('mx_')[0] + "analysis/ANA_" + mostRecentDir.split('/')[-2] + '.root', columns = ['rate','drate','time','channel','e'])
    df_ana = time_conversion(df_ana[(df_ana.channel==channel) & (df_ana.e < 665) & (df_ana.e > 655)])
    df = df.set_index('time').resample('10T').count().dropna().reset_index().rename(columns={'integral' : 'Count'})
    df = df.iloc[1:-1]
    fig,ax = plt.subplots(nrows=2,ncols=1)
    df.plot(x='time',y='Count',ax=ax[0])
    df_ana.plot(x='time',y='rate',ax=ax[1])
    plt.savefig(plotoutDir + '/time_series_channel'+str(channel)+'.png')
    plt.close()
    return

### Fetches processed and analyzed dataframes from root sources
dataframe = append_dfs()


## Energy filter ( energy > 0 KeV )
dataframe = dataframe[(dataframe.integral > 0) & (dataframe.error == 0) & (dataframe.channel==0)]

dataframe = time_conversion(dataframe[dataframe.channel == 1])
dataframe = dataframe.set_index('time').resample('10T').count().dropna().reset_index().rename(columns={'integral' : 'Count'})
print(dataframe)


"""
for chn in range(0,8):
    plot_spectra(dataframe,chn)
    plot_log_spectra(dataframe,chn)
    plot_time_series(dataframe,chn)
"""


