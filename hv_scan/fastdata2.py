import os
import sys
import configparser

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.dates as mdates

import pandas as pd
import numpy as np
from root_pandas import read_root

### Insert path to most recent directory with fast data

_,mostRecentDir = sys.argv

config = configparser.ConfigParser()
config.read('config.ini')

plotoutDir = config['Control']['plotoutdir']

channel_label = ['Background','Background','Ti-44','Ti-44','Cs-137','Cs-137','Co-60','Co-60']

peaks = [[0],
         [0],
         [511,1157,1668],
         [511,1157,1668],
         [661.7],
         [661.7],
         [1173,1332],
         [1173,1332]]

source_kev_ranges = (2000,1500,1500,3000)

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
"""
def get_mc_arr(channel):
    if channel == 4 or channel == 5:
        return pd.read_csv(os.getcwd()+'/monte_carlo/cs137.csv', header = None, names = ['energy','Count']),(580,750)
    elif channel == 6 or channel == 7:
        return pd.read_csv(os.getcwd()+'/monte_carlo/co60.csv', header = None, names = ['energy','Count']),(1050,1450)
    elif channel == 2 or channel == 3:
        return pd.read_csv(os.getcwd()+'/monte_carlo/ti44.csv', header = None, names = ['energy','Count']),(430,1720)
"""

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
    df = df.groupby(pd.cut(df['integral'],bins=600)).agg('count').rename(columns={'integral' : 'Count'}).reset_index()
    df = df.rename(columns={'integral' : 'energy'})
    df['energy'] = df['energy'].astype('str')
    df['energy'] = df['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df = df[['energy','Count']]
    df['Count'] = df['Count']/(td)
    return df

def plot_log_spectra(df,channel):
    df_spec = specdf(df[(df.error==0) & (df.integral < source_kev_ranges[channel//2])],channel)
    bins = str((df_spec['energy'].max()-df_spec['energy'].min())/500)
    colors = ['red','purple','orange']
    plt.figure()
    ax = df_spec.plot(x='energy',y='Count',  color = 'k', logy=True)
    for i in [1,2,4]:
        df_spec_e = specdf(df[df.error==i],channel)
        if i==4:
            df_spec_e.plot(x='energy',y='Count',logy=True, c=colors[i-2] , ax = ax, title = channel_label[channel]+' - channel '+str(channel))
        else:
            df_spec_e.plot(x='energy',y='Count',logy=True, c=colors[i-1] , ax = ax, title = channel_label[channel]+' - channel '+str(channel))
    ax.legend(['Error free emissions','Overflow error','Baseline rms error', 'Double peak error'])
    plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
    plt.xlabel('Energy (KeV)')
    plt.xlim(0,source_kev_ranges[channel//2])
    plt.savefig(plotoutDir+'/log_spectrum_channel'+str(channel)+'.png')
    plt.close()
    return

def plot_spectra(df,channel):
    df_spec = specdf(df[(df.error==0) & (df.integral < source_kev_ranges[channel//2])],channel)
    bins = str((df_spec['energy'].max()-df_spec['energy'].min())/500)
    plt.figure()
    df_spec.plot(x='energy',y='Count', color = 'k')
    plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
    plt.xlabel('Energy (KeV)')
    plt.savefig(plotoutDir+'/spectrum_channel'+str(channel)+'.png')
    plt.close()
    return

"""
def plot_spectra(df,channel):
    df_spec = specdf(df[df.error==0],channel)
    if channel < 2:
        bins = str((df_spec['energy'].max()-df_spec['energy'].min())/500)
        plt.figure()
        df_spec.plot(x='energy',y='Count', kind='scatter', color = 'k')
        plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
        plt.xlabel('Energy (KeV)')
        plt.savefig(plotoutDir+'/spectrum_mc_channel'+str(channel)+'.png')
        plt.close()
    else:
        df_mc,lims = get_mc_arr(channel)
        bins = str((df_spec['energy'].max()-df_spec['energy'].min())/500)
        df_spec = df_spec[(df_spec['energy'] >= lims[0]) & (df_spec['energy'] <= lims[1])]
        df_mc = df_mc[(df_mc['energy'] >= lims[0]) & (df_mc['energy'] <= lims[1])]
        ### Scales the simulation data by peaks to match the experimental peaks
        df_mc['Count'] = (df_spec['Count'].max()/df_mc['Count'].max())*df_mc['Count']
        ###
        plt.figure()
        ax = df_spec.plot(x='energy',y='Count', kind='scatter', color = 'k', title = channel_label[channel]+' - channel '+str(channel))
        df_mc.plot(x='energy',y='Count', ax = ax)
        ax.legend(['Experimental Peaks','Simulated Peaks'])
        plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
        plt.xlabel('Energy (KeV)')
        plt.savefig(plotoutDir+'/spectrum_mc_channel'+str(channel)+'.png')
        plt.close()
        return
"""

def plot_rates(channel):
    df_ana = read_root(mostRecentDir.split('mx_')[0] + "analysis/ANA_" + mostRecentDir.split('/')[-2] + '.root', columns = ['rate','drate','time','channel','e'])
    df_ana['e'] = np.where( ((df_ana['e'] > 655) & (df_ana['e'] < 666)), 'cs', np.where( ((df_ana['e'] > 1165) & (df_ana['e'] < 1179)), 'co1', np.where( ((df_ana['e'] > 1325) & (df_ana['e'] < 1339)),'co2', np.where( ((df_ana['e'] > 2500) & (df_ana['e'] < 2510)), 'co3','null' ) ) ) )
    df_ana['colors'] = np.where(( df_ana['e'] == 'cs') | (df_ana['e'] == 'co1'),0,np.where(df_ana['e']== 'co2', 1,2))
    df_ana = df_ana[df_ana.channel==channel]
    #fig,ax = plt.subplots(nrows=2,ncols=1)
    #plt.plot([],[])
    #plt.errorbar(x=df['time'].values, y=df['Count'].values, yerr=df['error'].values, fmt='.',c='k')
    #df.plot(x='time',y='Count', color="k", marker='.', yerr='error', ax=ax[1])
    colors = ['black','blue','orange']
    plt.errorbar(x=df_ana['time'].values,y=df_ana['rate'].values,yerr=df_ana['drate'].values,zorder=0, fmt="none",marker="none")
    plt.scatter(x=df_ana['time'].values,y=df_ana['rate'].values,c=df_ana['colors'].values, cmap=matplotlib.colors.ListedColormap(colors))
    #, ax=ax[0],fmt='.',yerr
    plt.xlabel('Time [s]')
    plt.ylabel('Peak Rate [Events in peak/s]')
    #plt.text('Peak rates calculated every '+str()+'s')
    plt.title(channel_label[channel]+' - channel '+str(channel))
    plt.savefig(plotoutDir + '/rate_channel'+str(channel)+'.png')
    plt.close()
    return

def plot_activity(df,channel):
    df = df[df.error == 0]
    df = time_conversion(df[df.channel==channel])
    ##time bin in seconds
    timebin = 600
    df = df.set_index('time').resample(str(timebin) +'S').count().dropna().reset_index().rename(columns={'integral' : 'Count'})
    df = df.iloc[1:-1]
    df['error'] = np.sqrt(df['Count'])
    plt.figure()
    #df.plot(x='time',y='Count', color="k", marker='.', yerr='error', title = channel_label[channel]+' - channel '+str(channel) )
    plt.plot([],[])
    plt.errorbar(x=df['time'].values, y=df['Count'].values, yerr=df['error'].values, fmt='.',color='k')
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.title(channel_label[channel]+' - channel '+str(channel))
    plt.ylabel('Count')
    plt.xlabel('Time')
    plt.savefig(plotoutDir + '/time_series_channel'+str(channel)+'.png')
    plt.close()
    return



### Fetches processed dataframe from root source
dataframe = append_dfs()

## Energy filter ( energy > 0 KeV )
dataframe = dataframe[dataframe.integral > 0]


for chn in range(0,8):
    plot_spectra(dataframe,chn)
    plot_log_spectra(dataframe,chn)
    plot_activity(dataframe,chn)
    if chn > 4:
        plot_rates(chn)
    print('Channel '+str(chn)+ ' plotted.')

