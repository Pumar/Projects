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

energy_bins = 500

source_kev_ranges = (2500,2500,1500,3000)

def append_dfs():
    rootf=[]
    for file in os.listdir(mostRecentDir):
	    if file.endswith('.root'):
		    rootf.append(file)
    if len(rootf)>1:
	    framesroot = [read_root(mostRecentDir+'/'+rf, columns=['channel','time','error','integral']) for rf in rootf]
	    df = pd.concat(framesroot,axis=0)
    else:
	    df = read_root(mostRecentDir+'/'+rootf[0], columns=['channel','time','error','integral'])
    return df

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
    df = df.groupby(pd.cut(df['integral'],bins=energy_bins)).agg('count').rename(columns={'integral' : 'Count'}).reset_index()
    df = df.rename(columns={'integral' : 'energy'})
    df['energy'] = df['energy'].astype('str')
    df['energy'] = df['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df = df[['energy','Count']]
    df['Count'] = df['Count']/(td)
    return df

def plot_log_spectra(df,channel):
    df_spec = specdf(df[(df.error==0) & (df.integral < source_kev_ranges[channel//2])],channel)
    bins = str((df_spec['energy'].max()-df_spec['energy'].min())/energy_bins)
    colors = ['red','purple','orange','yellow','green','violet','gray']
    plt.figure()
    ax = df_spec.plot(x='energy',y='Count',  color = 'k', logy=True)
    for i in df[df['channel'] == channel]['error'].unique():
        df_spec_e = specdf(df[(df.error==i)],channel)
        if i == 0:
            continue
        else:
            df_spec_e.plot(x='energy',y='Count',logy=True, c=colors[i-1] , ax = ax, title = channel_label[channel]+' - channel '+str(channel))
    ax.legend(['Error free emissions','Overflow error','Baseline rms error','Overflow + Baseline', 'Double peak error','Double Peak + Overflow','Double Peak + Baseline','All Errors'])
    plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
    plt.xlabel('Energy (KeV)')
    plt.xlim(0,source_kev_ranges[channel//2])
    plt.savefig(plotoutDir+'/spectra/log_spectrum_channel'+str(channel)+'.png')
    plt.close()
    return

def plot_spectra(df,channel):
    df_spec = specdf(df[(df.error==0) & (df.integral < source_kev_ranges[channel//2])],channel)
    bins = str((df_spec['energy'].max()-df_spec['energy'].min())/energy_bins)
    plt.figure()
    df_spec.plot(x='energy',y='Count', color = 'k')
    plt.ylabel('Rate (Hz) / ' +bins+' KeV' )
    plt.xlabel('Energy (KeV)')
    plt.title(channel_label[channel] +' Spectrum - Channel '+str(channel))
    plt.xlim(0,source_kev_ranges[channel//2])
    plt.legend('')
    plt.savefig(plotoutDir+'/spectra/spectrum_channel'+str(channel)+'.png')
    plt.close()
    return


def plot_rates(channel):
    df_ana = read_root(mostRecentDir.split('mx_')[0] + "analysis/ANA_" + mostRecentDir.split('/')[-2] + '.root', columns = ['rate','drate','time','channel','e'])
    df_ana = df_ana[df_ana.channel==channel]
    df_ana['e'] = np.where( ((df_ana['e'] > 655) & (df_ana['e'] < 666)), 'cs', np.where( ((df_ana['e'] > 1165) & (df_ana['e'] < 1179)
& ((df_ana['channel'] == 6) | (df_ana['channel'] == 7))), 'co1', np.where( ((df_ana['e'] > 1320) & (df_ana['e'] < 1340)),'co2', np.where( ((df_ana['e'] > 2495) & (df_ana['e'] < 2515)), 'co3',np.where( ((df_ana['e'] > 500) & (df_ana['e'] < 522)),'ti1', np.where( ((df_ana['e'] > 1145) & (df_ana['e'] < 1168) & ((df_ana['channel'] == 2) | (df_ana['channel'] == 3))),'ti2' ,np.where( ((df_ana['e'] > 1657) & (df_ana['e'] < 1679) & ((df_ana['channel'] == 2) | (df_ana['channel'] == 3))),'ti3','null' ) ) ) ) ) ) )
    df_ana['colors'] = np.where(( df_ana['e'] == 'cs') | (df_ana['e'] == 'co1') | (df_ana['e'] == 'ti1'),0,np.where( (df_ana['e']== 'co2') | (df_ana['e'] == 'ti2'), 1,2))
    colors = ['black','blue','orange']
    for i in df_ana['colors'].unique():
        df = df_ana[df_ana['colors']==i]
        plt.scatter(x=df['time'].values,y=df['rate'].values,c=colors[i])
        plt.errorbar(x=df['time'].values,y=df['rate'].values,yerr=df['drate'].values,zorder=0, fmt="none",marker="none")     
    plt.xlabel('Time [s]')
    plt.ylabel('Peak Rate [Events in peak/s]')
    #plt.text('Peak rates calculated every '+str()+'s')
    plt.title(channel_label[channel]+' Photopeak rate - channel '+str(channel))
    plt.savefig(plotoutDir + '/activity/rate_channel'+str(channel)+'.png')
    plt.close()
    return

def plot_activity(df,channel):
    ##error free dataframe
    df = df[df.error == 0]
    df = time_conversion(df[df.channel==channel])
    df = df.set_index('time').resample(str(600) +'S').count().dropna().reset_index().rename(columns={'integral' : 'Count'})
    df = df.iloc[1:-1]
    df['error'] = np.sqrt(df['Count'])
    plt.figure()
    plt.plot([],[])
    plt.errorbar(x=df['time'].values, y=df['Count'].values, yerr=df['error'].values, fmt='.',color='k')
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.title(channel_label[channel]+' activity - channel '+str(channel))
    plt.ylabel('Count')
    plt.xlabel('Time')
    plt.savefig(plotoutDir + '/activity/time_series_channel'+str(channel)+'.png')
    plt.close()
    return

def plot_err_activity(df,channel,err):
    df = df[df.error == err]
    df = time_conversion(df[df.channel==channel])
    ##time bin in seconds
    timebin = 600
    df = df.set_index('time').resample(str(timebin) +'S').count().dropna().reset_index().rename(columns={'integral' : 'Count'})
    df = df.iloc[1:-1]
    df['error'] = np.sqrt(df['Count'])
    colors = ['red','purple','orange','yellow','green','violet','gray']
    plt.figure()
    plt.plot([],[])
    plt.errorbar(x=df['time'].values, y=df['Count'].values, yerr=df['error'].values, fmt='.',color=colors[err-1])
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.title(channel_label[channel]+'- Error '+str(err)+' activity - channel '+str(channel))
    plt.ylabel('Count')
    plt.xlabel('Time')
    plt.savefig(plotoutDir + '/err_plots/error_'+str(err)+'_activity_channel'+str(channel)+'.png')
    plt.close()
    return

### Fetches processed dataframe from root source
dataframe = append_dfs()


## Energy filter to mask calibration error( energy > 0 KeV )
dataframe = dataframe[dataframe.integral > 0]


for chn in range(0,8):
    plot_spectra(dataframe,chn)
    plot_log_spectra(dataframe,chn)
    plot_activity(dataframe,chn)
    for err_code in range(1,8):
        plot_err_activity(dataframe,chn,err_code)
    if chn > 1:
        plot_rates(chn)
    print('Channel '+str(chn)+ ' plotted.')

