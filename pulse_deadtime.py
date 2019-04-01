import os
import sys
import configparser
import datetime

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.dates as mdates

import pandas as pd
import numpy as np
from root_pandas import read_root
from scipy.optimize import curve_fit

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

def timedelta(df):
    df['time'] = df['time'] - 2208988800
    df['time'] = pd.to_datetime(df['time'], unit = 's')
    return (df['time'].min()-df['time'].max()).seconds

def time_conversion(df):
    df['time'] = df['time'] - 2208988800
    df['time'] = pd.to_datetime(df['time'], unit = 's')
    return df





dataframe = append_dfs()
dataframe['time'] = dataframe['time'].astype(np.float64)
dataframe = dataframe[(dataframe.integral > 0) & (dataframe.channel == 3) & (dataframe.error == 0)][['time','integral']]

dataframe['time'].sort_values()

dataframe['time'] = dataframe['time'] - 2208988800.

dataframe['time'] = dataframe['time'] - dataframe['time'].min()

dataframe = dataframe['time']
#print(dataframe)
def activity_cons(df):
    labels = pd.cut(df,bins = 600, right = True)
    df = df.groupby(labels).count().rename(columns = {'time' : 'count'}).reset_index()
    df = df.rename(columns={ df.columns[1]: "count" })
    df['time'] = df['time'].astype('str')
    df['time'] = df['time'].str.split(',').str[0].astype('str')
    df['time'] = df['time'].str[1:]
    df['time'] = df['time'].astype('float')
    low_act = df.iloc[:300]
    low_act = low_act[low_act['count'] < 12500] 
    reg_act = df.iloc[245:261]
    return low_act['time'].min(), low_act['time'].max(), reg_act['time'].min(), reg_act['time'].max()

low_min,low_max,reg_min,reg_max = activity_cons(dataframe)



df_low = dataframe[( dataframe > low_min ) & ( dataframe < low_max )]
df_reg = dataframe[( dataframe > reg_min ) & ( dataframe < reg_max )]

df_low = df_low*1000000000
df_low = df_low.astype('int')

df_reg = df_reg*1000000000
df_reg = df_reg.astype('int')

#print(df_low.head())
#print(df_reg.head())

low_act_l = [(df_low.iloc[i+1] - df_low.iloc[i]) for i in range(df_low.shape[0]-1)]
#if (df_low.iloc[i+1] - df_low.iloc[i]) < 10000

reg_act_l = [(df_reg.iloc[i+1] - df_reg.iloc[i]) for i in range(df_reg.shape[0]-1)]
#if (df_reg.iloc[i+1] - df_reg.iloc[i]) < 10000

y_low,x_low = np.histogram(low_act_l,bins = 1000)
y_reg,x_reg = np.histogram(reg_act_l[:len(low_act_l)],bins = 1000)

fig, ax = plt.subplots()
ax.scatter(x_low[:-1], y_low,color='k',s=1)
plt.xlabel('Time (ns)')
plt.ylabel('Count')
plt.savefig('/media/sf_Modulation/Modulation_data/plots/deadtime/v3/pulse_deadtime_low_activity.png')
plt.close()

fig, ax = plt.subplots()
ax.scatter(x_reg[:-1], y_reg,color='k',s=1)
plt.xlabel('Time (ns)')
plt.ylabel('Count')
plt.savefig('/media/sf_Modulation/Modulation_data/plots/deadtime/v3/pulse_deadtime_reg_activity.png')
plt.close()


def exp_fit(x, *p):
    ##########################################################
    ### Exponential function
    #
    ### Input:
    #
    ## a, b - fit params
    #
    ##########################################################
    a, b = p
    return a*np.exp(-b*x)

p0 = [2000., 0.01]
coeff, var_matrix = curve_fit(exp_fit,x_low[:-1],y,p0=p0)



"""
plt.figure()
plt.hist(low_act_l, bins = 100)
plt.savefig('/media/sf_Modulation/Modulation_data/plots/deadtime/v2/low_noerr_pulse_deadtime_low_activity.png')
plt.close()

plt.figure()
plt.hist(reg_act_l[:len(low_act_l)], bins = 100)
plt.savefig('/media/sf_Modulation/Modulation_data/plots/deadtime/v2/low_noerr_pulse_deadtime_regular_activity.png')
plt.close()
"""
"""
plt.figure()
dataframe.plot(x='time',y='count')
plt.savefig('/media/sf_Modulation/Modulation_data/plots/deadtime/v2/pulse_deadtime.png')
plt.close()
"""
"""
def low_activity_period(df):
	df = df.set_index('time').resample('10T').count().dropna().reset_index().rename(columns={'integral' : 'Count'})
	df = df.iloc[2:-1]
	df['flag'] = np.where(df['Count']>65000,1,0)
	return df[df['flag']==0]['time'].min(),df[ (df['flag']==0) & (df['time'].dt.hour < 7) ]['time'].max()


def normal_activity_period(df):
	df = df.set_index('time').resample('10T').count().dropna().reset_index().rename(columns={'integral' : 'Count'})
	df = df.iloc[2:-1]
	return df[df['time'].dt.hour > 6]['time'].min(),df[df['time'].dt.hour < 7]['time'].max()

dataframe = time_conversion(dataframe)

low_min,low_max = low_activity_period(dataframe)
norm_min,norm_max = normal_activity_period(dataframe)

dataframe = dataframe['time']

low_act_series = dataframe[(dataframe > low_min) & (dataframe < low_max + datetime.timedelta(minutes = 10))]
low_act_series.sort_values()
###
#print( (low_act_series.iloc[500]-low_act_series.iloc[501]).microseconds )

low_pulse_deltat = [(low_act_series.iloc[i+1]-low_act_series.iloc[i]).microseconds for i in range(0,int(low_act_series.shape[0]-1)) if (low_act_series.iloc[i+1]-low_act_series.iloc[i]).microseconds < 10000]

norm_act_series = dataframe[(dataframe > norm_min) & (dataframe < norm_max + datetime.timedelta(minutes = 10))]
norm_act_series.sort_values()

norm_pulse_deltat = [(norm_act_series.iloc[i+1]-norm_act_series.iloc[i]).microseconds for i in range(0,int(norm_act_series.shape[0]-1)) if (norm_act_series.iloc[i+1]-norm_act_series.iloc[i]).microseconds < 10000]

plt.figure()
plt.hist(low_pulse_deltat[:1000], bins = 200)
plt.savefig('/media/sf_Modulation/Modulation_data/plots/pulse_deadtime_low_activity.png')
plt.close()

plt.figure()
plt.hist(norm_pulse_deltat[:1000], bins = 200)
plt.savefig('/media/sf_Modulation/Modulation_data/plots/pulse_deadtime_regular_activity.png')
plt.close()
"""



