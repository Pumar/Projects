import os
import sys
import configparser

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from root_pandas import read_root

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

    
def plot_ana_rates(channel):
    plt.figure()
    df_ana[df_ana.channel==channel].plot(x='time',y='rate',kind='scatter')
    plt.savefig(plotoutDir+'/rate_ana_channel'+str(channel)+'.png')
    plt.close()

def plot_energy(channel):
    plt.figure()
    df_ana[df_ana.channel==channel].plot(x='time',y='e',kind='scatter')
    plt.savefig(plotoutDir+'/energy_ana_channel'+str(channel)+'.png')
    plt.close()

df_ana = read_root(mostRecentDir.split('mx_b')[0] + "analysis/ANA_" + mostRecentDir.split('/')[-2] + '.root', columns = ['rate','drate','time','channel','e'])


for chn in range(0,8):
    if chn == 4 or chn == 5 :
        continue
    plot_ana_rates(chn)
    plot_energy(chn)

