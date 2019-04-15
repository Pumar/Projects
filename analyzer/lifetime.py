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
from scipy.optimize import curve_fit
import ROOT
import scipy.signal as signal


matplotlib.rcParams.update({'errorbar.capsize': 2})


df_ana = read_root('/home/modulation_cbpf/Documents/analysis/ANA_mx_b_20190323_0305.root', columns=['t0','time','channel','rate','drate','e'])
df_ana['time'] = df_ana['time'] + df_ana['t0'].min() - 2208988800
df_ana = df_ana[['time','channel','rate','drate','e']]

def add_df(df):
    dataframe = read_root('/home/modulation_cbpf/Documents/analysis/'+lines, columns=['t0','time','channel','rate','drate','e'])
    dataframe['time'] = dataframe['time'] + dataframe['t0'].min() - 2208988800
    dataframe = dataframe[['time','channel','rate','drate','e']]
    df = pd.concat([df,dataframe])
    return df

def exp_fit(x, *p):
    a, b = p
    return a*np.exp(-b*x)

def modulated_exp(x, *p):
    A,a,w,phi,l = p
    return A*(1+ a* np.cos(w*x+phi))*np.exp(-l*x)

def split_combine(df):
    df_cs = df[(df['channel'] == 4) | (df['channel'] == 5)][['time','channel','rate','drate']]
    df_ti = df[(df['channel'] == 2) | (df['channel'] == 3)][['time','channel','rate','drate','e']]
    df_co = df[df['channel'] > 5][['time','channel','rate','drate','e']]
    df_cs['peak'] = 'cs'
    ## Fix for rate dropping for cs, should be done better
    df_cs = df_cs[df_cs['rate'] > 70]
    #  Find peak rates for ti
    bins = pd.IntervalIndex.from_tuples([(505,517),(1139,1163),(1662,1680)])
    df_ti = df_ti.set_index(pd.cut(df_ti['e'],bins)).drop(labels='e',axis=1).reset_index()
    df_ti['e'] = df_ti['e'].astype('str')
    df_ti['e'] = np.where(df_ti['e'] == "(505, 517]", 'ti1',np.where(df_ti['e'] == "(1139, 1163]",'ti2','ti3'))
    df_ti = df_ti.rename(columns = {'e' : 'peak'})
    # Find peak rates for cobalt
    bins = pd.IntervalIndex.from_tuples([(1160,1186),(1320,1345),(2495,2520)])
    df_co = df_co.set_index(pd.cut(df_co['e'],bins)).drop(labels='e',axis=1).reset_index()
    df_co['e'] = df_co['e'].astype('str')
    df_co['e'] = np.where(df_co['e'] == "(1160, 1186]", 'co1',np.where(df_co['e'] == "(1320, 1345]",'co2','co3'))
    df_co = df_co.rename(columns = {'e' : 'peak'})
    #Append all split dfs
    return pd.concat([df_ti,df_cs,df_co])


with open('ana_files.txt') as f:
    for lines in f:
        lines = lines.rstrip()
        df_ana = add_df(df_ana)

df_ana['time'] = df_ana['time'] - df_ana['time'].min()
df_ana['time'] = df_ana['time']/(3600*24.)


df_ana = split_combine(df_ana)


for chn in df_ana['channel'].unique():
    df_c = df_ana[df_ana['channel'] == chn]
    for peak in df_c['peak'].unique():
        print(chn,peak)
        df = df_c[(df_c['peak'] == peak)]
        ### Exponential fit & half-life
        p0 = [df['rate'].max(),0.001]
        coeff, var_matrix = curve_fit(exp_fit, df['time'], df['rate'], p0=p0, maxfev = 5000)
        perr = np.sqrt(np.diag(var_matrix)) 
        fit = exp_fit(df['time'],*coeff)
        plt.figure()
        plt.errorbar(df['time'].values,df['rate'].values,yerr=df['drate'].values,fmt='.')
        plt.plot(df['time'].values,fit)
        plt.xlabel('Time [days]')
        plt.ylabel('Activity [Becquerel]')
        plt.text(1,df['rate'].min(),'$\lambda$ = '+str('{:.2e}'.format(coeff[1])) +'\n' + '$T_{1/2}$ = '+str(round(np.log(2)/(coeff[1]*365),2)) + ' $\pm$'+str(round((np.log(2)*perr[1])/(np.power(coeff[1],2)*365),2)) + ' years',bbox=dict(facecolor='red',alpha=0.5))
        plt.savefig('plots/expfit/'+str(peak)+'_expfit_chn'+str(chn)+'.png')
        plt.close()
        ##
        p0 = [df['rate'].max(),0.001,0.1,1,0.001]
        coeff, var_matrix = curve_fit(modulated_exp, df['time'], df['rate'], p0=p0, maxfev = 50000)
        perr = np.sqrt(np.diag(var_matrix))
        print(coeff)
        fit = modulated_exp(df['time'],*coeff)
        plt.figure()
        plt.errorbar(df['time'].values,df['rate'].values,yerr=df['drate'].values,fmt='.')
        plt.plot(df['time'].values,fit)
        plt.xlabel('Time [days]')
        plt.ylabel('Activity [Becquerel]')
        plt.savefig('plots/modulated_exp/'+str(peak)+'_chn'+str(chn)+'.png')
        #
        df_ls = df['time']
        df_ls['rate'] = fit
        df_ls = df.sample(frac=0.5).sort_values(by='time')
        ### lomb scargle
        f = np.linspace(0.01, 15, 1500)
        pgram = signal.lombscargle(df_ls['time'].values, df_ls['rate'].values,f, normalize=True)
        plt.figure()
        plt.plot(f,pgram)
        plt.savefig('plots/lomb_scargle/ls_'+str(peak)+'_'+str(chn)+'.png')
        plt.close()
        
