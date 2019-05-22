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
from scipy import stats
from scipy.stats import linregress
import ROOT
import scipy.signal as signal
from statsmodels.stats.diagnostic import kstest_normal
from statsmodels.stats.gof import gof_chisquare_discrete
from statsmodels.stats.stattools import durbin_watson
from statsmodels.stats.diagnostic import kstest_normal

matplotlib.rcParams.update({'errorbar.capsize': 2})


df_ana = read_root('/home/modulation_cbpf/Documents/12h_analysis/ANA_mx_b_20190323_0305.root', columns=['t0','time','channel','rate','drate','e'])
df_ana['time'] = df_ana['time'] + df_ana['t0'].min() - 2208988800
df_ana = df_ana[['time','channel','rate','drate','e']]

def add_df(df):
    dataframe = read_root('/home/modulation_cbpf/Documents/12h_analysis/'+lines, columns=['t0','time','channel','rate','drate','e'])
    dataframe['time'] = dataframe['time'] + dataframe['t0'].min() - 2208988800
    dataframe = dataframe[['time','channel','rate','drate','e']]
    df = pd.concat([df,dataframe])
    return df

def exp_fit(x, *p):
    a, b = p
    return a*np.exp(-b*x)

def linear_fit(x,*p):
    a, b = p
    return a-b*x

def double_linear(x,*p):
    a,b,c = p
    return a +b*x + c*x**2

def modulated_exp(x, *p):
    A,a,w,phi,l = p
    return A*(1+ a* np.cos(w*x+phi))*np.exp(-l*x)

def quadrature(x):
    return np.sqrt(np.sum(np.power(x,2)))

def split_combine(df):
    df_cs = df[(df['channel'] == 4) | (df['channel'] == 5)][['time','channel','rate','drate']]
    df_ti = df[(df['channel'] == 2) | (df['channel'] == 3)][['time','channel','rate','drate','e']]
    df_co = df[df['channel'] > 5][['time','channel','rate','drate','e']]
    df_cs['peak'] = 'cs'
    ## Fix for rate dropping for cs, should be done better 
    if df_cs.shape[0] > 0: 
        df_cs = df_cs.groupby(['time','peak']).agg({'channel' : 'min', 'rate' : 'sum', 'drate' : quadrature})
        df_cs.reset_index(inplace=True)
        df_cs.rename(columns = {'peak' : 'source'},inplace=True)
        df_cs.drop(labels='channel',axis=1,inplace=True)
        df_cs = df_cs[df_cs['rate'] > 160]
    #  Find peak rates for ti
    if df_ti.shape[0] > 0:
        bins = pd.IntervalIndex.from_tuples([(505,517),(1139,1163),(1662,1680)])
        df_ti = df_ti.set_index(pd.cut(df_ti['e'],bins)).drop(labels='e',axis=1).reset_index()
        df_ti['e'] = df_ti['e'].astype('str')
        df_ti['e'] = np.where(df_ti['e'] == "(505, 517]", 'ti1',np.where(df_ti['e'] == "(1139, 1163]",'ti2','ti3'))
        df_ti = df_ti.rename(columns = {'e' : 'peak'})
        df_ti = df_ti.groupby(['time','peak']).agg({'channel' : 'min', 'rate' : 'sum', 'drate' : quadrature})
        df_ti.reset_index(inplace=True)
        df_ti.drop(labels='channel',axis=1,inplace=True)
        df_ti['peak'] = df_ti['peak'].str[:-1]
        df_ti.rename(columns = {'peak' : 'source'}, inplace=True)
        df_ti = df_ti.groupby('time').agg({'source' : 'min', 'rate' : 'sum', 'drate' : quadrature})
        df_ti.reset_index(inplace=True)
    # Find peak rates for cobalt
    if df_co.shape[0] > 0:
        bins = pd.IntervalIndex.from_tuples([(1160,1186),(1320,1345),(2495,2520)])
        df_co = df_co.set_index(pd.cut(df_co['e'],bins)).drop(labels='e',axis=1).reset_index()
        df_co['e'] = df_co['e'].astype('str')
        df_co['e'] = np.where(df_co['e'] == "(1160, 1186]", 'co1',np.where(df_co['e'] == "(1320, 1345]",'co2','co3'))
        df_co = df_co.rename(columns = {'e' : 'peak'})
        df_co = df_co.groupby(['time','peak']).agg({'channel' : 'min', 'rate' : 'sum', 'drate' : quadrature})
        df_co.reset_index(inplace=True)
        df_co.drop(labels='channel',axis=1,inplace=True)
        df_co['peak'] = df_co['peak'].str[:-1]
        df_co.rename(columns = {'peak' : 'source'}, inplace=True)
        df_co = df_co.groupby('time').agg({'source' : 'min', 'rate' : 'sum', 'drate' : quadrature})
        df_co.reset_index(inplace=True)
    #Append all split dfs
    return pd.concat([df_ti,df_cs,df_co])


with open('ana_files.txt') as f:
    for lines in f:
        lines = lines.rstrip()
        df_ana = add_df(df_ana)

df_ana['time'] = df_ana['time'] - df_ana['time'].min()
df_ana['time'] = df_ana['time']/(3600*24.)


#df_ana = split_combine(df_ana)


df_co6 = split_combine(df_ana[df_ana.channel == 6])
df_co7 = split_combine(df_ana[df_ana.channel == 7])
df_co7 = df_co7.drop(labels=['channel','e','peak'],axis=1)
df_co6 = df_co6.drop(labels=['channel','e','peak'],axis=1)


df = df_co6
p0 = [df['rate'].max(),0.001]
coeff, var_matrix = curve_fit(exp_fit, df['time'], df['rate'], p0=p0, maxfev = 50000)
perr = np.sqrt(np.diag(var_matrix))
x = df['time'].values
y = df['rate'].values
y_err = df['drate']

#
p0 = [max(y),-1]
coeff2, var_matrix2 = curve_fit(linear_fit, x, y, p0=p0, maxfev = 50000)
sl, inter, rv, p, stderr = linregress(x=x,y=y)
#
e_fit = exp_fit(x,*coeff)
"""
slope = coeff[1]
intercept = np.log(coeff[0])
coeff = [intercept,slope]
df['lin'] = np.vectorize(linear_fit)(df['time'].values,*coeff)
df['residue'] = df['lin'] - y
df['exp_residue'] = e_fit-df['rate']
r_squared = 1 - (np.sum(np.power(df['residue'],2)))/(np.sum(np.power(y-np.mean(y),2)))
print(r_squared)
fit = linear_fit(x, *coeff)
"""
plt.figure()
plt.errorbar(x,y,yerr=y_err,fmt='.')
plt.plot(x,e_fit)
#
plt.title('The activity of the $^{60}$Co source, channel 6')
#plt.text(0,3.87,'$R^2 = $'+str(round(r_squared,4)),bbox=dict(facecolor='red',alpha=0.5))
plt.text(0,47.82,'$\lambda$ = '+str('{:.2e}'.format(coeff[1])) +'\n' + '$T_{1/2}$ = '+str(round(np.log(2)/(coeff[1]*365),2)) + ' $\pm$'+str(round((np.log(2)*perr[1])/(np.power(coeff[1],2)*365),2)) + ' years',bbox=dict(facecolor='red',alpha=0.5))
plt.ylabel('Activity [Becquerel]')
plt.xlabel('Time (days)')
#
plt.savefig('plots/agg_peaks/co6.png')
plt.close()
"""
plt.figure()
plt.hist(x=df['exp_residue'],histtype='step')

plt.xlabel('Residue (Absolute fit error)')
plt.ylabel('Counts')
#
plt.title('Residue distribution for $^{60}Co$, channel 6')
plt.savefig('plots/hyp_testing/hist_residues_co6.png')
plt.close()
plt.figure()
plt.scatter(x=df['lin'].values,y=df['residue'].values,c='r')
plt.axhline(y=0,c='k')
plt.xlim(min(df['lin'])-0.001,max(df['lin'])+0.001)
plt.ylim(min(df['residue'])-0.001,max(df['residue'])+0.001)
#
plt.title('Versus Fits for the log activity of $^{60}$Co, channel 6')
plt.xlabel('Fitted Value')
plt.ylabel('Residual')
plt.savefig('plots/hyp_testing/residue_co6.png')
plt.close()
ks_stat,p = kstest_normal(df['residue'],dist='norm',pvalmethod='approx')
print(p)
print(durbin_watson(df['residue']))
"""
sys.exit()
"""
### Exponential fits and half-lives
for source in df_ana['source'].unique():
    if source == 'co':
        continue
    df = df_ana[df_ana['source'] == source]
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
    if source == 'cs':
        plt.title('Full activity of $^{137}$Cs photopeaks')
        x,y = (-0.2,179.55)
        plt.text(x,y,'$\lambda$ = '+str('{:.2e}'.format(coeff[1])) +' $\pm$'+str('{:.2e}'.format(perr[1]))+' [days]$^{-1}$'+'\n' + '$T_{1/2}$ = '+str(round(np.log(2)/(coeff[1]*365),2)) + ' $\pm$'+str(round((np.log(2)*perr[1])/(np.power(coeff[1],2)*365),2)) + ' years',bbox=dict(facecolor='red',alpha=0.5))
    else:
        plt.title('Full activity of $^{44}$Ti photopeaks')
        x,y = (-0.2,472.01)
        plt.text(x,y,'$\lambda$ = '+str('{:.2e}'.format(coeff[1])) +' $\pm$'+str('{:.2e}'.format(perr[1]))+' [days]$^{-1}$'+'\n' + '$T_{1/2}$ = '+str(round(np.log(2)/(coeff[1]*365),2)) + ' $\pm$'+str(round((np.log(2)*perr[1])/(np.power(coeff[1],2)*365),2)) + ' years',bbox=dict(facecolor='red',alpha=0.5))   
    plt.savefig('plots/agg_peaks/'+str(source)+'.png')
    plt.close()
"""
"""
### Hypothesis testing framework
for source in df_ana['source'].unique():
    if source == 'co':
        continue
    df = df_ana[df_ana['source'] == source]
    p0 = [df['rate'].max(),0.001]
    coeff, var_matrix = curve_fit(exp_fit, df['time'], df['rate'], p0=p0, maxfev = 50000)
    perr = np.sqrt(np.diag(var_matrix))
    x = df['time'].values
    y = np.log(df['rate'])
    y_err = df['drate']/df['rate']
    
    #
    p0 = [max(y),-1]
    coeff2, var_matrix2 = curve_fit(linear_fit, x, y, p0=p0, maxfev = 50000)
    sl, inter, rv, p, stderr = linregress(x=x,y=y)
    #
    e_fit = exp_fit(x,*coeff)
    slope = coeff[1]
    intercept = np.log(coeff[0])
    coeff = [intercept,slope]
    df['lin'] = np.vectorize(linear_fit)(df['time'].values,*coeff)
    df['residue'] = df['lin'] - y
    df['exp_residue'] = e_fit-df['rate']
    r_squared = 1 - (np.sum(np.power(df['residue'],2)))/(np.sum(np.power(y-np.mean(y),2)))
    fit = linear_fit(x, *coeff)
    plt.figure()
    plt.errorbar(x,y,yerr=y_err,fmt='.')
    plt.plot(x,fit)
    if source == 'ti':
        plt.title('The logarithmic activity of the $^{44}$Ti source')
        plt.text(0,6.157,'$R^2 = $'+str(round(r_squared,4)),bbox=dict(facecolor='red',alpha=0.5))
    else:
        plt.title('The logarithmic activity of the $^{137}$Cs source')
        plt.text(0,min(y),'$R^2 = $'+str(round(r_squared,4)),bbox=dict(facecolor='red',alpha=0.5))
    plt.ylabel('log(A(t))')
    plt.xlabel('Time (days)')
    plt.savefig('plots/hyp_testing/'+str(source)+'_linfit.png')
    plt.close()
    plt.figure()
    plt.hist(x=df['residue'],histtype='step')
    plt.xlabel('Residue (Absolute fit error)')
    plt.ylabel('Counts')
    if source == 'ti':
        plt.title('Residue distribution for $^{44}$Ti')
    else:
        plt.title('Residue distribution for $^{137}$Cs')
    plt.savefig('plots/hyp_testing/hist_residues_'+str(source)+'.png')
    plt.close()
    plt.figure()
    plt.scatter(x=df['lin'].values,y=df['residue'].values,c='r')
    plt.axhline(y=0,c='k')
    plt.xlim(min(df['lin'])-0.0001,max(df['lin'])+0.0001)
    plt.ylim(min(df['residue'])-0.0001,max(df['residue'])+0.0001)
    if source == 'ti':
        plt.title('Versus Fits for the log activity of $^{44}$Ti')
    else:
        plt.title('Versus Fits for the log activity of $^{137}$Cs')
    plt.xlabel('Fitted Value')
    plt.ylabel('Residual')
    plt.savefig('plots/hyp_testing/residue_'+str(source)+'.png')
    plt.close()
    ks_stat,p = kstest_normal(df['residue'],dist='norm',pvalmethod='approx')
    print(source,(p))

"""
"""
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
        plt.xlabel('Frequency [Days]')
        plt.ylabel('Normalized Spectral Power')
        plt.xlim([0.1,15])
        plt.savefig('plots/lomb_scargle/ls_'+str(peak)+'_'+str(chn)+'.png')
        plt.close()
        """
