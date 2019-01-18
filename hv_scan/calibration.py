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


plotoutDir = "/home/caio/Documents/Plots/"

channel_label = ['Cs-137','Cs-137','Co-60','Co-60']

hv_file_list = ['/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181210_1857',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181210_1914',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181210_1930',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181210_1946',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181211_1628',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181211_1644']

"""
### Care to change to corrrect hv0
hv_file_list = ['/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20190115_1737']
hv0 = (645, 635, 875, 725)
###
"""
hv0 = (720, 710, 950, 800)
hv_step = 25


##############################################################################################################################

def append_dfs(root_dir):
    if len(rootf)>1:
	    framesroot = [read_root(root_dir+'/'+rf, columns=['channel','time','error','integral']) for rf in rootf]
	    df = pd.concat(framesroot,axis=0)
    else:
	    df = read_root(root_dir+'/'+rootf[0], columns=['channel','time','error','integral'])
    return df

def specdf(df,channel):
    df = df[(df['channel']==channel)][['time','integral']]
    df = df.groupby(pd.cut(df['integral'],bins=500)).agg('count').rename(columns={'integral' : 'Count'}).reset_index()
    df = df.rename(columns={'integral' : 'energy'})
    df['energy'] = df['energy'].astype('str')
    df['energy'] = df['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df = df[['energy','Count']]
    return df

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def peak_width(df,channel):
    df_spec = specdf(df[df.error==0],channel)
    Emax = float(df_spec[df_spec['Count']==df_spec['Count'].max()]['energy'].astype('float'))
    ### p0 has the estimates for fitting coefficients (Height of peak{A}, Mean{mu}, Standard Deviation{sigma})
    p0 = [df_spec['Count'].max().astype('float'), Emax, 20.]
    bounds = ([-np.inf, p0[1]*0.9, p0[2]*0.4],[np.inf,p0[1]*1.1, p0[2]*2.5])
    coeff, var_matrix = curve_fit(gauss, df_spec['energy'], df_spec['Count'], p0=p0, bounds = bounds)
    print(coeff[1])
    perr = np.sqrt(np.diag(var_matrix))
    dm = coeff[2]/coeff[1]
    ### f = A/B, used approx from https://www.sagepub.com/sites/default/files/upm-binaries/6427_Chapter_4__Lee_(Analyzing)_I_PDF_6.pdf
    #check wiki for error propagation
    return coeff,dm,np.absolute(dm)*np.sqrt(np.power(perr[2]/coeff[2],2)+np.power(perr[1]/coeff[1],2))

def plot_calib(hv,de,err,channel):
    title = 'deltaE_X_highvoltage_channel'+str(channel)+'.png'
    d = {'hv' : hv, 'de' : de, 'err' : err}
    df = pd.DataFrame(data=d)
    plt.figure()
    df.plot(x='hv',y='de',yerr = 'err', kind='scatter',title = ' Unitless Energy Peak Width X High Voltage ')
    plt.savefig(plotoutDir + '/hv_scan/'+title)
    plt.close()

def plot_spectra(df,channel,hv,pf):
    title = 'hv'+str(hv)
    df_spec = specdf(df[df.error==0],channel)
    Emax = int((df_spec[df_spec['Count']==df_spec['Count'].max()]['energy']).astype('int'))
    Cmax = df_spec['Count'].max()
    fit = gauss(df_spec['energy'],*pf)
    plt.figure()
    df_spec.plot(x='energy',y='Count', title = 'detector'+str(channel))
    plt.plot(df_spec['energy'], fit)
    plt.text(Emax*1.8, Cmax/2,'Epeak = '+str(Emax)+'keV')
    plt.savefig(plotoutDir + '/hv_scan/spectra/'+'detector'+str(channel)+'_'+title+'.png')
    plt.close()


for chn in range(0,4):
    hv_l = []
    de_l = []
    err_l = []
    for mod_dir in hv_file_list:
        rootf=[]
        for file in os.listdir(mod_dir):
	        if file.endswith('.root'):
		        rootf.append(file)
        ## Create dataframe for every data acquisition and overwrite previous one
        dataframe = append_dfs(mod_dir)
        hv_l.append(hv0[chn]-hv_file_list.index(mod_dir)*hv_step)
        coeff,de_mu,de_mu_err = peak_width(dataframe,chn)
        err_l.append(de_mu_err)
        de_l.append(de_mu)
        hv_s = hv0[chn]-hv_file_list.index(mod_dir)*25
        plot_spectra(dataframe,chn,hv_s,coeff)
    plot_calib(hv_l,de_l,err_l,chn)
    print('channel'+str(chn)+' plotted')
        
    