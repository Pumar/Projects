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

Cs_peaks = (661.7)
Co_peaks = (1173.2, 1332.5)

hv_file_list = ['/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181210_1857',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181210_1914',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181210_1930',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181210_1946',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181211_1628',
                '/home/caio/Documents/processing_data/HV_Scan_2018_12/Processed/mx_b_20181211_1644']

sig = r'$\sigma$'
mu = r'$\mu$'
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
    df = df.groupby(pd.cut(df['integral'],bins=700)).agg('count').rename(columns={'integral' : 'Count'}).reset_index()
    df = df.rename(columns={'integral' : 'energy'})
    df['energy'] = df['energy'].astype('str')
    df['energy'] = df['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df = df[['energy','Count']]
    #Poisson error for each "bin"
    df['y_err'] = np.sqrt(df['Count'])
    return df

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def exp_gauss(x, *p):
    A, mu, sigma, a, b = p
    return (A*np.exp(-(x-mu)**2/(2.*sigma**2)) +a*np.exp(-b*x))

def exp_double_gauss(x, *p):
    A1, mu1, sigma1, A2, mu2, sigma2, a, b = p
    return ((A1*np.exp(-(x-mu1)**2/(2.*sigma1**2))) + (A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))) + (a*np.exp(-b*x)))

def exp_fit(x, *p):
    a, b = p
    return a*np.exp(-b*x)

def peak_width(df,channel):
    ###########################################################################
    #
    ### Input:
    # df - dataframe with data to be fit (pd.DataFrame)
    # channel - channel number corresponding to the detector (int)
    ##
    ### Obs:
    ### f = A/B, used approx from https://www.sagepub.com/sites/default/files/upm-binaries/6427_Chapter_4__Lee_(Analyzing)_I_PDF_6.pdf
    ###check wiki page for error propagation
    #
    ###########################################################################
    df_spec = specdf(df[df.error==0],channel)
    if channel == 0 or channel == 1:
        df_spec = df_spec[df_spec['energy']>470]
        df_spec = df_spec[df_spec['energy']<(Cs_peaks+200)]
        ### p0 has the estimates for fitting coefficients (Height of peak{A}, Mean{mu}, Standard Deviation{sigma})
        p0 = [df_spec['Count'].max().astype('float'), Cs_peaks, 20., 2000, 0.001]
        bounds = ([-np.inf,-np.inf,-np.inf,-np.inf,0],[np.inf,np.inf,np.inf,np.inf,np.inf])
        coeff, var_matrix = curve_fit(exp_gauss, df_spec['energy'], df_spec['Count'], p0=p0, bounds=bounds, maxfev = 5000)
        perr = np.sqrt(np.diag(var_matrix))
        dm = coeff[2]/coeff[1]
        return coeff, perr, dm, np.absolute(dm)*np.sqrt(np.power(perr[2]/coeff[2],2)+np.power(perr[1]/coeff[1],2))
    else:
        df_spec = df_spec[df_spec['energy']>(1000)]
        df_spec = df_spec[df_spec['energy']<(Co_peaks[1]+350)]  
        ### p0 has the estimates for fitting coefficients (Height of peak{A}, Mean{mu}, Standard Deviation{sigma})x2 + (exp C1,C2) 
        p0 = [df_spec['Count'].max().astype('float'), Co_peaks[0], 20., df_spec['Count'].max().astype('float')*0.8, Co_peaks[1], 20., 1000, 0.00001]
        bounds = ([-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,0],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
        coeff, var_matrix = curve_fit(exp_double_gauss, df_spec['energy'], df_spec['Count'], p0=p0, maxfev = 5000)
        perr = np.sqrt(np.diag(var_matrix))
        dm = coeff[2]/coeff[1]
        ### f = A/B, used approx from https://www.sagepub.com/sites/default/files/upm-binaries/6427_Chapter_4__Lee_(Analyzing)_I_PDF_6.pdf
        #check wiki page for error propagation
        return coeff, perr, dm, np.absolute(dm)*np.sqrt(np.power(perr[2]/coeff[2],2)+np.power(perr[1]/coeff[1],2))

def plot_calib(hv,de,err,channel):
    title = 'deltaE_X_highvoltage_channel'+str(channel)+'.png'
    d = {'hv' : hv, 'de' : de, 'err' : err}
    df = pd.DataFrame(data=d)
    plt.figure()
    df.plot(x='hv',y='de',yerr = 'err', kind='scatter',title = 'Normalized Energy Peak Width (standard deviation/mean) X High Voltage ')
    plt.xlabel('Voltage (V)')
    plt.ylabel(r'$\sigma$'+'/'+'$\mu$')
    plt.savefig(plotoutDir + '/hv_scan/'+title)
    plt.close()

def plot_spectra(df,channel,hv,pf, pf_err):
    title = 'hv'+str(hv)
    df_spec = specdf(df[df.error==0],channel)
    Emax = int((df_spec[df_spec['Count']==df_spec['Count'].max()]['energy']).astype('int'))
    Cmax = df_spec['Count'].max()
    if channel == 0 or channel == 1:
        e_range = pd.Series([i for i in range(470,int(Cs_peaks)+200)])
        fit = exp_gauss(e_range,*pf)
        df_spec = df_spec[df_spec['energy']<Cs_peaks*2]
        mult_val = 1.2
    else:
        e_range = pd.Series([i for i in range(1000,int(Co_peaks[1])+350)])
        fit = exp_double_gauss(e_range,*pf)
        df_spec = df_spec[df_spec['energy']<sum(Co_peaks)*1.1]
        mult_val = 1.4
    Emax = int((df_spec[df_spec['Count']==df_spec['Count'].max()]['energy']).astype('int'))
    Cmax = df_spec['Count'].max()
    plt.figure()
    df_spec.plot(x='energy',y='Count',yerr= 'y_err', kind='scatter',s=4,color='k', title = channel_label[channel] + ' - Detector '+str(channel)+' with Voltage = '+str(hv) + ' V')
    plt.plot(e_range, fit,color='r')
    plt.text(Emax*mult_val, Cmax/2,'Main Energy Peak'+'\n'+r'$\mu$ = '+str(int(pf[1]))+r' $\pm$'+str(pf_err[1])[:3]+' keV'+ '\n' + r'$\sigma$ = '+str(coeff[2])[:5]+r' $\pm$'+str(pf_err[2])[:4] +' keV')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Event Count')
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
        coeff, coeff_err, de_mu,de_mu_err = peak_width(dataframe,chn)
        err_l.append(de_mu_err)
        de_l.append(de_mu)
        hv_s = hv0[chn]-hv_file_list.index(mod_dir)*25
        plot_spectra(dataframe,chn,hv_s,coeff, coeff_err)
    plot_calib(hv_l,de_l,err_l,chn)
    print('channel'+str(chn)+' plotted')
        
    
