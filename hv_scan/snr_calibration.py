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


##############################################################################################################################
#
#
###                                                GLOBAL VARIABLES
#
#
plotoutDir = "/home/caio/Documents/Plots/hv_scan_2019_02"

channel_label = ['Cs-137','Cs-137','Co-60','Co-60']

Cs_peaks = (661.7)
Co_peaks = (1173.2, 1332.5)
Ti_peaks = (511.,1157.02)

hv_file_list = ['/home/caio/Documents/processing_data/hv_scan/HV_Scan_2019_02/Processed/mx_b_20190225_1734',
                '/home/caio/Documents/processing_data/hv_scan/HV_Scan_2019_02/Processed/mx_b_20190225_1745',
                '/home/caio/Documents/processing_data/hv_scan/HV_Scan_2019_02/Processed/mx_b_20190225_1756',
                '/home/caio/Documents/processing_data/hv_scan/HV_Scan_2019_02/Processed/mx_b_20190225_1807',
                '/home/caio/Documents/processing_data/hv_scan/HV_Scan_2019_02/Processed/mx_b_20190225_1818',
                '/home/caio/Documents/processing_data/hv_scan/HV_Scan_2019_02/Processed/mx_b_20190225_1829',
                '/home/caio/Documents/processing_data/hv_scan/HV_Scan_2019_02/Processed/mx_b_20190225_1840',
                '/home/caio/Documents/processing_data/hv_scan/HV_Scan_2019_02/Processed/mx_b_20190225_1851']

hv0 = (760, 730, 660, 680)
hv_step = 20
#
#
#
##############################################################################################################################

def append_dfs(root_dir):
    ###########################################################################
    ### Description:
    #
    ## Takes each indivisual root file that's been processed (Processing - Joran version), creates a Pandas dataframe and appends all of them into a single Dataframe.
    #
    ### Input:
    #
    ## root_dir - File directory (list of strings)
    #
    ### Returns:
    #
    ## df - appended Dataframe (pd.DataFrame)
    #
    ###########################################################################
    if len(rootf)>1:
        framesroot = [read_root(root_dir+'/'+rf, columns=['channel','time','error','integral']) for rf in rootf]
        df = pd.concat(framesroot,axis=0)
    else:
        df = read_root(root_dir+'/'+rootf[0], columns=['channel','time','error','integral'])
    return df

def specdf(df,channel):
    ###########################################################################
    ### Description:
    #
    ## Bins events into evenly sized energy ranges that depend on the overall range
    #
    ### Input:
    #
    ## df - dataframe which will have its data binned to have a spectral format (pd.DataFrame)
    ## channel - channel number corresponding to the detector (int)
    #
    ### Returns:
    #
    ## df - transformed dataframe (pd.DataFrame)
    #
    ###########################################################################
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
    ##########################################################
    ### Gaussian function
    #
    ## Check https://en.wikipedia.org/wiki/Gaussian_function
    #
    ### Input:
    #
    ## A, mu, sigma - fit params
    #
    ##########################################################
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def exp_gauss(x, *p):
    ##########################################################
    ### Gaussian + exponential function
    #
    ### Input:
    #
    ## A, mu, sigma, a, b - fit params
    #
    ##########################################################
    A, mu, sigma, a, b = p
    return (A*np.exp(-(x-mu)**2/(2.*sigma**2)) +a*np.exp(-b*x))

def exp_double_gauss(x, *p):
    ##########################################################
    ### Gaussian + Gaussian + exponential function
    #
    ### Input: 
    #
    ## A1, mu1, sigma1, A2, mu2, sigma2, a, b - fit params
    #
    ##########################################################
    A1, mu1, sigma1, A2, mu2, sigma2, a, b = p
    return ((A1*np.exp(-(x-mu1)**2/(2.*sigma1**2))) + (A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))) + (a*np.exp(-b*x)))

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

def peak_width(df,channel):
    ###########################################################################
    #
    ### Input:
    #
    ## df - dataframe with data to be fit (pd.DataFrame)
    ## channel - channel number corresponding to the detector (int)
    #
    ## Returns:
    #
    ## coeff - coefficients of the function that was fit to the data inside df (list of floats)
    ## perr - error associated with each coefficient fit (list of floats)
    ## dm - optimize parameter, built as sigma/mu (float)
    ## error propagator - calculates and returns the propagation of the error of the function sigma/mu (float)
    #
    ### Obs:
    ## f = A/B, used approx from https://www.sagepub.com/sites/default/files/upm-binaries/6427_Chapter_4__Lee_(Analyzing)_I_PDF_6.pdf
    ## check wiki page for error propagation
    #
    ###########################################################################
    df_spec = specdf(df[df.error==0],channel)
    if channel == 4 or channel == 5:
        df_spec = df_spec[df_spec['energy']>470]
        df_spec = df_spec[df_spec['energy']<(Cs_peaks+200)]
        ### p0 has the estimates for fitting coefficients (Height of peak{A}, Mean{mu}, Standard Deviation{sigma})
        p0 = [df_spec['Count'].max().astype('float'), Cs_peaks, 20., 2000, 0.001]
        bounds = ([-np.inf,-np.inf,-np.inf,-np.inf,0],[np.inf,np.inf,np.inf,np.inf,np.inf])
        coeff, var_matrix = curve_fit(exp_gauss, df_spec['energy'], df_spec['Count'], p0=p0, bounds=bounds, maxfev = 5000)
        perr = np.sqrt(np.diag(var_matrix))
        dm = coeff[2]/coeff[1]
        #
        print('cs_exp = '+str(coeff[-2]))
        #
        ### f = A/B, used approx from https://www.sagepub.com/sites/default/files/upm-binaries/6427_Chapter_4__Lee_(Analyzing)_I_PDF_6.pdf
        #check wiki page for error propagation
        #return coeff, perr, dm, np.absolute(dm)*np.sqrt(np.power(perr[2]/coeff[2],2)+np.power(perr[1]/coeff[1],2))
        return df_spec['Count'].max().astype('float'),coeff, perr, dm, np.absolute(dm)*np.sqrt(np.power(perr[2]/coeff[2],2)+np.power(perr[1]/coeff[1],2))
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
        #return coeff, perr, dm, np.absolute(dm)*np.sqrt(np.power(perr[2]/coeff[2],2)+np.power(perr[1]/coeff[1],2))
        return df_spec['Count'].max().astype('float'),coeff, perr, dm, np.absolute(dm)*np.sqrt(np.power(perr[2]/coeff[2],2)+np.power(perr[1]/coeff[1],2))

def exponential_background(exp_coeff, lin_coeff, min_lim, max_lim):
    return ((lin_coeff/exp_coeff)*(np.exp(-exp_coeff*min_lim) - np.exp(-exp_coeff*max_lim)))

def snr(signal, background):
    return signal/np.sqrt(background)

def snr_test(signal,N):
    return signal/np.sqrt(N)
    
 
def plot_calib(hv,de,err,channel):
    ###########################################################################
    ### Description:
    #
    ## Plots the accuracy metric (sigma/mu) with error bars against the high voltage values, that correspond to each to a different data acquisition file. 
    ## File name ex: mx_b_20181210_1857
    #
    ### Inputs:
    #
    ## hv - high voltage values (list of ints)
    ## de - sigma/mu (list of floats)
    ## err - error propagated through the function sigma/mu (list of floats)
    ## channel - channel number corresponding to the detector (int)
    #
    ### Returns:
    #
    ## None
    #
    ###########################################################################
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
    ###########################################################################
    ### Description:
    #
    ## Plots the spectra for each channel in a specific hv value and its fit
    ## It's mostly a visual to check goodness of fit.
    #
    ### Inputs:
    #
    ## df - Dataframe that had its data fit, it is utilized to be compared against the plot of the fit function (pd.DataFrame)
    ## hv - high voltage values (list of ints)
    ## channel - channel number corresponding to the detector (int)
    ## pf - coefficients for the fit function (list of float)
    ## pf_err - error associated with each parameter pf[i] (list of float)
    #
    ### Returns:
    #
    ## None
    #
    ###########################################################################
    title = 'hv'+str(hv)
    df_spec = specdf(df[df.error==0],channel)
    Emax = int((df_spec[df_spec['Count']==df_spec['Count'].max()]['energy']).astype('int'))
    Cmax = df_spec['Count'].max()
    if channel == 4 or channel == 5:
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
    df_spec.plot(x='energy',y='Count',yerr= 'y_err', kind='scatter', s=4, color='k', title = channel_label[channel-4] + ' - Detector '+str(channel)+' with Voltage = '+str(hv) + ' V')
    plt.plot(e_range, fit,color='r')
    plt.text(Emax*mult_val, Cmax/2,'Main Energy Peak'+'\n'+r'$\mu$ = '+str(int(pf[1]))+r' $\pm$'+str(pf_err[1])[:3]+' keV'+ '\n' + r'$\sigma$ = '+str(coeff[2])[:5]+r' $\pm$'+str(pf_err[2])[:4] +' keV')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Event Count')
    plt.savefig(plotoutDir + '/hv_scan/spectra/'+'detector'+str(channel)+'_'+title+'.png')
    plt.close()

def plot_snr(hv,snr,channel):
    title = 'snr_hv_channel' + str(channel) + '.png'
    d = {'hv' : hv, 'snr' : snr}
    df = pd.DataFrame(data=d)
    plt.figure()
    df.plot(x='hv',y='snr', kind='scatter',title = 'Signal to Noise Ratio X High Voltage ')
    plt.xlabel('Voltage (V)')
    plt.ylabel('SNR (unitless)')
    plt.savefig(plotoutDir + '/hv_scan/'+title)
    plt.close()




###########################################################################
#
### Main Loop that calls the functions
#
###########################################################################
for chn in range(4,8):
    hv_l = []
    de_l = []
    err_l = []
    snr_l = []
    for mod_dir in hv_file_list:
        rootf=[]
        for file in os.listdir(mod_dir):
            if file.endswith('.root'):
                rootf.append(file)
        ## Create dataframe for every data acquisition (separate file names) and overwrite previous one
        dataframe = append_dfs(mod_dir)
        hv_l.append(hv0[chn-4]-hv_file_list.index(mod_dir)*hv_step)
        n_max,coeff, coeff_err, de_mu,de_mu_err = peak_width(dataframe,chn)
        #snr_l.append(snr(coeff[0],exponential_background(coeff[-1], coeff[-2], coeff[1] - coeff[2], coeff[1] + coeff[2])))
        snr_l.append(snr_test(coeff[0],n_max))
        err_l.append(de_mu_err)
        de_l.append(de_mu)
        hv_s = hv0[chn-4]-hv_file_list.index(mod_dir)*hv_step
        plot_spectra(dataframe,chn,hv_s,coeff, coeff_err)
    plot_calib(hv_l,de_l,err_l,chn)
    plot_snr(hv_l,snr_l,chn)
    print('channel'+str(chn)+' plotted')
        
    
