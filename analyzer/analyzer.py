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



simulation_files = ('NULL','MC_ti44_modulation.root', 'MC_cs137_modulation.root', 'MC_co60_modulation.root')

Cs_peaks = (661.7)
Co_peaks = (1173.2, 1332.5)
Ti_peaks = (511.,1157.02)

def append_dfs(root_dir):
    rootf=[]
    for file in os.listdir(root_dir):
	    if file.endswith('.root'):
		    rootf.append(file)
    if len(rootf)>1:
	    framesroot = [read_root(root_dir+'/'+rf, columns=['channel','time','error','integral']) for rf in rootf]
	    df = pd.concat(framesroot,axis=0)
    else:
	    df = read_root(root_dir+'/'+rootf[0], columns=['channel','time','error','integral'])
    return df

def timedelta(df):
    return (df['time'].max()-df['time'].min())

def time_info(df):
    df = df - 2208988800
    df = pd.to_datetime(df, unit = 's')
    return str(min(df).year) + '-' + str(min(df).month) + '-' +str(min(df).day)

def gauss(x, *p):
    A, mu, sigma = p
    return (A/(sigma*np.sqrt(2*np.pi)))*np.exp(-(x-mu)**2/(2.*sigma**2))

def triple_gauss(x, *p):
    A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3 = p
    return (A1/(sigma1*np.sqrt(2*np.pi)))*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + (A2/(sigma2*np.sqrt(2*np.pi)))*np.exp(-(x-mu2)**2/(2.*sigma2**2)) + (A3/(sigma3*np.sqrt(2*np.pi)))*np.exp(-(x-mu3)**2/(2.*sigma3**2))

def specdf(df,channel):
    df = df[(df['channel']==channel)][['time','integral']]
    td = timedelta(df)
    binnage = [i for i in range( int( int(df['integral'].min()) - ( int(df['integral'].min())%5)),int( int(df['integral'].max()) - ( int(df['integral'].max())%5)), 5)]
    count, bin_edges = np.histogram(df['integral'].values,bins=binnage)
    unpacked_hist = np.true_divide(count.astype('float64'),(td))
    return unpacked_hist,bin_edges[:-1]

class HistogramFile(object):

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self.file = ROOT.TFile.Open(self.filename, 'read')
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.file.Close()

    def get_histogram(self, name):
        """Return the histogram identified by name from the file.
        """
        # The TFile::Get() method returns a pointer to an object stored in a ROOT file.
        hist = self.file.Get(name)
        if hist:
            return hist
        else:
            raise RuntimeError('Unable to retrieve histogram named {0} from {1}'.format(name, self.filename))
    
def compton_simulation(filename):
    with HistogramFile(filename) as f:
        hist = f.get_histogram('h2;1')
        unpacked_hist = []
        for i in range(0,hist.GetNbinsX()):
             unpacked_hist.append(hist.GetBinContent(i))
        return np.array(unpacked_hist), np.array([5*e for e in range(600)])


def subtract_compton(df,channel):
    sim_count, sim_energy = compton_simulation(simulation_files[channel//2])
    fast_count,fast_energy = specdf(df,channel)
    bg_count, bg_energy = specdf(df,0)
    ### Range fixing for simulation data
    sim_energy = np.array(filter(lambda x: x >= min(fast_energy),sim_energy))
    if len(sim_count) != len(sim_energy):
        sim_count = sim_count[len(sim_count)-len(sim_energy):]
    sim_energy = np.array(filter(lambda x: x <= max(fast_energy),sim_energy))
    if len(sim_count) != len(sim_energy):
        sim_count = sim_count[:(len(sim_energy)-len(sim_count))]
    ##
    ### Range fixing for background data
    if len(bg_energy) >= len(fast_energy):
        bg_energy = np.array(filter(lambda x: x >= min(fast_energy),bg_energy))
        if len(bg_count) != len(bg_energy):
            bg_count = bg_count[len(bg_count)-len(bg_energy):]
        bg_energy = np.array(filter(lambda x: x <= max(fast_energy),bg_energy))
        if len(bg_count) != len(bg_energy):
            bg_count = bg_count[:(len(bg_energy)-len(bg_count))]
    else:
        #print(min(bg_energy),max(bg_energy),min(fast_energy),max(fast_energy))
        if bg_energy[-1] != fast_energy[-1]:
            bg_count = np.append(bg_count,np.multiply([i for i in range(max(bg_energy),max(fast_energy),5)],0))
            bg_energy = np.append(bg_energy,[i for i in range(max(bg_energy),max(fast_energy),5)])
        if bg_energy[0] != fast_energy[0]:    
            bg_count = np.append(np.multiply([i for i in range(min(fast_energy),min(bg_energy),5)],0),bg_count)
            bg_energy = np.append([i for i in range(min(fast_energy),min(bg_energy),5)],bg_energy)
    ##
    ### Fitting compton simulation to data close to the photopeaks
    if channel == 2 or channel == 3:
        ratio = max(fast_count[int(np.where(fast_energy ==570)[0]):int(np.where(fast_energy == 800)[0])])/max(sim_count[int(np.where(sim_energy == 570)[0]):int(np.where(sim_energy == 800)[0])])
    elif channel == 4 or channel == 5:
        ratio = max(fast_count[int(np.where(fast_energy ==150)[0]):int(np.where(fast_energy == 250)[0])])/max(sim_count[int(np.where(sim_energy == 150)[0]):int(np.where(sim_energy == 250)[0])])
    else:
        ratio = max(fast_count[int(np.where(fast_energy ==850)[0]):int(np.where(fast_energy == 1000)[0])])/max(sim_count[int(np.where(sim_energy == 850)[0]):int(np.where(sim_energy == 1000)[0])])
    sim_count = sim_count*ratio
    ###
    assert len(bg_count) == len(fast_count) , 'Standardization of energy bins and ranges failed.'
    #
    ### Activity from background + compton scattering is removed
    fast_count = np.subtract(fast_count,bg_count)
    return np.subtract(fast_count,sim_count),fast_energy,sim_count,np.add(fast_count,bg_count),bg_count



def fit_routine(counts,energy_range,channel):
    if channel == 4 or channel == 5:
        p0 = [max(counts)*20*np.sqrt(2*np.pi),Cs_peaks, 20.]
        coeff, var_matrix = curve_fit(gauss, energy_range, counts, p0=p0,absolute_sigma=False, maxfev = 5000)
        perr = np.sqrt(np.diag(var_matrix))
        area = coeff[0]/5.
        area_err = perr[0]/5.
        return np.array(area), np.array(area_err),coeff
    elif channel < 4:
        p0 = [max(counts)*20*np.sqrt(2*np.pi),Ti_peaks[0], 20.,max(counts)*2*np.sqrt(2*np.pi),Ti_peaks[1], 20.,max(counts)*2*np.sqrt(2*np.pi),Ti_peaks[0]+Ti_peaks[1], 20.]
        coeff, var_matrix = curve_fit(triple_gauss, energy_range, counts, p0=p0,absolute_sigma=False, maxfev = 5000)
        perr = np.sqrt(np.diag(var_matrix))
        area = np.true_divide(np.array([coeff[i] for i in range(0,7,3)]),5)
        area_err = np.true_divide(np.array([perr[i] for i in range(0,7,3)]),5)
        return area,area_err,coeff
    else:
        p0 = [max(counts)*20*np.sqrt(2*np.pi),Co_peaks[0], 20.,max(counts)*0.8*20*np.sqrt(2*np.pi),Co_peaks[1], 30.,max(counts)*2*np.sqrt(2*np.pi),Co_peaks[0]+Co_peaks[1], 40.]
        coeff, var_matrix = curve_fit(triple_gauss, energy_range, counts, p0=p0,absolute_sigma=False, maxfev = 5000)
        perr = np.sqrt(np.diag(var_matrix))
        area = np.true_divide(np.array([coeff[i] for i in range(0,7,3)]),5)
        area_err = np.true_divide(np.array([perr[i] for i in range(0,7,3)]),5)
        return area,area_err,coeff

def error_propagation(A,B,sigmaA,sigmaB,sigmaAB):
    return np.abs(A*B)*np.sqrt(np.power(sigmaA/A,2) + np.power(sigmaB/B,2) + 2*sigmaAB/(A*B))

def gaussian_integral(N,sigma):
    return N*sigma*np.sqrt(2*np.pi)


def plot_fits(energy_range,counts,pf,channel):
    plt.figure()
    if channel == 4 or channel == 5:
        fit = gauss(energy_range,*pf)
    else:
        fit = triple_gauss(energy_range,*pf)
    plt.plot(energy_range,fit,c='r')
    plt.scatter(energy_range,counts,c='k',s=1)
    plt.savefig('fit.png')
    plt.close()
    

"""
dataframe = append_dfs('/media/sf_Modulation/Modulation_data/Processed/mx_b_20190401_0305/')
dataframe = dataframe[(dataframe['integral'] >= 0) & (dataframe['integral'] <= 3000) & (dataframe['error'] == 0)][['time','integral','channel']]

chn = 2
final, energy, simulation_count, count, bg_count = subtract_compton(dataframe,chn)
final2,energy2, _, _, _ = subtract_compton(dataframe,chn+1)
final = np.where(final < 0, 0, final)
final2 = np.where(final2 < 0, 0, final2)
print(sum(final),sum(final2))
rate,drate,coeff = fit_routine(final,energy,chn)
rate2,drate2,coeff2 = fit_routine(final2,energy2,chn)
rate = np.add(rate,rate2)
drate = np.sqrt(np.power(drate,2) + np.power(drate2,2))
print('Fit Coefficients 1: '+str(coeff))
print('Channel ' +str(chn)+'. Rates: ' + str(rate) + ' Rates Error: '+str(drate))
plt.figure()
plt.scatter(energy,simulation_count,c='k',s=1)
plt.scatter(energy,count,c='r',s=1)
plt.scatter(energy,final,c='green',s=1)
plt.scatter(energy,bg_count,c='blue',s=1)
plt.axhline(y=0,c='k')
if chn < 4:
    plt.xlim(250,1900)
elif chn < 6:
    plt.xlim(400,800)
else:
    plt.xlim(800, 2800)
plt.savefig(simulation_files[2//2].split('_')[1]+'_compton_bg.png')
plt.close()

sys.exit()
"""

ti_lifetime = []
ti_lifetime_error = []
cs_lifetime = []
cs_lifetime_error = []
co_lifetime = []
co_lifetime_error = []
time = []
cols = ['time','ti1','ti1_err','ti2','ti2_err','ti3','ti3_err','cs','cs_err','co1','co1_err','co2','co2_err','co3','co3_err']
with open('proc_dirs.txt') as f:
    for lines in f:
        lines = lines.rstrip()
        dataframe = append_dfs('/media/sf_Modulation/Modulation_data/Processed/'+lines)
        dataframe = dataframe[(dataframe['integral'] >= 0) & (dataframe['integral'] <= 3000) & (dataframe['error'] == 0)][['time','integral','channel']]
        time.append(min(dataframe['time']))
        for chn in range(2,7,2):
            final, energy, _, _, _ = subtract_compton(dataframe,chn)
            final2,energy2, _, _, _ = subtract_compton(dataframe,chn+1)
            final = np.where(final < 0, 0, final)
            final2 = np.where(final2 < 0, 0, final2)
            #print(sum(final),sum(final2))
            rate,drate,coeff = fit_routine(final,energy,chn)
            rate2,drate2,coeff2 = fit_routine(final2,energy2,chn)
            rate = np.add(rate,rate2)
            drate = np.sqrt(np.power(drate,2) + np.power(drate2,2))
            #print('Fit Coefficients 1: '+str(coeff))
            #print('Channel ' +str(chn)+'. Rates: ' + str(rate) + ' Rates Error: '+str(drate))
            if chn == 2:
                ti_lifetime.append(rate.tolist())
                ti_lifetime_error.append(drate.tolist())
            elif chn == 4:
                cs_lifetime.append(rate)
                cs_lifetime_error.append(drate)
            else:
                co_lifetime.append(rate.tolist())
                co_lifetime_error.append(drate.tolist())


df_rates = pd.DataFrame(list(zip(time,[ti_lifetime[i][0] for i in range(len(ti_lifetime))], [ti_lifetime_error[i][0] for i in range(len(ti_lifetime_error))] , [ti_lifetime[i][1] for i in range(len(ti_lifetime))], [ti_lifetime_error[i][1] for i in range(len(ti_lifetime_error))] , [ti_lifetime[i][2] for i in range(len(ti_lifetime))], [ti_lifetime_error[i][2] for i in range(len(ti_lifetime_error))] ,cs_lifetime, cs_lifetime_error, [co_lifetime[i][0] for i in range(len(co_lifetime))], [co_lifetime_error[i][0] for i in range(len(co_lifetime_error))] , [co_lifetime[i][1] for i in range(len(co_lifetime))], [co_lifetime_error[i][1] for i in range(len(co_lifetime_error))] , [co_lifetime[i][2] for i in range(len(co_lifetime))], [co_lifetime_error[i][2] for i in range(len(co_lifetime_error))])), columns = cols)

df_rates.to_csv('rates.csv', index=False)

plt.figure()
plt.errorbar(x=time,y=cs_lifetime,yerr=cs_lifetime_error,fmt='o')
plt.savefig('lifetime.png')
plt.close()

