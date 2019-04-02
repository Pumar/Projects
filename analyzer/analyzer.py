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


### Insert path to most recent directory with fast data

#_,mostRecentDir = sys.argv


simulation_files = ('NULL','MC_ti44_modulation.root', 'MC_cs137_modulation.root', 'MC_co60_modulation.root')


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
        return unpacked_hist, [5*e for e in range(600)]

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
    df = df.groupby(pd.cut(df['integral'],bins=600)).agg('count').rename(columns={'integral' : 'Count'}).reset_index()
    df = df.rename(columns={'integral' : 'energy'})
    df['energy'] = df['energy'].astype('str')
    df['energy'] = df['energy'].str.split(',').str[0].str.split('.').str[0].str[1:].astype('int')
    df = df[['energy','Count']]
    #Poisson error for each "bin"
    df['y_err'] = np.sqrt(df['Count'])
    return df

count,energy = compton_simulation(simulation_files[3])


plt.figure()
plt.scatter(energy,count,c='k',s=1)
plt.savefig('compton.png')
plt.close()

