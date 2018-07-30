import os
import sys
import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from root_pandas import read_root

_,mostRecentDir = sys.argv
plotoutDir = '/home/caio/Documents/Plots'

for file in os.listdir(mostRecentDir):
	if file.endswith('.sroot'):
		sroot = file

slowdata = read_root(mostRecentDir+'/'+sroot)
slowdata['stime'] = slowdata['stime'] - 2208988800 # Convert from Mac to UNIX time
slowdata['stime'] = pd.to_datetime(slowdata['stime'], unit = 's')
slowdata = slowdata.rename(columns = {'stime' : 'time'})

dataIndices = ['temp', 'pres', 'humid', 'btot']
hvs = ['hv'+str(i) for i in range(0,8)]
axisLabels = ["Temperature [Deg]", "Pressure [Pa]", "Humidity [%]", "Magnetic field [Xgauss]"]
plotTitles = ["Temperature", "Pressure", "Humidity", "Magnetic field"]
fileNameStems = ['temperature', 'pressure', 'humidity', 'magfield']

#Outlier filtering (when arduino and labview are out of sync, wild values for slow control are recorded)
slowdata = slowdata[(slowdata.temp > 25) & (slowdata.temp < 35)]
slowdata = slowdata[(slowdata.pres > 101000) & (slowdata.pres < 102000)]
slowdata = slowdata[(slowdata.btot > 0) & (slowdata.btot < 1000)]
slowdata = slowdata[(slowdata.humid > 40) & (slowdata.humid < 50)]
for hv in hvs:
    slowdata = slowdata[(slowdata[hv] > slowdata[hv].mean() - 3) & (slowdata[hv] < slowdata[hv].mean() + 3)]


for hv in hvs:
    plot= slowdata.plot(x='time',y=hv)
    fig = plot.get_figure()
    fig.savefig(plotoutDir+'/'+hv+"_output.png")


