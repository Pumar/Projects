import os
import sys
import configparser

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from root_pandas import read_root

### Insert path to most recent directory with slow data
_,mostRecentDir = sys.argv

config = configparser.ConfigParser()
config.read('config.ini')

plotoutDir = config['SlowControl']['plotoutdir']

def getconfigrange(section, value):
    return eval(config[section][value])

pressrangealarm = getconfigrange('AlarmRanges', 'pressure')
temprangealarm = getconfigrange('AlarmRanges', 'temperature')
humidrangealarm = getconfigrange('AlarmRanges', 'humidity')
magrangealarm = getconfigrange('AlarmRanges', 'magfield')

"""
not needed
pressrange = getconfigrange('PlotRanges', 'pressure')
temprange = getconfigrange('PlotRanges','temperature')
humidrange = getconfigrange('PlotRanges','humidity')
magrange = getconfigrange('PlotRanges','magfield')
"""
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


"""
for vl in slowdata.columns:
	print vl+',  std:  '+str(slowdata[vl].std())
"""

#Outlier filtering (when arduino and labview are out of sync, wild values for slow control are recorded)
slowdata = slowdata[(slowdata.temp > 25) & (slowdata.temp < 35)]
slowdata = slowdata[(slowdata.pres > 100000) & (slowdata.pres < 103000)]
slowdata = slowdata[(slowdata.btot > 0) & (slowdata.btot < 1500)]
slowdata = slowdata[(slowdata.humid > 42) & (slowdata.humid < 50)]
"""
for hv in hvs:
    slowdata = slowdata[(slowdata[hv] > slowdata[hv].mean() - 5) & (slowdata[hv] < slowdata[hv].mean() + 5)]
"""

"""
for hv in hvs:
    plot= slowdata.plot(x='time',y=hv)
    fig = plot.get_figure()
    fig.savefig(plotoutDir+'/'+hv+"_output.png")
"""

for di,fn in zip(dataIndices,fileNameStems):
	plot = slowdata.plot(x='time',y=di)
	ax = plt.gca()
	plt.ylim(getconfigrange('PlotRanges',fn))
	fig = plot.get_figure()
	fig.savefig(plotoutDir+'/'+di+"_output.png")
