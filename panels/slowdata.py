import os
import sys
import configparser

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import pandas as pd
import numpy as np
from root_pandas import read_root

### Insert path to most recent directory with slow data
_,mostRecentDir = sys.argv

config = configparser.ConfigParser()
config.read('config.ini')

plotoutDir = config['Control']['plotoutdir']

def getconfigrange(section, value):
    return eval(config[section][value])

alarmranges = [getconfigrange('AlarmRanges', 'temperature'),
               getconfigrange('AlarmRanges', 'pressure'),
               getconfigrange('AlarmRanges', 'humidity'),
               getconfigrange('AlarmRanges', 'magfield')]


dataIndices = ['temp', 'pres', 'humid', 'btot']
hvs = ['hv'+str(i) for i in range(0,8)]
axisLabels = ["Temperature [Deg]", "Pressure [Pa]", "Humidity [%]", "Magnetic field [Xgauss]"]
plotTitles = ["Temperature", "Pressure", "Humidity", "Magnetic field"]
fileNameStems = ['temperature', 'pressure', 'humidity', 'magfield']



###Get slow data
for rawfile in os.listdir(mostRecentDir):
	if rawfile.endswith('.sroot'):
		sroot = rawfile



slowdata = read_root(mostRecentDir+'/'+sroot)
slowdata['stime'] = slowdata['stime'] - 2208988800 # Convert from Mac to UNIX time
slowdata['stime'] = pd.to_datetime(slowdata['stime'], unit = 's')
slowdata = slowdata.rename(columns = {'stime' : 'time'})



#Outlier filtering (when arduino and labview are out of sync, wild values for slow control are recorded)
slowdata = slowdata[(slowdata['temp'] > slowdata['temp'].median() - 5) & (slowdata['temp'] < slowdata['temp'].median() + 5)]
slowdata = slowdata[(slowdata['pres'] > slowdata['pres'].median() - 2000) & (slowdata['pres'] < slowdata['pres'].median() + 2000)]
slowdata = slowdata[(slowdata['btot'] > (slowdata['btot'].median() - 50)) & (slowdata.btot < (slowdata['btot'].median() + 50))]
#slowdata = slowdata[(slowdata['humid'] > slowdata['humid'].median() - 15) & (slowdata['humid'] < slowdata['humid'].median() + 15)]


for hv in hvs:
    plt.figure()
    slowdata.plot(x='time',y=hv,c='k')
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.ylim(slowdata[hv].median() - 1.5, slowdata[hv].median() + 1.5)
    plt.title('High Voltage '+hv[-1])
    plt.savefig(plotoutDir+'/'+'high_voltage'+hv[-1]+".png")
    plt.close()

for i in range(4):
    plt.figure()
    slowdata.plot(x='time',y=dataIndices[i],c='k')
    plt.axhline(y=alarmranges[i][0],color='r')
    plt.axhline(y=alarmranges[i][1],color='r')
    plt.fill_between(slowdata['time'].values,alarmranges[i][0],alarmranges[i][1],facecolor='green',alpha=0.2)
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.ylim(alarmranges[i][0] - (alarmranges[i][1]-alarmranges[i][0]),alarmranges[i][1] + (alarmranges[i][1]-alarmranges[i][0]) )
    plt.ylabel(axisLabels[i])
    plt.title(plotTitles[i])
    plt.savefig(plotoutDir+'/'+fileNameStems[i]+'.png')

