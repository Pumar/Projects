import pandas as pd
import matplotlib
import numpy as np 
import pylab as pl 


df = pd.read_csv('turnstile_weather_v2.csv')

df_r = np.array(df['ENTRIESn_hourly'][df['rain']==1])
df_nr = np.array(df['ENTRIESn_hourly'][df['rain']==0])
muR = np.mean(df_r)
mu0R = np.mean(df_nr)
"""
print muR
print mu0R

pl.figure(1)
pl.subplot(211)
pl.hist(df_r,60,facecolor='blue')
pl.text(5000,2000,r'$\mu=2029$')
pl.title('Ridership during a rainy day')

pl.subplot(212)
pl.hist(df_nr,60,facecolor='gray')
pl.text(5000,6000,r'$\mu=1846$')
pl.title('Ridership on a clear day')
pl.show()
"""
dfE = np.array(df['ENTRIESn_hourly'][df['ENTRIESn_hourly']!=0])
dfday = np.array(df['day_week'][df['ENTRIESn_hourly']!=0])
pl.plot(dfday,dfE,'ko')
pl.title('Ridership by day of the week')
pl.xlabel('Monday = 0 ... Sunday = 6')
pl.show()