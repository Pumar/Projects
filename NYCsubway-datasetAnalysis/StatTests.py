import pandas as pd
from scipy import stats
import matplotlib
import pylab as pl
import statsmodels.api as sm 
import numpy as np 
import datetime
df = pd.read_csv('turnstile_weather_v2.csv')

ent_rain = df['ENTRIESn_hourly'][df['rain']==1.0][df['ENTRIESn_hourly']!=0]
ent_norain = df['ENTRIESn_hourly'][df['rain']==0.0][df['ENTRIESn_hourly']!=0]
#dta = df[['DATEn','TIMEn']][ df['UNIT']]
print (df['UNIT'].str.contains('R464'))
nd = pd.DataFrame()
print nd
"""
U, p = stats.mannwhitneyu(np.array(ent_rain),np.array(ent_norain))
print U,p
W,P = stats.shapiro(np.array(df['ENTRIESn_hourly']))
print P
"""

"""
pl.figure()
pl.hist(np.array(ent_rain),60,alpha=1,facecolor='blue')
pl.hist(np.array(ent_norain),60,alpha=0.5,facecolor='r')
pl.title('Ridership on the NYC subway: rainy and clear days')
pl.show()
pl.close()


kt = np.array(df['ENTRIESn_hourly'][df['ENTRIESn_hourly']!=0])
pl.hist(kt,60, normed=1, facecolor='blue', alpha=0.8)
pl.title('Total ridership, normalised')
pl.show()"""