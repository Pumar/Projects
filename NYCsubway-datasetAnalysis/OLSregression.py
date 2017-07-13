import pandas as pd
import numpy as np 
import statsmodels.api as sm 
import matplotlib
import pylab as pl

df = pd.read_csv('turnstile_weather_v2.csv')


df1 = df.drop('ENTRIESn_hourly',axis=1)
df1 = df1.drop('weather_lat',axis=1)
df1 = df1.drop('weather_lon',axis=1)
#df1 = df1.drop('latitude',axis=1)
#df1 = df1.drop('longitude',axis=1)
df1 = df1.drop('ENTRIESn',axis=1)
df1 = df1.drop('EXITSn',axis=1)
df1 = df1.select_dtypes(exclude=[object])

def lin_reg(dataframe, var_array):

	X=sm.add_constant(dataframe)
	model = sm.OLS(var_array,X)
	results = model.fit()
	return results


def compute_r_squared(data, predictions):
    r_squared = 1 - np.sum(np.square(data-predictions))/np.sum(np.square(data - np.mean(data)))
    return r_squared


results = lin_reg(df1,df['ENTRIESn_hourly'])

thetas = results.params[1:]
residuals=[]
for index,row in df1.iterrows():
	residuals.append(np.dot(thetas,np.array(row)) + results.params[0])

print 'The R^2 value is: '+str(compute_r_squared(df['ENTRIESn_hourly'],residuals))
print thetas