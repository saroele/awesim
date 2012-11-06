# -*- coding: utf-8 -*-

"""
Utility functions for awesim

Created 20120911 by RDC
"""
from __future__ import division
import pandas as pd
import numpy as np
from scipy.integrate import cumtrapz
from scipy.stats import spearmanr
import pdb
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from datetime import datetime, timedelta

def make_datetimeindex(array_in_seconds, year):
    """
    Create a pandas DateIndex from a time vector in seconds and the year.
    """
    
    start = pd.datetime(year, 1, 1)
    datetimes = [start + pd.datetools.timedelta(t/86400.) for t in array_in_seconds]
    
    return pd.DatetimeIndex(datetimes)


def aggregate_by_time(signal, time, period=86400, interval=900, label='left', year=2012):
    """
    Function to calculate the aggregated average of a timeseries by 
    period (typical a day) in bins of interval seconds (default = 900s).
    
    label = 'left' or 'right'.  'Left' means that the label i contains data from 
    i till i+1, 'right' means that label i contains data from i-1 till i.    
    
    Returns an array with period/interval values, one for each interval
    of the period. 
    
    A few limitations of the method:
        - the period has to be a multiple of the interval
            
    This function can be used in the post-processing too.
    """
    
 
    dr = make_datetimeindex(time, year)
    df = pd.DataFrame(data=signal, index=dr, columns=['signal'])
    
    return aggregate_dataframe(df, period, interval, label).values
     
    
def aggregate_dataframe(dataframe, period=86400, interval=3600, label='left'):
    """
    Function to calculate the aggregated average of a timeseries by 
    period (typical a day) in bins of interval seconds (default = 3600s).
    
    label = 'left' or 'right'.  'Left' means that the label i contains data from 
    i till i+1, 'right' means that label i contains data from i-1 till i.    
    
    Returns a new dataframe with period/interval values, one for each interval
    of the period. 
    
    A few limitations of the method:
        - the period has to be a multiple of the interval
            
    Example of usefulness: if the timeseries has 15-minute values for 1 year of
    eg. the electricity consumption of a building.  
    - You want to know how a typical daily profile looks like, by 15 minuts 
      ==> period=86400, interval=900
    - you want to know how a typical weekly profile looks like, by hour:
      ==> period = 7*86400, interval=3600
      
    """
    pdb.set_trace()
    # first, create cumulative integrated signals for every column, put these
    # in a new dataframe called cum

    cum = pd.DataFrame(index=dataframe.index)
    for c in dataframe.columns:
        # we need to remove the empty values for the cumtrapz function to work
        ts = dataframe[c].dropna()
        cum[c] = cumtrapz(ts.values, ts.index.asi8/1e9, initial=0)
  
    # first, resample the dataframe by the given interval   
    # We convert it to milliseconds in order to obtain integer values for most cases
    interval_string = str(int(interval*1000)) + 'L'    
    df_resampled = cum.resample(interval_string, how='last', 
                                      closed=label, label=label)
                                      
    
    df_diff = pd.DataFrame(index=df_resampled.index)    
    for c in df_resampled.columns:
        # diffdata is the average signal during each interval
        reshaped_array = df_resampled[c].values.reshape(len(df_resampled))    
        diffdata = np.zeros(len(reshaped_array))
        diffdata[:-1] = np.diff(reshaped_array)
        df_diff[c] = diffdata
    
    
    
    # now create bins for the groupby() method
    # time in seconds    
    time_s = df_diff.index.asi8/1e9
    time_s -= time_s[0]
    try:
        df_diff['bins'] = np.mod(time_s, period)
    except(KeyError):
        df_diff = pd.DataFrame(df_resampled)
        df_diff['bins'] = np.mod(time_s, period)
        
    
    df_aggr = df_diff.groupby('bins').mean()
    
    # replace the bins by a real datetime index    
    df_aggr.index = df_diff.index[:len(df_aggr)]
    
    return df_aggr


def analyse_cputime_single(cputime, time, var, cumulative=False, interval=900, plot=True):
    """Analyse the relation between cputime and a trajectory.
    
    It is required that the time array contains an entry for every multiple of
    the interval.

    Parameters:
    -----------
    * cputime: array with (trajectory of the) cputime (cumulative)
    * time: simulation time (array in seconds, can contain events)
    * var: variable for which the correlation is to be found
    * cumulative: if True, the variable is a cumulative (integrated) one
    * interval: the desired interval for resampling, integer in seconds
    * plot: if True, a time plot will be made

    Output:
    -------
    Correlation coefficient as float, a screen dump and plot
    """

    if not cumulative:
        var_cum = cumtrapz(var, time)
    else:
        var_cum = var
        
    x = np.arange(time[0], time[-1], interval)
    index = np.searchsorted(time, x)
    cpu_smpl = cputime[index]
    var_smpl = var_cum[index]
    
    cpu_diff = np.diff(cpu_smpl)
    var_diff = np.diff(var_smpl)
    
    corr = spearmanr(var_diff, cpu_diff)[0]
    
    print 'The correlation is %g' % (corr)
    
    if plot:
        start = datetime(2011, 1, 1)
        datetimes = [start + timedelta(t/86400.) for t in x]
        time4plots = date2num(datetimes)
        plt.figure()
        plt.plot_date(time4plots[:-1], cpu_diff/cpu_diff.max(), 'r', label='cpu-time')
        plt.plot_date(time4plots[:-1], var_diff/var_diff.max(), 'b', label='variable')
        plt.legend(loc='best')
        
    return corr
     
     

        
        