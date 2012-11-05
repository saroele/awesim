# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 21:36:52 2012

@author: RDC
"""

import numpy as np
#import os
#import scipy.io
#import re
import copy
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
#import cPickle as pickle
#import bisect
#import tables as tbl
from datetime import datetime, timedelta
from .utilities import make_datetimeindex
import pandas
import pdb
import pickle


class Result(object):
    """
    Class containing the result for one single variable but for different 
    simulations.  An instance from this class is returned from Simdex.get()
    
    This class also contains the plot functionality and methods to apply 
    basic operations and functions to it. 
    """
    
    def __init__(self, values, time=None, identifiers=None, **kwargs):
        """
        Instantiate a Result object. 
        
        Variables
        ---------
        
        values = dictionary {sid:values}
        time (optional) = dictionary {sid:time}
        identifiers (optional) = dictionary {sid:identifier}
        **kwargs are converted into attributes.  This is useful to pass eg. the 
        year for the data (use for example year=2010)
        
        """
        
        self.val = values
        if time is not None:
            self.time = time
            self.time4plots = {}
        if identifiers is not None:
            self.identifiers = identifiers
        self.simulations = sorted(self.val.keys())
        
        for k,v in kwargs.items():
            setattr(self, k, v)

    def save(self, filename):
        """
        save(filename)
        
        Save the Simdex object by pickling it with cPickle
        
        To unpickle (= load) use the following command:
            objectname = pickle.load(open(filename,'rb'))
            # 'rb' stands for 'read, binary'
            
        """
        
        f = file(filename,'wb') 
        # wb stands for 'write, binary'
        pickle.dump(self, f)
        f.close()
        
        return filename + ' created'
            
    def values(self):
        """
        Return a list with as elements, the values of the variable in the order 
        of the sid's
        
        It does not seem a good idea to return an array by default, cause the 
        variables for different SID's can have different lengths. 
        Exception: when the length of each of the variables is 1, a reshaped 
        array is returned.
        
        If a value in a single length array is None, it is replaced by NaN in the
        returned array.
        """
        
        result = [self.val[sid] for sid in self.simulations]
        
        lengths=[]
        for i,x in enumerate(result):
            try:
                l = len(x)
            except TypeError:
                # x is a value, it is a numpy.float64, so length=1
                if x is None:
                    result[i]=np.NaN
                l = 1
            lengths.append(l)
                    
        ls = np.array(lengths)
        
        if np.all(ls==1):
            return np.array(result).reshape(len(ls))
        else:
            return result
        
    def trapz(self):
        """
        Integrate the values(time) using the composite trapezoidal rule
        Returns an array with the integrated values, in sorted order
        """
        
        if not hasattr(self, 'time'):
            raise AttributeError("This Result object has no attribute 'time'")
        
        result = []
        for sid in self.simulations:
            result.append(np.trapz(self.val[sid], x=self.time[sid]))
        
        return np.array(result)
        

    def aggregate(self, period=86400, interval=3600):
        
        def aggregate_by_time(signal, time, period=86400, interval=900, label='left'):
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
            
            def make_datetimeindex(array_in_seconds, year):
                """
                Create a pandas DateIndex from a time vector in seconds and the year.
                """
                
                start = pandas.datetime(year, 1, 1)
                datetimes = [start + pandas.datetools.timedelta(t/86400.) for t in array_in_seconds]
                
                return pandas.DatetimeIndex(datetimes)
            
            interval_string = str(interval) + 'S'    
            dr = make_datetimeindex(time, 2012)
            df = pandas.DataFrame(data=signal, index=dr, columns=['signal'])
            df15min = df.resample(interval_string, closed=label, label=label)
            
            # now create bins for the groupby() method
            time_s = df15min.index.asi8/1e9
            time_s -= time_s[0]
            df15min['bins'] = np.mod(time_s, period)
            
            df_aggr = df15min.groupby(['bins']).mean()
            
            return df_aggr       

        result = {}        
        for sid in self.simulations:
            result[sid] = aggregate_by_time(self.val[sid], self.time[sid], 
                                            period=period, interval=interval)
            
        return result

    def smooth(self, interval=300):
        """
        Calculate the running average of a timeseries
        
        Parameters
        ----------
        interval: interval for the running average, in seconds (default = 300s)
        
        Returns
        -------
        
        returns a result object with smoothened values and adapted time
        """
        
        def smooth_by_time(signal, time, interval=300, label='left'):
            """
            Function to calculate the running average of a timeseries 
            in bins of interval seconds (default = 300s).
            
            """
            #pdb.set_trace()
            ratio = interval/(time[1]-time[0])
            
            time = np.arange(0, time[-1]+interval, interval)
            
            data = np.zeros(len(time))
            data[0] = signal[0]
            for i in range(1,len(data)):
                data[i]=np.mean(signal[((i-1)*ratio):(i*ratio)])
                
            return data, time       
    
        value = {}
        time = {}
        for sid in self.simulations:
            value[sid], time[sid] = smooth_by_time(signal=self.val[sid], time=self.time[sid], interval=interval)

        result = Result(values=value, time=time)
        
        return result


    def to_dataframe(self):
        """
        Return a pandas dataframe from this result
        """
        pdb.set_trace()
        # check existence of attributes
        if not hasattr(self, 'year'):
            print 'We suppose the data is for 2011'
            self.year=2011

        for i, sid in enumerate(sorted(self.val.keys())):
            # create a df from this single 'column'
            index = make_datetimeindex(self.time[sid], self.year)
            df_right = pandas.DataFrame(data=self.val[sid], index=index, columns=[sid])

            if i==0:
                df = copy.deepcopy(df_right)
            else:
                df = df.join(df_right, how='outer', sort=True)

        return df
            



    def plot(self, ylabel=None):
        """
        Creates a matplotlib figure with a simple plot of the timeseries for 
        each of the simulations in self.val
        
        A string can be passed (ylabel) that will be used to label the y-axis
        """
        
        # In order to plot the timeseries nicely with dates, we use plot_date()
        def create_time4plot(sid):
            """Convert time into matplotlib format"""
            start = datetime(self.year, 1, 1)
            datetimes = [start + timedelta(t/86400.) for t in self.time[sid]]
            self.time4plots[sid] = date2num(datetimes)
            

       
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hold = True        

        # we have to make a distinction between plotting timeseries and other results
        try:
            if len(self.time[self.simulations[0]])==len(self.val[self.simulations[0]]): 
                plot_type='plot_date'
            elif len(self.val[self.simulations[0]]) <> 1:
                # most probably an aggregated array. So plot the x values for 
                # each SID as a series 
                plot_type='aggregated'
            else:
                # length == 1, so single_value
                plot_type = 'single_value'
        except:
            try:
                if len(self.val[self.simulations[0]]) <> 1:          
                    # most probably an aggregated array. So plot the x values for 
                    # each SID as a series 
                    plot_type='aggregated'
                else:
                    # length == 1, so single_value
                    plot_type = 'single_value'
            except:
                # I get an exception when trying to get the length of a single
                # intgegrated value
                plot_type = 'single_value'
                
        if plot_type=='plot_date':
            for sid in self.simulations:
                try:
                    label = self.identifiers[sid]
                except KeyError:
                    label=sid
                if not self.time4plots.has_key(sid):
                    create_time4plot(sid)
                    
                ax.plot_date(self.time4plots[sid], self.val[sid], fmt='', ls = '-', 
                             label=label)            
        elif plot_type == 'single_value':
            ax.plot(range(len(self.simulations)), self.values(), 'D')
            ax.set_xticks(range(len(self.simulations)))
            ticklabels = [self.identifiers[sid] for sid in self.simulations]
            ax.set_xticklabels(ticklabels)
        else:
            # aggregated, plot lines with markers            
            for sid in self.simulations:
                try:
                    label = self.identifiers[sid]
                except KeyError:
                    label=sid
                
                ax.plot(self.val[sid], 'o-', label=label)            
                             
        leg = ax.legend(loc='best')
        lines = ax.get_lines()
        if plot_type == 'plot_date':
            ax.set_xlabel('time')
        ax.set_ylabel(ylabel)
        plt.grid()
        
        return [fig, lines, leg]
            
        
