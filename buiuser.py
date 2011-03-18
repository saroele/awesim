# -*- coding: utf-8 -*-
"""
Created on Mon Mar 07 09:40:50 2011

@author: RBa

This module contains tools to:
    1) smoothen user behaviour to time and create a new input file for sims
    based on simulated 1-min stochastic data in BWF.mo en BWFlib.mo
    2) analyse the relation between the depicted (determinstistic or stochatsic)
    user behaviour and the stochatsic electricity production, e.g. by PV

Most important methods:
    - smoothen: smoothen data according to a certain timestep and create a new 
    file to read in the new data
    - ...
"""

# import os
# import time
import numpy as np
# import matplotlib.pyplot as plt
# from simman import Simulation, Simdex

def get_smoothened_data(from_data, from_time, to_time, to_file):
    """
    get_smoothened_data(data,from,to)
    
    Smoothen the data of a certain datainput to a new timestep and write this
    to a .txt-file for direct use as input of dymola-simulations
    - from_data = np.array of data through time
    - from_time = int in seconds with the current time step or np.array of time
    - to_time = int in seconds whith the future time step or np;.array of time
    The function returns
    - to_file = .txt-file of data through time for reading in dymola
    """
    # first we check if the provided data has the correct structure to smoothen
    # secondly we define the current timestep of the provided data in data_time
    try:
        if type(from_data) == np.array or type(from_data) == np.ndarray or \
        type(from_data) == list:
            if type(from_time) == np.array or type(from_time) == np.ndarray \
            or type(from_time) == list:
                if type(to_time) == int and abs(to_time)>0:
                    from_step = abs(from_time[1]-from_time[0])
                    to_step = abs(to_time)
                    print 'The data is a %s with a step of %s s to be \
                    converted to a time step of %s s.' \
                    %(type(from_data),from_step,to_step)
                elif type(to_time) == np.array or type(to_time) == np.ndarray \
                or type(to_time) == list:
                    from_step = abs(from_time)
                    to_step = abs(to_time[1]-to_time[0])                 
                    print 'The data is a %s with a step of %s s to be \
                    converted to a time step of %s s.' \
                    %(type(from_data),from_step,to_step)
    except:
        raise TypeError('The input data is incorrect.')

    # if to_time is an integer of the time step, to_time is transferred to an 
    # new array consisting the actual time line of the new data
    if type(to_time) == int:
        ratio = to_step / from_step
        if ratio < 1:
            raise TypeError('One can not smoothen to a smaller timestep ...')
        size = 1 + (len(from_data) - 1) / ratio    
        to_time = np.zeros(size)
        to_time[0] = from_time[0]
        for i in range(1, len(to_time)):
            to_time[i] = from_time[i*ratio]
    
    # now get the new smoothened data by means of a simple moving average
    to_data = []
    to_data = moving_average(from_data, from_step, to_step)    

    # and write the new data to a .txt-file with two colums (time,data)
    filename = ''.join([to_file, '.txt'])
    text_file = open(filename, 'w')
    for line in range(len(to_data)):
        text_file.write('%s %s \n' %(to_time[line], to_data[line]))
    text_file.close()
    
    return (to_time, to_data)

def moving_average(from_data, from_step, to_step):
    """
    moving_average(data,from,to)
    
    Smoothen the data of a certain datainput to a new timestep
    - from_data = np.array of data through time
    - from_step = int in seconds with the current time step
    - to_step = int in seconds whith the future time step
    The function returns
    - to_data = np.array of data through time
    """
    
    # first we check if the provided data has the correct structure to smoothen
    # secondly we define the current timestep of the provided data in data_time
    try:
        if type(from_data) == np.array or type(from_data) == np.ndarray or \
        type(from_data) == list:
            if type(from_step) == int and abs(from_step)>0:
                if type(to_step) == int and abs(to_step)>0:
                    print 'The input data is correct.'
    except:
        raise TypeError('The input data is incorrect.')

    # the ratio of smoothening is defined    
    ratio = to_step / from_step
    if ratio < 1:
        raise TypeError('One can not smoothen to a smaller timestep ...')
    
    # an empty list is declared for storage of the smoothened data    
    size = 1 + (len(from_data) - 1) / ratio    
    to_data = np.zeros(size)

    # the list to_data is filled by the averaged data of from_data
    to_data[0] = from_data[0]
    for i in range(1, len(to_data)):
        data = from_data[i*ratio-ratio+1:i*ratio+1]
        to_data[i] = np.average(data)
    
    return to_data


