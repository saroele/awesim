# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 08:51:02 2011

@author: RBa
"""

__author__ = "Ruben Baetens"
__version__ = "0.0.1"



import matplotlib.pyplot as plt
import numpy as np
from buiuser import *
import re
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import mpl_toolkits.axisartist as AA
from simman import Simulation

def plt_comfort(self):
    """
    plt_comfort(self)
        
    plot a graphical unit of the indoor temperatures of #zones
    - self = simman.Simulation-object
    """

    # now get the data of the zone temperatures
    time = self.get_value(u'Time')
    print 'Time data,'

    # now get the data of the outdoor temperature
    dickie = self.exist(u'sim.Te')
    outdoor = self.get_value(dickie[-1]) - 273.15
    print 'the outdoor temperature'
            
    # now get the data of the zone temperatures
    dickie = self.exist(u'summary.Top')
    temps = np.zeros((len(time),len(dickie)))
    for colum in range(len(dickie)):            
        temps[:,colum] = self.get_value(dickie[colum])
    temps = temps - 273.15
    print 'and %s indoor temperatures are succesfully found in the Sim.object' \
    %(len(dickie))

    # now get the statistical data of the temperatures ...
    median, p_stdev, m_stdev, ppp_stdev, mmm_stdev = get_stdev(time, temps)

    to_time, outdoor = get_smoothened_data(outdoor, time, 86400, 'to_file')
    to_time, median = get_smoothened_data(median, time, 86400, 'to_file')
    to_time, p_stdev = get_smoothened_data(p_stdev, time, 86400, 'to_file')
    to_time, m_stdev = get_smoothened_data(m_stdev, time, 86400, 'to_file')
    to_time, ppp_stdev = get_smoothened_data(ppp_stdev, time, 86400, 'to_file')
    to_time, mmm_stdev = get_smoothened_data(mmm_stdev, time, 86400, 'to_file')

    # ... and plot them
    gu_yeardata(to_time, outdoor, median, p_stdev, m_stdev, ppp_stdev, mmm_stdev)

def get_stdev(time, temps):
    """
    get_stdev(time, temps)
    
    Get the standard deviations sigma and 3*sigma from a np.array
    - time = reference np.array of the timeline of the provided data
    - temps = np.array of data from which the standard deviation has to be known
    The function returns
    - median, p_stdev, m_stdev, ppp_stdev, mmm_stdev
    """

    try: 
        if len(time) == np.shape(temps)[0]:
            print 'stdev will be calculated'
        elif len(time) == np.shape(temps)[1]:
            temps = temps.T            
            print 'stdev will be calculated though after being transponed'
    except:
        raise IOError
    
    # basic definition of standard deviation figures
    one_std = int(np.shape(temps)[1] * 0.682689492137 / 2)
    three_std = int(np.shape(temps)[1] * 0.997300203937 / 2)
    print 'one_std is %s and three_std is %s' %(one_std, three_std)
    # sort all provided data in a new array with the same dimensions
    sorted_temps = np.zeros((np.shape(temps)[0],np.shape(temps)[1]))    
    for line in range(np.shape(temps)[0]):
        sorted_line = np.sort(temps[line, :])
        sorted_temps[line,:] = sorted_line
    # get the deviations of the provided data as five np.arrays
    median = sorted_temps[:, np.shape(temps)[1]/2]
    p_stdev = sorted_temps[:, np.shape(temps)[1]/2+one_std]
    m_stdev = sorted_temps[:, np.shape(temps)[1]/2-one_std]
    ppp_stdev = sorted_temps[:, np.shape(temps)[1]/2+three_std]
    mmm_stdev = sorted_temps[:, np.shape(temps)[1]/2-three_std]
    
    return median, p_stdev, m_stdev, ppp_stdev, mmm_stdev   

def gu_yeardata(time, ref, median, p_stdev, m_stdev, ppp_stdev, mmm_stdev):
    """
    gu_buitemp(time, outdoor, temps)
    
    plot a graphical unit of the indoor temperatures of #zones
    - time = reference np.array of the timeline of the provided data
    - outdoor = np.array of temperature
    - temps = np.array of data from which the plot is wanted
    """

    if time[-1] - time [0] > 10000000:
        xticks = (0, 2678400, 5097600, 7776000, 10368000, 13046400, \
        15638400, 18316800, 20995200, 23587200, 26265600, 28857600, 31536000)
        xticknames =  ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', \
        'Aug', 'Sep', 'Oct', 'Nov', 'Dec', '')
    
    fig = plt.figure(1, figsize=(10, 4))
        # graphical unit for the large graph
    gu1 = AA.Subplot(fig, 1, 1, 1)
    fig.add_subplot(gu1)
    gu1.set_ylabel('Operative temperature, Celsius')
    gu1.set_ylim( (-10, 35) )
    gu1.set_yticks( (-10, -5, 0, 5, 10, 15, 20, 25, 30, 35) )
    # then we put the x-axes right
    # here it is possible to do _xmajortick and _xminortick
    gu1.axis['bottom', 'top', 'right'].set_visible(False)
    gu1.axis['timeline'] = gu1.new_floating_axis(nth_coord=0, \
    value=0, axis_direction = 'bottom')
    gu1.axis['timeline'].toggle(all=True)
    gu1.set_xlim( (0, 31536000) )
    gu1.set_xticks(xticks)
    gu1.set_xticklabels(xticknames, ha = 'left', size = 'small', rotation = 45)
    # gu_history.set_xminorticks[(-0, 2678400, 5097600, 7776000)]
    
    # then we put the title correct
    gu1.set_title('Temperature', ha = 'right', position = (0.17, 1))
    # then we plot the data
    p1 = gu1.plot(time, ref, color='red', lw = 1)
    p2 = gu1.plot(time, median, '--', color='black', lw = 1)
    p3 = gu1.plot(time, p_stdev, color='grey', lw = 1)
    p4 = gu1.plot(time, m_stdev, color='grey', lw = 1)
    gu1.plot(time, ppp_stdev, color='lightgrey', lw = 1)
    gu1.plot(time, mmm_stdev, color='lightgrey', lw = 1)
    # then we fill the plotted data
    gu1.fill_between(time, p_stdev, m_stdev, where = p_stdev >= m_stdev, \
    facecolor = 'grey', interpolate = True)    
    gu1.fill_between(time, ppp_stdev, p_stdev, where = ppp_stdev >= p_stdev, \
    facecolor = 'lightgrey', interpolate = True)    
    gu1.fill_between(time, m_stdev, mmm_stdev, where = m_stdev >= mmm_stdev, \
    facecolor = 'lightgrey', interpolate = True)  
    
    gu1.legend( (p1), ('Outdoor'), 'lower left')    
    gu1.legend( (p2), ('median'), 'lower center')    
    gu1.legend( (p3, p4), ('sigma','3*sigma'), loc = (0,-0.2))    
    
    fig.show()
