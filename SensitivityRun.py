# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 11:42:03 2011

@author: RBa
"""

from __future__ import division
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pp
from shutil import rmtree

import pymosim
from simman import Simulation, Simdex

sys.path.append(os.path.abspath(r'D:\Ruben_BWK239\GIT_python'))

# Script settings #############################################################
setup_and_run = False   
get_results = True
analyse_results = True

# number of cpu's for the sensitivity run
ncpus=3

# directory of the sensitivity run.  Make sure it contains a dsin.txt and 
# dymosim.exe and a subfolder 'inputs_pr' with all input files 
work_dir = r'C:\Workspace\BS2011\PR2'

# the sensitivity for the sensitivity analysis. For all variables, a 
# run will be executed with (1-sensitivity) * default and 
# (1+sensitivity) * default.   
sensitivity = 0.1

# Specify here the parameters and values for the sensitivity run
# format: list of ['long_name', reference value, 'short_name']
# attention: for integers, the values will be rounded to the nearest integer
# and for default integers smaller than 1/sensitivity, the sensitivity will 
# be studied taking default +/ 1
pars = [
    [r'building_forGrid.heaSys.betaFactorHeatPump',0.8, 'beta'],
#    [r'building_forGrid.heaSys.volumeTank', 0.25, 'vol_tank'],
#    [r'building_forGrid.heaSys.HPControl.dTSafetyTop', 3, 'dTSafetyTop'],
    [u'building_forGrid.heaSys.timeFilter', 43200, 'timeFilter'],
    [u'building_forGrid.heaSys.dTSupRetNom', 6, 'dTSupRetNom'],
    [u'building_forGrid.heaSys.FHChars[1].T', 0.2, 'FHChars_T1'],
    [u'building_forGrid.heaSys.FHChars[2].T', 0.2, 'FHChars_T2']
        ]
# Specify the variables required for the analysis
variables = dict(TMixed = u'building_forGrid.heaSys.TDHW',
                 mDHW = u'building_forGrid.heaSys.mDHW',
                 PEl = u'building_forGrid.heaSys.P',
                 QHP = u'building_forGrid.heaSys.QHP',
                 QHeatTotal = 'building_forGrid.heaSys.QHeatTotal')
                 
"""
Remark: the x- and y-axis for the plots can be defined in the 'analyse_results'
section (you can search for to_plot_x to get there).
""" 

###############################################################################
os.chdir(os.path.abspath(work_dir))

if setup_and_run:
    if os.path.exists(os.path.join(work_dir, 'results_pr')):
        remove_results = raw_input('Remove the results in results_pr ? (y/n) \n ==> ')
        if remove_results == 'y' or remove_results == 'Y':
            rmtree(os.path.join(work_dir, 'results_pr'))
        else:
            raise NotImplementedError("Remove the results_pr folder first !")

    # first we must make a dictionary of the parameteric we want to do
    # and we add the parameters and their values we want to run
    pars_to_run = {}
    for par in pars:
        pars_to_run[par[0]] = par[1]
    
    # now we set everything ready for parallel python, the return consists of a 
    # list of cd-paths where all sub-sets are located as for parallel python
    sub_dir = pymosim.set_sensitivity_run(work_dir, pars_to_run)
    
    # initialise parallel python
    ppservers=()
    job_server = pp.Server(ppservers = ppservers, ncpus=ncpus)
    # define the inputs and start the simulations
    inputs = tuple(sub_dir)
    jobs = [(input, job_server.submit(pymosim.start_parametric_run, args = (input,), modules = ("subprocess","os"))) for input in inputs]
    job_server.wait() 
    job_server.destroy()
    
    pymosim.close_parametric_run(work_dir, sub_dir)

###############################################################################
if get_results:
    print 'Getting results for:'
    # make a list of logfiles
    os.chdir(os.path.join(work_dir, 'results_pr'))    
    fileList = os.listdir(os.getcwd())
    log_files = [f for f in fileList if f.find('dslog_run_') > -1]
    # mat_files = [f for f in fileList if f.find('result_run_') > -1]
    
    # results is list with dictionaries containing all useful info from the runs
    results = []
    
    # first, get info from the log files    
    for l in log_files:
        summary = pymosim.analyse_log(l)
        summary['mat_file'] = l.replace('.txt','.mat').replace('dslog', 'result')
        results.append(summary)
    # suc6 = [l['successful'] for l in logs]    checked, all are successful!  
    
    # second, get the used parameters and simulation results
    
    # paramaeteres = contains short and long names of the parameters we need.
    # The data will be extracted from the .mat files and put in a dictionary
    # The short names will be the keys in that dictionary     
    parameters = {}
    for par in pars:
        parameters[par[2]] = par[0]
    
    for r in results:
        print ''.join(['\t - ', r['mat_file']])
        # r is the dictionary containing all results for this specific run
        sim = Simulation(r['mat_file'])
         # we add info to the dictionary: key = var, values = value for var        
        for k,v in parameters.iteritems():
            r[k] = sim.get_value(v)
        for k,v in variables.iteritems():
            r[k] = sim.get_value(v)
        # we also get Time 
        r['Time'] = sim.get_value('Time')

###############################################################################        
if analyse_results:  
    print 'Analyzing results for:'      
    for r in results:          
        print ''.join(['\t - ', r['mat_file']])
        # DHW analysis
        TDHW_too_low = np.nonzero(r['TMixed'] < (273.15+45.))[0]
        m_flow_discomfort = r['mDHW'][TDHW_too_low]
        m_flow_discomfort_weighted = (r['mDHW'] * (273.15+45.-r['TMixed']))[TDHW_too_low]
        DHW_discomfort = sum(m_flow_discomfort)/sum(r['mDHW'])
        dT_discomfort = sum(m_flow_discomfort_weighted)/sum(m_flow_discomfort)
        r['DHW_discomfort'] = DHW_discomfort
        r['dT_discomfort'] = dT_discomfort
        r['TDHW_min'] = min(r['TMixed'])
        
        # Heating analysis
        # this could be more efficient: getting the results should happen in the 
        # previous section in order to avoid recreating sim objects
        total = 0
        sim = Simulation(r['mat_file'])
        for z in range(int(sim.get_value(u'building_forGrid.heaSys.n_C'))):
            var_name = ''.join(['discomfort_zone_', str(z+1)])
            dt = sim.get_value(''.join([u'building_forGrid.TopAsked[', str(z+1), u']'])) \
                 - sim.get_value(''.join([u'building_forGrid.Top[', str(z+1), u']']))
            dt[np.nonzero(dt<0)]=0
            dt_int = np.trapz(dt, r['Time'])/3600 # in Kh
            r[var_name] = dt_int.sum()
            total += r[var_name]
        r['Heating_discomfort'] = total    
        
        # Energy consumption and SPF
        r['PEl_yr'] = np.trapz(r['PEl'], r['Time'])/1e6 # MJ
        r['QHeatTotal_yr'] = np.trapz(r['QHeatTotal'], r['Time'])/1e6 # MJ
        r['QHP_yr'] = np.trapz(r['QHP'], r['Time'])/1e6 # MJ
        r['SPF_HP'] = r['QHP_yr'] / r['PEl_yr']
        r['SPF_system'] = r['QHeatTotal_yr'] / r['PEl_yr']
    
    # Make the graph(s) ########################################################
    
    # make numpy arrays with the values of the things we'd like to plot
    # these arrays are named according to the variable they contain!!
    # to make the script safe, the arrays are added as attributes to a dummy 
    # class called arrays. 
    
    class Dummy():
        pass
    arrays = Dummy()
    
    # x and y axis of the plots
    to_plot_x = ['DHW_discomfort', 'Heating_discomfort', 'CPU_time']
    to_plot_y = ['PEl_yr', 'SPF_system']
    
    to_plot = []
    to_plot.extend(to_plot_x)
    to_plot.extend(to_plot_y)    
    # beside results, we also need all parameters
    to_plot.extend(parameters.keys())
        
    for k in to_plot:
        ar = np.array([np.double(d[k]) for d in results])
        exec(''.join(['arrays.', k, '= ar']))        
        
    # setting the colors, 1 color for each parameter
    cm = plt.cm.jet
    colors = [cm.__call__(f) for f in np.arange(0,1,1./len(parameters))]
   
    plots_to_make = [(x,y) for x in to_plot_x for y in to_plot_y]
    for x,y in plots_to_make:
        # make 1 plot for each x-axis / y-axis combination        
        array_x = eval(''.join(['arrays.', x]))
        array_y = eval(''.join(['arrays.', y]))
        fig = plt.figure()
        ax = plt.subplot(111)
        
        color_index = 0
        for k, v in parameters.items():
            c = colors[color_index]
            # plot the sensitivity results for each parameter
            # 1. get x and y values
            par_values = eval(''.join(['arrays.', k]))
            # index should become a list with 3 indices pointing to 
            # the lower value, ref value and upper value for the considered parameter            
            try:
                index = [np.nonzero(par_values < pars_to_run[v])[0][0]]
            except(IndexError):
               index = [0]                    
            index.append(0)
            try:             
                index.append(np.nonzero(par_values > pars_to_run[v])[0][0])
            except(IndexError):
                index.append(0)
            
            array_x_sel = array_x[index]
            array_y_sel = array_y[index] 
            
            # 2. Plot lower, ref and upper as a colored line, no markers
            ax.plot(array_x_sel, array_y_sel, color=('0.8'))

            # 3. Plot lower as a hollow round and upper as solid round
            ax.plot(array_x_sel[0], array_y_sel[0], linestyle = 'None', 
                    marker = 'o', markersize = 10,  
                    markeredgecolor = c, markerfacecolor = c, label = k)
            ax.plot(array_x_sel[2], array_y_sel[2], linestyle = 'None', 
                    marker = '*', markersize = 13, markerfacecolor = c,
                    markeredgecolor = c)   
            
            color_index += 1
                    
        # 4. Plot the ref? Not needed maybe: = crossing of the lines
        ylim = ax.get_ylim()
        ax.set_ylim(top = ylim[1] + 1.2 * (ylim[1] - ylim[0]))
            
        leg = ax.legend(loc = 'upper right')
        for t in leg.get_texts():
            t.set_fontsize('small')    # the legend text fontsize
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        fig.suptitle(''.join([y, ' versus ', x]))

# free the results_pr directory    
os.chdir(os.path.abspath(work_dir))