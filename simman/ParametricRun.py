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
setup = False
run_it = False   
get_results = True
analyse_results = True

if setup:
    # number of cpu's for the parametric run
    ncpus=3
    
    # directory of the parametric run.  Make sure it contains a dsin.txt and 
    # dymosim.exe and a subfolder 'inputs_pr' with all input files 
    work_dir = 'C:\\Workspace\\Modelica\\DeVloei\\Work\\ParRun'
    
    # Specify here the parameters and values for the parametric run
    # format: list of ['long_name', [val1, ..., valn], 'short_name']
    pars = [
    #  [r'building_forGrid.heaSys.onOffDelay', [900, 2000], 'onOffDelay'],
      [r'building_forGrid.heaSys.betaFactorHeatPump', [0.5,0.8], 'beta'],
      [r'building_forGrid.whichUser',[3,30], 'user'],
      [r'building_forGrid.heaSys.volumeTank',[0.2,0.3, 0.4], 'vol_tank'],
      [r'building_forGrid.heaSys.HPControl.dTSafetyTop',[3,4,5], 'dTSafetyTop']
            ]
    # Specify the variables required for the analysis
    variables = dict(TMixed = u'building_forGrid.heaSys.dHW.TMixed',
                     m_flowTotal = u'building_forGrid.heaSys.dHW.m_flowTotal',
                     PEl = u'building_forGrid.heaSys.heatPump.PEl' )
                     
    """
    Remark: the x- and y-axis for the plots can be defined in the 'analyse_results'
    section (you can search for to_plot_x to get there).
    """ 

###############################################################################
os.chdir(os.path.abspath(work_dir))

if run_it:
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
        pars_to_run[par[0]] = np.array(par[1])
    
    # now we set everything ready for parallel python, the return consistst of a 
    # list of cd-paths where all sub-sets are located as for parallel python
    sub_dir = pymosim.set_parametric_run(work_dir, pars_to_run)
    
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
    for r in results:          
        # analyse what you want and add results to the dictionary        
        TDHW_too_low = np.nonzero(r['TMixed'] < (273.15+45.))[0]
        m_flow_discomfort = r['m_flowTotal'][TDHW_too_low]
        m_flow_discomfort_weighted = (r['m_flowTotal'] * (273.15+45.-r['TMixed']))[TDHW_too_low]
        DHW_discomfort = sum(m_flow_discomfort)/sum(r['m_flowTotal'])
        dT_discomfort = sum(m_flow_discomfort_weighted)/sum(m_flow_discomfort)
        r['DHW_discomfort'] = DHW_discomfort
        r['dT_discomfort'] = dT_discomfort
        r['TDHW_min'] = min(r['TMixed'])
        r['PEl_total'] = np.trapz(r['PEl'], r['Time'])/1e6 # MJ
    

    # make numpy arrays with the values of the things we'd like to plot
    # these arrays are named according to the variable they contain!!
    # to make the script safe, the arrays are added as attributes to a dummy 
    # class called arrays. 
    class Dummy():
        pass
    arrays = Dummy()
    
    # x and y axis of the plots
    to_plot_x = ['DHW_discomfort', 'dT_discomfort', 'CPU_time']
    to_plot_y = ['PEl_total', 'dT_discomfort']
    
    to_plot = []
    to_plot.extend(to_plot_x)
    to_plot.extend(to_plot_y)    
    # beside results, we also need all parameters
    to_plot.extend(parameters.keys())
        
    for k in to_plot:
        ar = np.array([np.double(d[k]) for d in results])
        exec(''.join(['arrays.', k, '= ar']))        
        
    # Make the graph(s) ########################################################   
   
    nb_axes = len(parameters)
    if nb_axes < 4:
        nb_rows = 1
    elif nb_axes < 7:
        nb_rows = 2
    elif nb_axes < 14:
        nb_rows = 3       
    else: 
        nb_rows = 4
    nb_cols = int(np.ceil(nb_axes/nb_rows))

    plots_to_make = [(x,y) for x in to_plot_x for y in to_plot_y]
    for x,y in plots_to_make:
        # make 1 plot for each x-axis / y-axis combination        
        array_x = eval(''.join(['arrays.', x]))
        array_y = eval(''.join(['arrays.', y]))
        fig = plt.figure()
        for a, par in enumerate(parameters):
            # make the subplots, 1 per parameter variation
            ax = plt.subplot(nb_rows, nb_cols, a+1)
            all_res = ax.plot(array_x, array_y,
                              color = '0.5', linestyle = '', marker = 'o')
            for par_value in pars_to_run[parameters[par]]:
                # make the colored subplots for each value of the parameter                
                rel_diff = (eval(''.join(['arrays.', par])) - par_value)/par_value            
                index = np.nonzero(np.abs(rel_diff) < 0.001)[0]            
                array_x_sel = array_x[index]
                array_y_sel = array_y[index]           
                ax.plot(array_x_sel, array_y_sel, linestyle = '', marker = 'o',
                        label = ''.join([par, '=', str(par_value)]))
            ax.set_ylim(top = 1.2 * ax.get_ylim()[1])
            
            leg = ax.legend(loc = 'upper right')
            for t in leg.get_texts():
                t.set_fontsize('small')    # the legend text fontsize
            ax.set_xlabel(x)
            ax.set_ylabel(y)
        fig.suptitle(''.join([y, ' versus ', x]))

# free the results_pr directory    
os.chdir(os.path.abspath(work_dir))