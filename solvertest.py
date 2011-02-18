# -*- coding: utf-8 -*-
"""
Created on Wed Feb 09 22:05:39 2011

@author: RDC
"""
from __future__ import division
import time
import os
from shutil import copyfile
import matplotlib.pyplot as plt
import numpy as np

from simman import Simulation, Simdex
from pymosim import *


# Script settings #############################################################

algos = range(1,29)
time_out = 200 # seconds, or None
run_simulations = False
analyse_logs = True
analyse_results = True
vars_to_check = [u'pq[32].vi.i.re', u'dwellings[17].pq.p']
max_perc_wrong = 2. # base 100
max_rel_rms = 1. # base 100
reference_solver = '8'      # as string!! 8 = Dassl


###############################################################################

if run_simulations:
    for algo in [str(a) for a in algos]:
        # first, create filenames for new input, result and log file
        dsin_file = ''.join(['dsin_algo_', algo, '.txt'])    
        result_file = ''.join(['Result_algo_', algo, '.mat'])
        log_file = ''.join(['dslog_algo_', algo, '.txt'])
        
        # then, create input file with the new algorithm number
        set_solver(algo, copy_to = dsin_file)
        
        # choose the appropriate dymosim (= function of the solver)
        if int(algo) < 15:
            dymosim = 'dymosim_dassl'
        else:
            dymosim = 'dymosim_rk'
        
        # start the simulation, kill process if needed and measure time        
        start_time = time.time()
        proc = run_ds(dymosim, dsin_file, result_file)
        timed_out = kill_after(proc, time_out)
        stop_time = time.time()
        copyfile('dslog.txt', log_file)

        # give some user feed back through screen dump and add info to the log
        if timed_out:
            print 'Simulation with solver ID %s timed out after %.1f seconds' \
                   % (algo, time_out)
            # if timed out, add this info to the log file
            lg = open(log_file, 'a')
            lg.write('\n\n')
            lg.write(''.join(['Corresponding result file (if any): ', 
                              result_file, '\n']))
            lg.write('This simulation timed out and was killed after %.1f seconds' \
                     % time_out)
            lg.close()
        else:
            print 'Simulation with algo_id %s took %.1f s' % \
                    (algo, stop_time - start_time)
            lg = open(log_file, 'a')
            lg.write('\n\n')
            lg.write(''.join(['Corresponding result file: ', 
                              result_file, '\n']))
            lg.close()

###############################################################################
if analyse_logs:
    # results is list with dictionaries.  Because of the for loop, this list
    # is sorted according to ascending solver ID
    results = []
    for algo in [str(a) for a in algos]:
        log_file = ''.join(['dslog_algo_', algo, '.txt'])
        summary = analyse_log(log_file)
        summary['solver ID'] = algo
        results.append(summary)
        # we create ref_sum, the dict of the refernce_solver (used in plots)        
        if algo == reference_solver:
            ref_sum = summary

# Analyse the successful simulations ##########################################

if analyse_results:
    # the principle of the data organisation is as follows
    # every simulation has a simulation summary which is a dictionary
    # ALL information that is linked to a simulation is added to this 
    # dictionary
    # All these dictionaries are grouped in lists.  There can be different
    # lists, grouping eg the summaries according to successfulness, or sorting
    # them according to a certain criteria
    
    results_successful = [s for s in results if s['successful'] == True]
    
    # now create a sim object, and add the info we need from the result file
    # to the summary dictionary
    vars_to_check.append('Time')
    for summary in results_successful:
        sim = Simulation(summary['result file'])        
        for var in vars_to_check:
            summary[var] = sim.get_value(var)
            
    res_suc_sortedCPU = sorted(results_successful, 
                          key = lambda x: float(x['CPU time']))
    
    algos = [s['solver ID'] for s in res_suc_sortedCPU]
    CPU_time = [s['CPU time'] for s in res_suc_sortedCPU]    
    perc_wrong = [s['perc_wrong'] for s in res_suc_sortedCPU]
    alg = [s['algorithm'] for s in res_suc_sortedCPU]
    algorithm = [''.join(['ID = ', ID, ' - ', a]) for ID, a in zip(algos, alg)]

    # now, make graph of ALL successfull simulations, whether they have many
    # rejected timesteps or not
        
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax2 = ax1.twinx()
    lines2 = ax2.plot(range(len(CPU_time)), perc_wrong, 'go', markersize = 5, 
             label = 'RK-only: rejected steps')
    ax2.set_ylabel('Percentage rejected steps [%]')
    ax2.set_xlim((-1, len(CPU_time)))
    ax2.set_ylim(0, 1.3 * ax2.get_ylim()[1])
    leg2 = ax2.legend(loc='upper right')
    for t in leg2.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
    
    lines1 = ax1.plot(range(len(CPU_time)), CPU_time, 'rD', markersize = 10, 
             label = 'CPU time')
    ax1.set_xlim((-1, len(CPU_time)))
    ax1.set_ylabel('CPU time in seconds')
    ticks = ax1.set_xticks(range(len(CPU_time)))
    labels = ax1.set_xticklabels(algorithm, rotation = 'vertical')
    ax1.set_xlim((-1, len(CPU_time)))
    ax1.set_ylim(0, 1.3 * ax1.get_ylim()[1])
    leg1 = ax1.legend(loc='upper left')
    for t in leg1.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
    
    # make plots of the vars_to_check compared for each solver with the reference
    # solver
    
    nb_axes = len(res_suc_sortedCPU) - 1
    if nb_axes < 4:
        nb_rows = 1
    elif nb_axes < 7:
        nb_rows = 2
    elif nb_axes < 14:
        nb_rows = 3       
    else: 
        nb_rows = 4
    nb_cols = int(np.ceil(nb_axes/nb_rows))

    for var in vars_to_check[:-1]:
        fig = plt.figure()
        nb_ax = 0
        for x in range(nb_axes + 1):
            algo_id = res_suc_sortedCPU[x]['solver ID']            
            if algo_id != reference_solver:
                nb_ax += 1
                ax = plt.subplot(nb_rows, nb_cols, nb_ax)
                refline = ax.plot(ref_sum['Time'], 
                                  ref_sum[var],
                                  color = 'r', linewidth = 2)
                cmpline = ax.plot(res_suc_sortedCPU[x]['Time'], 
                                  res_suc_sortedCPU[x][var],
                                  label = res_suc_sortedCPU[x]['algorithm'],
                                  color = 'b')        
                leg = ax.legend()
                for t in leg.get_texts():
                    t.set_fontsize('small')    # the legend text fontsize
        fig.suptitle(var)


    # next, filter further and keep only simulations according to the 
    # max_perc_wrong setting
    res_good_sortedCPU = [s for s in res_suc_sortedCPU if \
                    s['perc_wrong'] <= max_perc_wrong]

    algos = [s['solver ID'] for s in res_good_sortedCPU]
    CPU_time = [s['CPU time'] for s in res_good_sortedCPU]    
    perc_wrong = [s['perc_wrong'] for s in res_good_sortedCPU]
    alg = [s['algorithm'] for s in res_good_sortedCPU]
    algorithm = [''.join(['ID = ', ID, ' - ', a]) for ID, a in zip(algos, alg)]

    # And make exactly the same graph, only with 'good' simulations this time
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax2 = ax1.twinx()
    lines2 = ax2.plot(range(len(CPU_time)), perc_wrong, 'go', markersize = 5, 
             label = 'RK-only: rejected steps')
    ax2.set_ylabel('Percentage rejected steps [%]')
    ax2.set_xlim((-1, len(CPU_time)))
    ax2.set_ylim(0, 1.3 * ax2.get_ylim()[1])
    leg2 = ax2.legend(loc='upper right')
    for t in leg2.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
    
    lines1 = ax1.plot(range(len(CPU_time)), CPU_time, 'rD', markersize = 10, 
             label = 'CPU time')
    ax1.set_xlim((-1, len(CPU_time)))
    ax1.set_ylabel('CPU time in seconds')
    ticks = ax1.set_xticks(range(len(CPU_time)))
    labels = ax1.set_xticklabels(algorithm, rotation = 'vertical')
    ax1.set_xlim((-1, len(CPU_time)))
    ax1.set_ylim(0, 1.3 * ax1.get_ylim()[1])
    leg1 = ax1.legend(loc='upper left')
    for t in leg1.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
    
  
    # make plots of the vars_to_check compared for each solver with the reference
    # solver
    
    nb_axes = len(res_good_sortedCPU) - 1
    if nb_axes < 4:
        nb_rows = 1
    elif nb_axes < 7:
        nb_rows = 2
    elif nb_axes < 14:
        nb_rows = 3       
    else: 
        nb_rows = 4
    nb_cols = int(np.ceil(nb_axes/nb_rows))

    for var in vars_to_check[:-1]:
        fig = plt.figure()
        nb_ax = 0
        for x in range(nb_axes + 1):
            algo_id = res_good_sortedCPU[x]['solver ID']            
            if algo_id != reference_solver:
                nb_ax += 1
                ax = plt.subplot(nb_rows, nb_cols, nb_ax)
                refline = ax.plot(ref_sum['Time'], 
                                  ref_sum[var],
                                  color = 'r', linewidth = 2)
                cmpline = ax.plot(res_good_sortedCPU[x]['Time'], 
                                  res_good_sortedCPU[x][var],
                                  label = res_good_sortedCPU[x]['algorithm'],
                                  color = 'b')        
                leg = ax.legend()
                for t in leg.get_texts():
                    t.set_fontsize('small')    # the legend text fontsize
        fig.suptitle(var)
    
    
    # finally, check RMS error for the given vars in vars_to_check
    # strangely, some simulations have no output on the last timestep
    min_length = min([len(s['Time']) for s in res_good_sortedCPU])
    # check length of our ref_sum
    for var in vars_to_check:
        if len(ref_sum[var]) > min_length:
            # this simulation has more timesteps than the minimum, slice
            ref_sum[var] = ref_sum[var][:min_length]
        if var != 'Time':
                # compute reference RMS value of the signal and add it to the summary
                rmsref = np.sqrt(ref_sum[var]**2).sum()
                ref_sum[var + '_rmsref'] = rmsref
    # compute rms of the difference and also relatively to the reference_solve
    # the relative RMS is expressed as percentage (base 100)
    for summary in res_good_sortedCPU:
        for var in vars_to_check:
            if len(summary[var]) > min_length:
                # this simulation has more timesteps than the minimum, slice
                summary[var] = summary[var][:min_length]
            if var != 'Time':
                # compute RMS value and add it to the summary
                rmsdif = np.sqrt((summary[var] - ref_sum[var])**2).sum()
                try:
                    rmsrel = 100 * rmsdif / ref_sum[var + '_rmsref']
                except ZeroDivisionError:
                    rmsrel = rmsdif
                summary[var + '_rmsdif'] = rmsdif
                summary[var + '_rmsrel'] = rmsrel
    
    # and now create an overview plot again
    RMS_rel = {}
    for var in vars_to_check[:-1]:
        cmd = ''.join(["[s['", var, '_rmsrel', "'] for s in res_good_sortedCPU]"])
        RMS_rel[var] = eval(cmd)
    
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax2 = ax1.twinx()
    for var in RMS_rel:
        # var is the key
        lines2 = ax2.plot(range(len(CPU_time)), RMS_rel[var], 'o', 
                          markersize = 5, 
                          label = ''.join(['Relative RMS for ', var]))
    ax2.set_ylabel('Relative RMS [%]')
    ax2.set_xlim((-1, len(CPU_time)))
    ax2.set_ylim(0, 1.3 * ax2.get_ylim()[1])
    leg2 = ax2.legend(loc='upper right')
    for t in leg2.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
    
    lines1 = ax1.plot(range(len(CPU_time)), CPU_time, 'rD', markersize = 10, 
             label = 'CPU time')
    ax1.set_xlim((-1, len(CPU_time)))
    ax1.set_ylabel('CPU time in seconds')
    ticks = ax1.set_xticks(range(len(CPU_time)))
    labels = ax1.set_xticklabels(algorithm, rotation = 'vertical')
    ax1.set_xlim((-1, len(CPU_time)))
    ax1.set_ylim(0, 1.3 * ax1.get_ylim()[1])
    leg1 = ax1.legend(loc='upper left')
    for t in leg1.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
    
    # next, filter further and keep only simulations according to the 
    # max_rel_rms setting
    res_perfect_sortedCPU = [s for s in res_good_sortedCPU if \
                    s[vars_to_check[0] + '_rmsrel'] <= max_rel_rms]

    algos = [s['solver ID'] for s in res_perfect_sortedCPU]
    CPU_time = [s['CPU time'] for s in res_perfect_sortedCPU]    
    perc_wrong = [s['perc_wrong'] for s in res_perfect_sortedCPU]
    alg = [s['algorithm'] for s in res_perfect_sortedCPU]
    algorithm = [''.join(['ID = ', ID, ' - ', a]) for ID, a in zip(algos, alg)]
    RMS_rel = {}
    for var in vars_to_check[:-1]:
        cmd = ''.join(["[s['", var, '_rmsrel', "'] for s in res_perfect_sortedCPU]"])
        RMS_rel[var] = eval(cmd)
        
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax2 = ax1.twinx()
    for var in RMS_rel:
        # var is the key
        lines2 = ax2.plot(range(len(CPU_time)), RMS_rel[var], 'o', 
                          markersize = 5, 
                          label = ''.join(['Relative RMS for ', var]))
    ax2.set_ylabel('Relative RMS [%]')
    ax2.set_xlim((-1, len(CPU_time)))
    ax2.set_ylim(0, 1.3 * ax2.get_ylim()[1])
    leg2 = ax2.legend(loc='upper right')
    for t in leg2.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
    
    lines1 = ax1.plot(range(len(CPU_time)), CPU_time, 'rD', markersize = 10, 
             label = 'CPU time')
    ax1.set_xlim((-1, len(CPU_time)))
    ax1.set_ylabel('CPU time in seconds')
    ticks = ax1.set_xticks(range(len(CPU_time)))
    labels = ax1.set_xticklabels(algorithm, rotation = 'vertical')
    ax1.set_xlim((-1, len(CPU_time)))
    ax1.set_ylim(0, 1.3 * ax1.get_ylim()[1])
    leg1 = ax1.legend(loc='upper left')
    for t in leg1.get_texts():
        t.set_fontsize('small')    # the legend text fontsize


    # finally finally, make plots of the vars_to_check compared for each solver 
    # with the reference solver for all 'perfect' simulations
    
    nb_axes = len(res_perfect_sortedCPU) - 1
    if nb_axes < 4:
        nb_rows = 1
    elif nb_axes < 7:
        nb_rows = 2
    elif nb_axes < 14:
        nb_rows = 3       
    else: 
        nb_rows = 4
    nb_cols = int(np.ceil(nb_axes/nb_rows))

    for var in vars_to_check[:-1]:
        fig = plt.figure()
        nb_ax = 0
        for x in range(nb_axes + 1):
            algo_id = res_perfect_sortedCPU[x]['solver ID']            
            if algo_id != reference_solver:
                nb_ax += 1
                ax = plt.subplot(nb_rows, nb_cols, nb_ax)
                refline = ax.plot(ref_sum['Time'], 
                                  ref_sum[var],
                                  color = 'r', linewidth = 2)
                cmpline = ax.plot(res_perfect_sortedCPU[x]['Time'], 
                                  res_perfect_sortedCPU[x][var],
                                  label = res_perfect_sortedCPU[x]['algorithm'],
                                  color = 'b')        
                leg = ax.legend()
                for t in leg.get_texts():
                    t.set_fontsize('small')    # the legend text fontsize
        fig.suptitle(var)
    
    plt.show()


"""
... "climateDataHourly.mat" loading (tables for interpolation)
... "climateData.mat" loading (tables for interpolation)
Model: TME.WIP.Quarter
Unknown method in Godess. Number: 16
Known ones:
15  - radau IIa (3 stages, order 5) order 3 dense output, fully implicit
16  - esdirk23a(4 stages, order 3) dense output, singly implicit
17 - esdirk34a(5 stages, order 4) dense output, singly implicit
18  - esdirk45a(7 stages, order 5) dense output, singly implicit
19  - dopri45  (7 stages, order 5) order 4 dense output, explicit
20  - dopri78 (13 stages, order 8) no dense output, explicit
21  - dopri853(16 stages, order 8) dense output, explicit
23  - sharp67 (12 stages, order 7) no dense output, explicit
24  - sdirk34hw (5 stages, order 4) order 3 dense output, singly implicit
25  - sdirk34var  (5 stages, order 4) order 3 dense output, singly implicit
26 - cerk23   (4 stages, order 3) dense output, explicit
27 - cerk34   (6 stages, order 4) dense output, explicit
28 - cerk45   (8 stages, order 5) dense output, explicit
 Dense output is required for for finding events.
 The order of the dense output is given if less than the usual order, it indicates that
 the error is slightly larger than normal for events (minor loss of accuracy).
 For stiff problems implicit methods must be used.
  Singly implicit methods allow successive solution of the stages (more stages)
  whereas fully implicit has to solve all stages simulatenously (more linear algebra).
 Higher order methods are more suited for strict tolerances.
Specify method from the above

"""


"""
Allowed integration methods are:

        |         | model |       |        | dense | state |
method  | imethod | typ   | stiff | order  | output| event |
--------+---------+-------+-------+--------+-------+-------+
deabm   |    1    |   ode |   no  |  1-12  |  yes  |   no  |
lsode1  |    2    |   ode |   no  |  1-12  |  yes  |   no  |
lsode2  |    3    |   ode |  yes  |  1-5   |  yes  |   no  |
lsodar  |    4    |   ode |  both |1-12,1-5|  yes  |  yes  |
dopri5  |    5    |   ode |   no  |   8    |   no  |   no  |
dopri8  |    6    |   ode |   no  |   5    |   no  |   no  |
grk4t   |    7    |   ode |  yes  |   4    |   no  |   no  |
dassl   |    8    |   dae |  yes  |  1-5   |  yes  |  yes  |
odassl  |    9    |  hdae |  yes  |  1-5   |  yes  |  yes  |
mexx    |   10    |  hdae |   no  |  2-24  |   no  |   no  |
euler   |   11    |   ode |   no  |   1    |   no  |  yes  |
rkfix2  |   12    |   ode |   no  |   2    |   no  |  yes  |
rkfix3  |   13    |   ode |   no  |   3    |   no  |  yes  |
rkfix4  |   14    |   ode |   no  |   4    |   no  |  yes  |
--------+---------+-------+-------+--------+-------+-------+
"""










"""
   #             | model|       |        | dense | state |
   # Algorithm   | typ  | stiff | order  | output| event |
   # ------------+------+-------+--------+-------+-------+
   #  1 | deabm  |  ode |   no  |  1-12  |  yes  |   no  |
   #  2 | lsode1 |  ode |   no  |  1-12  |  yes  |   no  |
   #  3 | lsode2 |  ode |  yes  |  1-5   |  yes  |   no  |
   #  4 | lsodar |  ode |  both |1-12,1-5|  yes  |  yes  |
   #  5 | dopri5 |  ode |   no  |   5    |   no  |   no  |
   #  6 | dopri8 |  ode |   no  |   8    |   no  |   no  |
   #  7 | grk4t  |  ode |  yes  |   4    |   no  |   no  |
   #  8 | dassl  |  dae |  yes  |  1-5   |  yes  |  yes  |
   #  9 | odassl | hdae |  yes  |  1-5   |  yes  |  yes  |
   # 10 | mexx   | hdae |   no  |  2-24  |   no  |   no  |
   # 11 | euler  |  ode |   no  |   1    |   no  |  yes  |
   # 12 | rkfix2 |  ode |   no  |   2    |   no  |  yes  |
   # 13 | rkfix3 |  ode |   no  |   3    |   no  |  yes  |
   # 14 | rkfix4 |  ode |   no  |   4    |   no  |  yes  |
   #>=14| others |  ode |yes/no |  2-5   |   yes |  yes  |
   # ---+--------+------+-------+--------+-------+-------+
   # euler and rkfix have fixed stepsize.

"""

"""
To do:
    - check if the simulation was successful
    - get all the resultfiles in a simdex or simulation objects and plot something
"""
