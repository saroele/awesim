# -*- coding: utf-8 -*-
"""
Created on Wed Feb 09 22:05:39 2011

@author: RDC
"""
from __future__ import division
import time
import os
from subprocess import Popen
from simman import Simulation, Simdex
from shutil import copyfile



def set_solver(solver, dsin = '', copy_to = None):
    """
    set_solver (solver, dsin = '')
    
    Set the solver in the simulation input file dsin
    - solver = int or anything for which int(solver) lies between (1-28)
    - dsin (default = dsin.txt): string, input file to change
    - copy_to (default = no copy of dsin.txt): string, name of new (copied) file
    
    Returns: nothing
    
    """
    
    if dsin == '':
        dsin = 'dsin.txt'
            
    try:
        if int(solver) not in range(1,29):
            raise ValueError('solver %d is not between 1 and 28' % int(solver))
    except:
        raise TypeError('solver %s could not be converted into an integer' \
                        % solver)
    solver = int(solver)                   
    
    dsin_file = open(dsin, 'r')
    dsin_adapted = open('dsin_temp.txt', 'w')
    for s in dsin_file:
        if s.find('Integration algorithm as integer') > -1:
            # get the currently used Algorithm and replace it by something else
            splitted = s.split()            
            old_solver = splitted[0]
            s = s.replace(old_solver, str(solver))
        dsin_adapted.write(s)
    dsin_file.close()
    dsin_adapted.close()
    
    # only if original dsin had to be replaced
    if copy_to != None:    
        copyfile('dsin_temp.txt', copy_to)
    
def set_ststst(start = 0, stop = 86400, step = 60, dsin = '', copy_to = None):
    """
    better name to be found.  set_times? 
    """
    pass

def set_par(parameter, value, dsin = '', copy_to = None):
    """
    anyone?
    """
    pass


def kill_after(proc, seconds):
    """
    Wait for a process to finish, else kill it after timeout (s)
    Returns True if process had to be killed
    
    Based on: http://stackoverflow.com/questions/1359383/python-run-a-process-and-kill-it-if-it-doesnt-end-within-one-hour    
    """

    start = time.time()
    
    if seconds is not None:
        # we wait or kill
        end = start + seconds
        interval = min(seconds / 1000.0, 0.25)
        while True:
            res = proc.poll()
            if res is not None:
                return False
            if time.time() >= end:
                proc.kill()
                return True
            time.sleep(interval)        
    else:
        print "Waiting for Godot"
        while True:
            res = proc.poll()
            if res is not None:
                return False
            time.sleep(1)

def analyse_log(log_file):
    """
    Check if the simulation ended successfully, which solver was used and
    how much time it took.  Optionally, show the number of this and that
    Returns a dictionary with the results
    
    """
    
    summary = {'successful':False}
    lf = open(log_file, 'r')
    lines = lf.readlines()    
    for line_number, line in enumerate(lines):
        if line.find('Integration terminated successfully at T =') > -1:
            summary['successful'] = True
        elif line.find('CPU time for integration') > -1 or \
             line.find('CPU-time for integration') > -1:
            summary['CPU time'] = line.split(' ')[-2]
        elif line.find('Number of (successful) steps') > -1:
            summary['steps_ok'] = line.split(' ')[-1]
        elif line.find('Number of rejected steps') > -1:
            summary['steps_nok'] = line.split(' ')[-1]
        elif line.find('Integration started at 0 using integration method:') > -1:
            summary['algorithm'] = lines[line_number + 1].strip('\n')
        elif line.find('Integration started at T = 0 using integration method') > -1:
            summary['algorithm'] = line.split(' ')[-1].strip('\n')
        elif line.find('This simulation timed out and was killed') > -1:
            summary['successful'] = False
            summary['timed out'] = True
        elif line.find('Corresponding result file') > -1:
            summary['result file'] = line.split(' ')[-1].strip('\n')
    lf.close()
    if summary.has_key('steps_nok'):    
        summary['perc_wrong'] = 100. * float(summary['steps_nok']) / \
                                float(summary['steps_ok'])
    else:
        summary['perc_wrong'] = 0
    return summary        
        
def run_ds(dymosim = '', dsin = '', result = ''):
    """
    run_ds(dymosim = '', file_in = '', result = '')    
    
    Call a dymosim.exe and execute it with the following arguments:
        - dymosim (default = dymosim.exe): filename of the dymosim.exe object to 
          be called
        - dsin (default = dsin.txt): input file
        - result (default = result.mat): result file
                          
    This function returns the Popen process it starts.
    
    """
    
    arguments = {'dymosim':'dymosim', 
                'dsin':'dsin.txt', 
                'result':'result.mat'}
    
    for arg in arguments:
        if arg is not '':
            arguments[arg] = eval(arg)
        
    
    oscmd = ' '.join([arguments['dymosim'], '-s', arguments['dsin'], 
                      arguments['result']])
        
        
    proc = Popen(oscmd)
    return proc    
    


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
