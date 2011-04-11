# -*- coding: utf-8 -*-
"""
Created on Wed Feb 09 22:05:39 2011

@author: RDC
"""
from __future__ import division
import time, os, itertools, shutil
from subprocess import Popen
import subprocess
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

def set_par(parameter, value, dsin='', copy_to = None):
    """
    This function depicts the change of a parametervalue in an existing dsin.txt
    """

    if dsin == '':
        dsin = 'dsin.txt'

    parameter_dsin = '# ' + str(parameter)
    dsin_file = open(dsin, 'r+')
    for s in dsin_file:
        if s.find(parameter) > -1:
            print 'The parameter', parameter, 'is found in', dsin
        if s.find(parameter_dsin) > -1:
            splitted = s.split()
            s = s.replace(splitted[1], value)
            print 'and is replace by', value, '.'
    dsin_file.close()
    
    if copy_to != None:
        copyfile(dsin, copy_to)

def start_parametric_run(path):
    
    import os, subprocess
    
    os.chdir(path)
    arguments = {'dymosim':'dymosim', 'dsin':'dsin.txt', 'result':'result.mat'}
    oscmd = ' '.join([arguments['dymosim'], '-s', arguments['dsin'], arguments['result']])
                      
    errormessage = ''
    try:    
        proc = subprocess.Popen(oscmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = proc.communicate()
        proc.wait()
    except:
        errormessage = 'Simulation not even started'
    
    print output
    
    return [errormessage, proc.returncode]    

def set_parametric_run(path, dickie):
    # we make a separte array of all parameter names
    parameters = dickie.keys()
    print 'The parameters are read from the dictionary,'
    # and one with all data
    iterables = []
    for i in range(len(dickie)):
        iterables.append(dickie[parameters[i]])
    else:
        print 'all iterables have been depicted'
    # from which we can calculate all combinations
    combination = []
    for i in itertools.product(*iterables):
        combination.append(i)
    else:
        print 'and all combinations defined.'
    # for CPU time saving measures, we transofrm 'parameters' once for reading
    # in dsin.txt
    parameters_dsin = []
    for i in range(len(parameters)):
        parameters_dsin.append('# ' + str(parameters[i]))
    else:
        print 'Parameter names are converted for reading in dsin.txt'
    # for which we can make a lot of new dsin.txt-files, copy them to a subfile
    # and make the file ready to run with a dymosim.exe and the input files.
    path_list = []    
    for i in range(len(combination)):
        path_run = set_simulation(path, parameters_dsin, combination[i], i)
        path_list.append(path_run)
        print 'Setting simulation %s out of %s' %(i,len(combination))
    else:
        print 'Mission accomplished.'
        
    return path_list

def close_parametric_run(workdir, subdir):
    """
    This function finishes a parametric run by reordening all subsets
    
    workdir = the main map where all subsets are located in
    subdir = a list of paths to all subsets
    """

    result_path = workdir + '\\results_pr'
    os.makedirs(result_path)

    for i in range(len(subdir)):
        resultfile_oldpath = subdir[i] + '\\result.mat'
        resultfile_newpath = result_path + '\\result_run_' + str(i) + '.mat'
        existing = os.access(resultfile_oldpath, os.F_OK)        
        if  existing:
            copyfile(resultfile_oldpath,resultfile_newpath)
        else:
            print 'Run_' + str(i) + ' failed to simulate.'

        logfile_oldpath = subdir[i] + '\\dslog.txt'
        logfile_newpath = result_path + '\\dslog_run_' + str(i) + '.txt'
        existing = os.access(resultfile_oldpath, os.F_OK)        
        if  existing:
            copyfile(logfile_oldpath, logfile_newpath)
        else:
            print 'Run_' + str(i) + ' did not even start ??!!'

    for i in range(len(subdir)):
        existing = os.access(subdir[i], os.R_OK)        
        if  existing:
            shutil.rmtree(subdir[i])

def set_simulation(path, parameters, values, copy_to = None, dsin = '', dymosim = ''):
    """
    This function sets everything ready for starting a new simulation based on
    the given parameters. 
    
    path = the path where a dymosim.exe and dsin.txt
    parameters = a list with the parameters to be changed in the input values
        these are already transformed for reading in dsin.txt by adding a '# '
    values = the values the depicted parameters need to get 
    copy_to = a follow number for creating subfiles
    """
    
    os.chdir(path)    
    
    if dsin == '':
        dsin = 'dsin.txt'
    
    if dymosim == '':
        dymosim = 'dymosim.exe'

    inputpath = path + '\\inputs_pr'
    inputfiles = os.listdir(inputpath)

    dsin_file = open(dsin, 'r')
    file_data = dsin_file.readlines()
    dsin_file.close()

    for i in range(len(file_data)):
        for j in range(len(parameters)):
            if file_data[i].find(parameters[j]) > -1:
                splitted = file_data[i].split()
                splitted[1] = str(values[j])
                splitted.append('\n')
                file_data[i] = ' '.join(splitted)
                print 'The parameter', parameters[j], 'is found in', dsin, 'and is replace by', values[j], '.'

    dsin_temp = open('dsin_temp.txt', 'w')
    dsin_temp.writelines(file_data)
    dsin_temp.close()
    
    if copy_to != None:
        dsin_to = os.getcwd() + '\\run_' + str(copy_to)
        dsin_file = os.getcwd() + '\\run_' + str(copy_to) + '\\dsin.txt'
        dymosim_file = os.getcwd() + '\\run_' + str(copy_to) + '\\dymosim.exe'
        os.makedirs(dsin_to)
        copyfile('dsin_temp.txt', dsin_file)
        copyfile(dymosim, dymosim_file)
        for i in range(len(inputfiles)):
            old_path = inputpath + '\\' + inputfiles[i]
            new_path = os.getcwd() + '\\run_' + str(copy_to) +'\\'+ inputfiles[i]
            copyfile(old_path,new_path)

    return dsin_to
    

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
        if eval(arg) is not '':
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
