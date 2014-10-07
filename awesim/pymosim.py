# -*- coding: utf-8 -*-
"""
Created on Wed Feb 09 22:05:39 2011

@author: RDC
"""
from __future__ import division
import time, os, itertools, shutil
from subprocess import Popen
import subprocess
import sys
import shutil
import copy
import numpy as np
import pdb


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
        shutil.copyfile('dsin_temp.txt', copy_to)
 
 
def set_init_from_file(source_path, destination_path, startTime=None, stopTime=None):
    """
	replace the 'destination' file with a 'source 'file. 
	Function to replace a Dymola simulation file (dsin.txt, with initial values for all variables) by 
	a Dymola simulation file (dsfinal.txt, with end values for every variable). 
		- startTime: reset the txtline indicating the start of the simulation with startTime (s)
		- stopTime: reset the txtline indicating the end of the simulation with stopTime (s)
	"""
    
    if startTime==None and stopTime==None:
        shutil.copyfile(source_path,destination_path)
    else:
        with    open(source_path,'r') as input_file, \
                open(destination_path, 'w') as output_file:

            if startTime==None:  # startTime is not None
                for line in input_file:
                    if '# StopTime' in line:
                        output_file.write( str(stopTime)   +\
                        '                       # StopTime     Time at which integration stops\n')             
                    else:
                        output_file.write(line)
                        
            if stopTime==None: # stopTime is not None
                for line in input_file:
                    if'# StartTime' in line:
                        output_file.write( str(startTime)   +\
                        '                       # StartTime     Time at which integration starts\n')
                    else:
                        output_file.write(line)
            
            else:               # both are not None
                for line in input_file:                        
                    if '# StopTime' in line:
                        output_file.write( str(stopTime)   +\
                        '                       # StopTime     Time at which integration stops\n')
                    elif'# StartTime' in line:
                        output_file.write( str(startTime)   +\
                        '                       # StartTime     Time at which integration starts\n')
                    else:
                        output_file.write(line)
    

def set_par(parameter, value, dsin='dsin.txt', copy_to=None):
    """
    Change a parametervalue in an existing dsin.txt
    
    Parameters
    ----------
    
    parameter: string with parameter name
    value: value to put for this parameter
    dsin: filename of the dsin.txt. 
    copy_to = filename of the resulting file, defaults to dsin
    
    
    """

    parameter_dsin = '# ' + str(parameter)
    orig_file = open(dsin, 'r')
    lines = orig_file.readlines() # list of strings, each ending with '\n'
    orig_file.close()
    
    for linenumber, s in enumerate(lines):
        if s.find(parameter_dsin) > -1:
            # first we check that the parameter is not an auxiliary
            # parameter (5th value of the 'array line' should be a 1)
            splitted = s.split()            
            try:
                if not splitted[-4] == '1':
                    raise ValueError("The parameter %s is of type 'auxiliary'.\n\
                    it cannot be set in the dymosim input file. " % (parameter))# check if the value to write is in this line, or the previous one
            except:
                print "The parameter %s is of type 'auxiliary'.\n\It cannot be set in the dymosim input file. " % (parameter)
                raise
                    
            # check structure of the file
            if len(splitted) != 8:
                # for some reason, the line is splitted.  We have to change the 
                # second value of the previous line
                prev_splitted = lines[linenumber-1].split()              
                old_value = copy.copy(prev_splitted[1])
                prev_splitted[1] = str(value)
                prev_splitted.append('\n')
                lines[linenumber-1] = ' '.join(prev_splitted)
            else:
                # all is nicely in one line
                old_value = copy.copy(splitted[1])
                splitted[1] = str(value)
                splitted.append('\n')  
                lines[linenumber] = ' '.join(splitted)
            # we don't need to search the rest of the file
            break            

        
    # Write the file
    
    if copy_to is None:
        copy_to = dsin
    
    writefile = file(copy_to, 'w')
    writefile.writelines(lines)
    writefile.close()
    
    print '%s found in %s: %s is replaced by %s' \
           % (parameter, dsin, old_value, value)
    

def start_parametric_run(path):
    
    import os, subprocess
    
    os.chdir(path)
    # trying to make it work on the TME server
    arguments = {'dymosim':'dymosim'}
    oscmd = ' '.join([arguments['dymosim'], '-s'])
                      
    errormessage = ''
    try:    
        proc = subprocess.Popen(oscmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = proc.communicate()
        # proc.wait()
    except:
        errormessage = 'Simulation not even started'
    
    print output
    
    return [errormessage, proc.returncode]    

def set_parametric_run(path, dickie):
    # we make a separate array of all parameter names
    parameters = dickie.keys()
    print 'The parameters are read from the dictionary,'
    # and one with all data
    iterables = dickie.values()
    # from which we can calculate all combinations
    combination = []
    for i in itertools.product(*iterables):
        combination.append(i)
    else:
        print 'and all combinations defined.'
    # for CPU time saving measures, we transofrm 'parameters' once for reading
    # in dsin.txt
    parameters_dsin = ['# ' + str(p) for p in parameters]

    # for which we can make a lot of new dsin.txt-files, copy them to a subfile
    # and make the file ready to run with a dymosim.exe and the input files.
    path_list = []    
    for i,c in enumerate(combination):
        path_run = set_simulation(path, parameters_dsin, c, i)
        path_list.append(path_run)
        print 'Setting simulation %s out of %s' %(i,len(combination))
    else:
        print 'All simulations are set.'
        
    return path_list

def set_sensitivity_run(path, dickie, deviation = 0.1):
    # we make a separate array of all parameter names
    # and one with all data
    print 'The parameters are read from the dictionary'
    parameters = dickie.keys()
    defaults = dickie.values()
    combination = [copy.copy(defaults)]
 
    for index, item in enumerate(defaults):
        if isinstance(item, int):
            if item < (1./deviation):
                lower = item-1
                upper = item+1
            else:
                lower = int(round((1-deviation) * item))
                upper = int(round((1+deviation) * item))
        else:
            lower = (1 - deviation) * item
            upper = (1 + deviation) * item
        new = copy.copy(defaults)
        new[index]=lower
        combination.append(copy.copy(new))
        new[index]=upper
        combination.append(new)
    
    # for CPU time saving measures, we transform 'parameters' once for reading
    # in dsin.txt
    parameters_dsin = ['# ' + str(p) for p in parameters]

    # for which we can make a lot of new dsin.txt-files, copy them to a subfile
    # and make the file ready to run with a dymosim.exe and the input files.
    path_list = []    
    for i in range(len(combination)):
        path_run = set_simulation(path, parameters_dsin, combination[i], i)
        path_list.append(path_run)
        print 'Setting simulation %s out of %s' %(i,len(combination))
    else:
        print 'All simulations are set.'
        
    return path_list

def cleanup_parrun(workdir, targetdir=None, subdir=None, remove=False):
    """
    Clean a folder with simulations as created by a parametric run.  
    
    Parameters
    ----------
    
    workdir = the main map where all subsets are located in
    targetdir = path to the folder that will contain all resulting files    
    subdir = a list of paths to all subsets.  If not provided, all subdirs
    containing 'run_' are treated.
    remove: boolean, if True, the subdirs are all removed.
    
    Result
    ------
    The .mat files, log files and dsin files are put in a folder results_pr
    and all subset folders (run_x) are removed.
    
    """

    #pdb.set_trace()
    if targetdir is None:
        targetdir = os.path.join(workdir, 'results_pr')
        
    os.makedirs(targetdir)
    files_copied = np.array([False, False, False])    
    
    if subdir is None:
        subdir_short = [f for f in os.listdir(workdir) if f.find('run_')>-1]
        subdir = [os.path.join(workdir, p) for p in subdir_short]
        
    for folder in subdir:
        print 'Processing ', folder
        run_id = os.path.split(folder)[-1]

        # find and copy the .mat file        
        resultfile_oldpath = os.path.join(folder, 'dsres.mat')
        if os.path.exists(resultfile_oldpath):
            resultfile_newpath = os.path.join(targetdir, run_id + '.mat')
            shutil.move(resultfile_oldpath, resultfile_newpath)
            files_copied[0] = True
        else:
            print resultfile_oldpath, ' not found.'
            
        # find and copy the log file
        logfile_oldpath = os.path.join(folder, 'dslog.txt')
        if os.path.exists(logfile_oldpath):
            logfile_newpath = os.path.join(targetdir, run_id + '.txt')
            shutil.move(logfile_oldpath, logfile_newpath)
            files_copied[1] = True
        else:
            print logfile_oldpath, ' not found.'
            
        # find and copy the dsin.txt file
        dsinfile_oldpath = os.path.join(folder, 'dsin.txt')
        if os.path.exists(dsinfile_oldpath):
            dsinfile_newpath = os.path.join(targetdir, 'dsin_' + run_id + '.txt')
            shutil.move(dsinfile_oldpath, dsinfile_newpath)
            files_copied[2] = True
        else:
            print dsinfile_oldpath, ' not found.'

        # remove the folder with al content
        if np.all(files_copied) and remove:        
            shutil.rmtree(folder)

         

def set_simulation(path, parameters, values, copy_to = None, dsin_old = '', dsin_new = '', dymosim = ''):
    """
    This function sets everything ready for starting a new simulation based on
    the given parameters. 
    
    path = the path where a dymosim.exe and dsin.txt are waiting for you
    parameters = a list with the parameters to be changed in the input values
        these are already transformed for reading in dsin.txt by adding a '# '
    values = the values the depicted parameters need to get 
    copy_to = a follow number for creating subfiles
    """
    my_wor_dir = os.getcwd()
    os.chdir(path)  
    if dsin_old == '':
        dsin = 'dsin.txt'
    else:
        dsin = dsin_old
        
    if dsin_new == '':
        dsin_new = 'dsin.txt'
        
    replace_original = False
    if dsin==dsin_new:
        replace_original=True
        dsin_new = dsin_new.strip('.txt')+'temp.txt'
        
    if dymosim == '':
        dymosim = 'dymosim.exe'
        
    my_pars = list(parameters)
    my_vals = list(values)
    # Get parameters from sim in order of appearance
    sorted_pars = []
    sorted_vals = []
    with open(dsin,'r') as old_file:
        for f in old_file:  
            #pdb.set_trace()
            if '# ' in f:
                if any(x in f.split() for x in my_pars):
                    for i, par in enumerate(my_pars):
                        if par in f.split():
                            sorted_pars.append(par)
                            sorted_vals.append(my_vals[i])
                            my_index = i
                    my_pars.pop(my_index)
                    my_vals.pop(my_index)
    if len(sorted_pars)<len(parameters):
        print 'the following parameters where not found in dsin_old: {}'.format(list(set(parameters) - set(sorted_pars)))
    # rewrite to new file and change lines of parameters in 'parameters'
    with open(dsin,'r') as old_file, open(dsin_new,'w') as new_file:
        pars_passed=0
        line_nb=0
        search_par = sorted_pars[pars_passed]
        search_val = sorted_vals[pars_passed]
        for f in old_file:  
            if line_nb==0:
                to_write = ""
            elif '# '+search_par + '\n' in f:
                to_write_tmp = prev_line.split()
                to_write = ' '.join(to_write_tmp[:1] + [str(search_val)] + to_write_tmp[2:] + ['\n'])
                pars_passed+=1
                try:
                    search_par = sorted_pars[pars_passed]
                    search_val = sorted_vals[pars_passed]
                except:
                    pass
            else:
                to_write = prev_line
            new_file.write(to_write)
            prev_line = f
            line_nb+=1
        new_file.write(prev_line)
    
    if replace_original:
        shutil.copyfile(dsin_new, dsin_old)
    
    if copy_to != None:
        dsin_to = os.getcwd() + '\\run_' + str(copy_to)
        dsin_file = os.getcwd() + '\\run_' + str(copy_to) + '\\dsin.txt'
        dymosim_file = os.getcwd() + '\\run_' + str(copy_to) + '\\dymosim.exe'
        os.makedirs(dsin_to)
        shutil.copyfile('dsin_temp.txt', dsin_file)
        shutil.copyfile(dymosim, dymosim_file)
    else:
        dsin_to = os.getcwd()
        
    os.chdir(my_wor_dir) 
    
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
    analyse_log(log_file)
    log_file = string with path to a dslog.txt file
    
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
            summary['cpu_time'] = line.split(' ')[-2]
        elif line.find('Number of (successful) steps') > -1:
            summary['successful_steps'] = line.split(' ')[-1]
        elif line.find('Number of (model) time events') > -1:
            summary['time_events_model'] = line.split(' ')[-1]
        elif line.find('Number of (U) time events') > -1:
            summary['time_events_U'] = line.split(' ')[-1]
        elif line.find('Number of state    events') > -1:
            summary['state_events'] = line.split(' ')[-1]
        elif line.find('Number of step     events') > -1:
            summary['step_events'] = line.split(' ')[-1]
        elif line.find('Minimum integration stepsize') > -1:
            summary['step_size_min'] = line.split(' ')[-1]
        elif line.find('Maximum integration stepsize') > -1:
            summary['step_size_max'] = line.split(' ')[-1]
        elif line.find('Maximum integration order') > -1:
            summary['int_order_max'] = line.split(' ')[-1]
        elif line.find('Number of rejected steps') > -1:
            summary['steps_nok'] = line.split(' ')[-1]
        elif line.find('Integration started at 0 using integration method:') > -1:
            summary['algorithm'] = lines[line_number + 1].strip('\n')
        elif line.find('Integration started at T = 0 using integration method') > -1:
            summary['algorithm'] = line.split(' ')[-1].strip('\n')
        elif line.find('This simulation timed out and was killed') > -1:
            summary['successful'] = False
            summary['timed_out'] = True
        elif line.find('Corresponding result file') > -1:
            summary['result file'] = line.split(' ')[-1].strip('\n')
    lf.close()
    if summary.has_key('steps_nok'):    
        summary['perc_wrong'] = 100. * float(summary['steps_nok']) / \
                                float(summary['steps_ok'])
    else:
        summary['perc_wrong'] = 0
    return summary        
        
def run_ds(dymosim = '', dsin = '', result = '', dsu=''):
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
                'result':'result.mat',
                'dsu': 'dsu.txt'}
    
    for arg in arguments:
        if eval(arg) is not '':
            arguments[arg] = eval(arg)
    print arguments    
    
    oscmd = ' '.join([arguments['dymosim'], '-s','-u', arguments['dsu'],
                    arguments['dsin'], arguments['result']])

    print oscmd
    
    proc = Popen(oscmd)
    return proc    
    

def create_input_file(data, filename, discrete=False, compress=True):
    """
    Create an input file for the TimeTables from the MSL.
    The input files are in ascii format.
    
    Parameters
    ==========    
    * data: array with time as first column.  All columns of this array
      will be written in the ascii file.
    * filename: filename (with extension).  If this file already exists, it 
      will be overwritten
    * discrete: if True, the data array will be modified to become a discrete
      profile. At each timestep, an additional line will be created.  See
      the documentation of the Modelica.Timetable.mo model for more info.
    * compress: if True, all reduntant lines will be removed from data.
    """
    
    if compress:
        # rows with only zeros are removed, UNLESS they come after a row
        # containing any value. 
        
        row_contains_data = data[:,1:].any(axis=1)
        row_contains_data_rolled = np.roll(row_contains_data, 1)
        row_contains_data[0] = True
        row_contains_data[-1] = True
        keep = np.logical_or(row_contains_data, row_contains_data_rolled)
        data = data[keep]
    
    l,w = data.shape

    if discrete:
        # create a second, shifted data array.  We'll write each row of this
        # shifted array after each row of the original one.
        data_shifted = data.copy()[:-1,:]
        data_shifted[:,0] = data[1:,0]
        shape_string = '(' + str(2*l-1) + ',' + str(w) + ')'
    else:
        shape_string = '(' + str(l) + ',' + str(w) + ')'
        
    
    
    f = open(filename, 'w')
    f.write(u'#1\n')
    f.write(''.join([u'double data', shape_string,  
                     u'# Profiles created by python script: ', sys.argv[0], 
                    '\n']))
    for i in range(data.shape[0] - 1):
        f.write('\t'.join([str(v) for v in data[i,:]]))
        f.write('\n')        
        if discrete:
            f.write('\t'.join([str(v) for v in data_shifted[i,:]]))
            f.write('\n')
            
    f.write('\t'.join([str(v) for v in data[-1,:]]))
    f.write('\n')
    f.close()

    


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
