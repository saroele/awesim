# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 11:42:03 2011

@author: RBa
"""

import os, pp, numpy, shutil
from pymosim import *

work_dir = r'C:\Workspace\BS2011\PR_Test'

os.chdir(work_dir)

# first we must make a dictionary of the parameteric we want to do
# and we add the parameters and their values we want to run
dickie = {}
dickie.update({r'PVTimeOff_nom' : numpy.arange(100,1000,100)})
# dickie.update({r'omvormer.eff' : numpy.array([0.9,0.85])})
# now we set everything ready for parallel python, the return consistst of a 
# list of cd-paths where all sub-sets are located as for parallel python
sub_dir = set_parametric_run(work_dir, dickie)

# initialise the parallel python stuff
ppservers=()
job_server = pp.Server(ncpus = 2, ppservers = ppservers)
# define the inputs and start the simulations
inputs = tuple(sub_dir)
jobs = [(input, job_server.submit(start_parametric_run, args = (input,), modules = ("subprocess","os"))) for input in inputs]
job_server.wait()
job_server.destroy()

close_parametric_run(work_dir, sub_dir)

