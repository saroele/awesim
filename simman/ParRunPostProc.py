# -*- coding: utf-8 -*-
"""
Created on Fri May 06 10:46:27 2011

@author: RDC
"""

from __future__ import division
import os
from shutil import copyfile
import matplotlib.pyplot as plt
import numpy as np

from simman import Simulation, Simdex
from pymosim import *


# Script settings #############################################################

analyse_logs = True
get_results = True
make_graphs = True

vars_to_get = [u'pq[32].vi.i.re', u'dwellings[17].pq.p']



###############################################################################

if analyse_logs:
    # make a list of logfiles    
    fileList = os.listdir(os.getcwd())
    log_files = [f for f in fileList if f.find('dslog_run_') > -1]
    mat_files = [f for f in fileList if f.find('result_run_') > -1]
    
    
    logs = []
    for l in log_files:
        logs.append(analyse_log(l))