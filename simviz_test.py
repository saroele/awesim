# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 09:39:55 2011

@author: RBa
"""


import os
import simviz
import numpy as np
from simman import Simulation, Simdex

path = r'E:\6_Teaching\2011_RobKenneth\result'
os.chdir(path)

s1 = Simulation('HydronischHoutK14')

simviz.plt_comfort(s1)

