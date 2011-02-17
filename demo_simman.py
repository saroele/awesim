# -*- coding: utf-8 -*-
"""

Demonstration of the Simdex and Simulation classes from
simman.py

Created on Thu Feb 17 11:36:47 2011

@author: RDC
"""

import os
from simman import Simulation, Simdex, load_simdex


# first, test the package
runfile(r'C:\Workspace\Python\SimulationManagement\test_simman.py', 
        wdir=r'C:\Workspace\Python\SimulationManagement')


# there are different ways to create a simdex

s1 = Simdex() # an empty simdex object
s1.scan() # scan current folder and add all found simulation files to the simdex

s2 = Simdex('SubfolderWithCrappyFiles')
    # create a simdex directly from a folder
    
print s1
print s2

# now, let's look at some attributes
s1.get_filenames() # method, so add ()
s1.parameters # attribute, so no ()
s1.parametermap
s1.variables
s1.variablemap

# remove a simulation from the simdex manually
s3 = s1.remove(0)
print s3
    # important: this method and some other methods return a NEW simdex
    # if you do not want this, you can do the following
s1 = s1.remove(0)
print s1
s1.parametermap

s1.exist('c2')  # returns a list with 2 lists by default: parameters matching 
                # the regex, and variables matching the regex. 
s1.get_parameters(u'c1.C')

# filter options ##############################################################
# based on identity
s3 = s1.get_identical(2)
    # returns a new simdex with all simulations that have the same set of 
    # parameters (even if they have different values)

# based on parametervalues
fltr = {'c1.C':800}
s4 = s1.filter(fltr)

# add a paramter to the filter
fltr2 = {'newParameter' : ''}
s4 = s1.filter(fltr2)

# plotting
s1.plot('c1.T')
s1.scatterplot('c1.T', 'c2.T')

# saving a specific simdex
s2 = s1.filter(fltr)
s2.save('simdex2')

del s2
print s2

s = load_simdex('simdex2')
s.plot('r.heatPort_b.Q_flow')

