# -*- coding: utf-8 -*-
"""

Demonstration of the Simdex and Simulation classes from
simman.py

Created on Thu Feb 17 11:36:47 2011

@author: RDC
"""

import numpy as np
import os
from simman import Simulation, Simdex, Process, load_simdex

# first, test the package
#runfile(r'C:\Workspace\Python\SimulationManagement\test_simman.py', 
#        wdir=r'C:\Workspace\Python\SimulationManagement')

# a Simulation is a python object for 1 simulation result file

sim=Simulation('LinkedCapacities_A') #with our without .mat extension
sim.separate() # optional, makes attributes for parameters and variables

sim.parameters
sim.variables

# in big simulation files, it's not always easy to find the right parameter 
# or variable
sim.exist('q_flow')

time = sim.get_value('Time')
Q = sim.get_value(u'r.heatPort_a.Q_flow')
Q_sum = np.trapz(Q, time)
print "The total energy that flowed through the resistance is %.1f J" %Q_sum

for p,v in zip(sim.parameters, sim.parametervalues):
    print ''.join([p, ' = ', str(v)])


# there are different ways to create a simdex

s1 = Simdex() # an empty simdex object
s1.scan() # scan current folder and add all found simulation files to the simdex

s2 = Simdex('SubfolderWithCrappyFiles')
    # create a simdex directly from a folder
    
print s1
print s2

# now, let's look at some attributes
s1.get_filenames() # method, so add ()
s1.get_filenames('path')
s1.parameters # attribute, so no ()
s1.parametermap
s1.parametervalues
s1.variables
s1.variablemap
s1.h5_path

# Now look at the h5 file (open it with ViTables)
# Demonstrate natural naming of pytables
s1.openh5()
s1.h5.root.SID0001.r_dot_heatPort_a_dot_Q_flow
s1.h5.root.SID0001.r_dot_heatPort_a_dot_Q_flow.read()
s1.h5.close()

# Next step: using a process

p1 = Process(variables={'T2':'c2.T', 'dt2':'c[2].der(T)'})
s1 = Simdex(os.getcwd(), process=p1)


mothers=['c1', 'c2']
parameters={'cap1':'c1.C', 'res':'r.R'}
sub_pars={'cap':'C'}
variables={}
sub_vars={'Qflow':'heatPort.Q_flow'}
pp=['Qflow10 = 10 * Qflow']
J2kWh = 1e-6/3.6
vars_to_integrate = {'Qflow':J2kWh}

p = Process(mothers=mothers, parameters=parameters, 
                    sub_pars=sub_pars, variables=variables, 
                    sub_vars=sub_vars, pp=pp, integrate=vars_to_integrate)


# remove a simulation from the simdex manually
s3 = s1.filter_remove(['SID0000'])
print s3
    # important: this method and some other methods return a NEW simdex
    # if you do not want this, you can do the following
s1 = s1.filter_remove(['SID0000'])
print s1
s1.parametermap

s1.exist('c2')  # returns a list with 2 lists by default: parameters matching 
                # the regex, and variables matching the regex. 
s1.get(u'c1.C')

# filter options ##############################################################
# based on identity
s3 = s1.filter_similar('SID0002')
    # returns a new simdex with all simulations that have the same set of 
    # parameters (even if they have different values)

# based on parametervalues
fltr = {'c1.C':800}
s4 = s1.filter(fltr)

# add a paramter to the filter
fltr2 = {'newParameter' : ''}
s4 = s1.filter(fltr2)

# plotting
s3.plot('c1.T')
s3.scatterplot('c1.T', 'c2.T')

# saving a specific simdex
s2 = s1.filter(fltr)
s2.save('simdex2')



s = load_simdex('simdex2')
s.plot('r.heatPort_b.Q_flow')