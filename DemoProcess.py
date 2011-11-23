# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 11:06:19 2011

@author: RDC
"""

from simman import Simulation, Simdex, Process

mothers=['c1', 'c2']
parameters={'cap1':'c1.C', 'res':'r.R'}
sub_pars={'cap':'C'}
variables={}
sub_vars={'Qflow':'heatPort.Q_flow'}
pp=['Qflow10 = 10 * Qflow']

process = Process(mothers=mothers, parameters=parameters, 
                    sub_pars=sub_pars, variables=variables, 
                    sub_vars=sub_vars, pp=pp)

class DummyProcess(object):
    def invoke_upper(self):
        print locals()
        


class MyClass(object):
    
    def __init__(self, process):
        self.process=process
        
    def my_method(self, a, b):
        print 'return the sum of ', a, ' and ', b        
        return a+b
        

dummyprocess=DummyProcess()

myclass=MyClass(dummyprocess)

c= myclass.my_method(3, 5.67)
print c

simdex=Simdex()