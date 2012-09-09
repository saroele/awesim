# -*- coding: utf-8 -*-
"""


@author: Roel De Coninck
"""
#import numpy as np
#import os
#import scipy.io
#import re
import copy
#import matplotlib.pyplot as plt
#from matplotlib.dates import date2num
#import cPickle as pickle
#import bisect
#import tables as tbl
#from datetime import datetime, timedelta
#import pandas
import pdb

class Process(object):
    """
    Class defining pre- and post processing of a simulation
    
    The variable 'Time' is automatically created
    
    This documentation needs completion, mainly the __init__() method
    
    
    """
    
    def __init__(self, mothers=None, parameters=None, sub_pars=None, variables=None,
                 sub_vars=None, pp=None, integrate=None):
        """Instantiate the Process object
        
        parameters
        ----------
        - mothers: a list of mothers objects, as strings.  Arrays are supported
          in that case the '[' will be converted into '_' and the ] will be 
          removed.  Example d[1] becomes d_1.  In combination with a sub_par
          or sub_var, the . will be replaced by a '_'.  
          So the shortname for d[12].TSto will be d_12_TSto
        - parameters: a dictionary with shortname:longname pairs for parameters.
          The shortnames can be used in postprocessing formulas (see pp)
        - sub_pars: a dictionary with shortname:longname pairs for parameters 
          that are children of the mothers.  
          Example: if you have a c1.TInitialValue and a c2.TInitialValue, 
          you can have mothers = ['c1', c2'] and 
          sub_vars = {'TStart': 'TInitialValue'}
        - variables: analoguous to parameters
        - sub_vars: analoguous to sub_pars
        - pp: a list of post-processing strings.  Each string will be executed
          and will generate a new variable that will be stored and that can be
          used in following post porcessing strings. If the variable names used
          in a pp string are sub_vars or sub_pars, the string will be executed
          for each mother.  
        
          It is crucial to surround the names of variables by spaces so they can be
          looked up in the parameters and variables dictionaries.

        - integrate: dictionary with variables to be integrated over time.  
          The dictionary has shortname:scaling pairs, the scaling factor will
          be multiplied to the result of np.trapz(shortname, time).  The 
          resulting value will be stored under shortname_Int.  If the shortname
          is a sub_var, the integration will be done for every mother.
          
          Example: mothers = ['c1', c2'], sub_vars = {'Q':'Q_flow'}, 
          integrate = {'Q':1e-6} will create c1_Q_Int and c2_Q_Int.
          
        """
        #pdb.set_trace()
        pp_int = []
               
        # make/complete the variables and parameter dicts, full paths
        if variables is None:
            self.variables = {}
        else:
            self.variables = copy.copy(variables)
               
        if parameters is None:
            self.parameters = {} 
        else:
            self.parameters = copy.copy(parameters)
               
        # make sure there are sub_vars and sub_pars, or create empty ones
        # otherwise, an exception will be thrown in simulation.postprocess.convert(p)
        if sub_vars is None:
            self.sub_vars = {}
        if sub_pars is None:
            self.sub_pars = {}
        if mothers is not None:
            self.mothers = copy.copy(mothers)
            for m in self.mothers:
                m_orig = copy.copy(m)
                m = m.replace('[', '_')
                m = m.replace(']', '')
                if sub_vars is not None:
                    self.sub_vars = copy.copy(sub_vars)                    
                    for shortname, longname in self.sub_vars.iteritems():
                        self.variables['_'.join([m, shortname])] = '.'.join([m_orig, longname])
                else:
                    self.sub_vars = {}
                if sub_pars is not None:
                    self.sub_pars = copy.copy(sub_pars)
                    for shortname, longname in self.sub_pars.iteritems():
                        self.parameters['_'.join([m, shortname])] = '.'.join([m_orig, longname])
                else:
                    self.sub_pars = {}
            self.mothers = [m.replace('[', '_') for m in self.mothers]
            self.mothers = [m.replace(']', '') for m in self.mothers]
        else:
            self.mothers = []
        # paramaeters = contains short and long names of the parameters we need.
        # The data will be extracted from the .mat files and put in a dictionary
        # The short names will be the keys in that dictionary     
        if self.variables is not {}:
            if not self.variables.has_key('Time'):
                self.variables['Time'] = 'Time'
        
        if integrate is not None:
            for name, conversion in integrate.iteritems():
                # maka a pp string for this integration action
                s = ''.join([name+'_Int',
                              ' = ',
                              'np.trapz( ',
                              name,
                              ' , ',
                              'Time',
                              ' ,axis=0)*',
                              str(conversion),
                              ' if ',
                              name,
                              ' .shape[0]== Time .shape[0] else np.array([0.0])'
                              ])
                pp_int.append(s)

        # now put the pp_int in front of the self.pp if any        
        self.pp = []
        self.pp.extend(pp_int)
        if pp is not None:
            self.pp.extend(pp)
            
        
        
    def __str__(self):
        """Return a print string"""
        
        s = '\n'+79*'-'+'\n'
        s += 'The content of this Process object is:\n'
        s += 'Parameters:'
        for shortname,longname in sorted(self.parameters.iteritems()):
            s += '\t'.join(['\n', shortname ,'=', longname])
        
        s += '\n\nVariables:'
        for shortname,longname in sorted(self.variables.iteritems()):
            s += '\t'.join(['\n', shortname ,'=', longname])
        
        s += '\n\nPost-processing:\n'        
        if self.pp is not None:
            ppprint = '\n'.join(self.pp)
        else:
            ppprint = '\nNo post-processing defined\n'
        s += ppprint
        
        return s
        