# -*- coding: utf-8 -*-
"""
class Simdex:
-------------
    A Simdex object is an index of all the parameters and variables in several 
    .mat files.  It can be used to select a set of .mat files from a large set
    in order to analyse the results of this set of simulations.
    
    Most important attributes are:
        - self.simulations (a list of filenames of the indexed simulations)
        - self.parameters (a list of all parameters in the set of simulations)
        - self.parametermap (numpy array mapping which simulation has 
          which parameters)        
        - self.parametervalues (numpy array containing the values of the 
          mapped parameters)
        - self.variables (a list of all variables in the set of simulations)
        - self.variablemap (numpy array mapping which simulation has 
          which variables)
        - self.filterset (dictionary tracking all executed filter options on 
          the current set) 
          
        There is no attribute containing the variable values themselves because 
        that would blow up the size of simdex.  These values are kept in the 
        original .mat files and only extracted when needed

Most important methods (* = implemented):

    * __init__(folder): create simdex based from all .mat simulation files in 
      the current folder
    - update(folderlist): update simdex with all .mat files found in a list of 
      folders. Checks if the found files are already in simulations, else they 
      are added and their parameters and ALL attributes are updated with the 
      info found in the new files
    * remove(sim_id): remove a simulation from the simdex
    * print(): gives a nice overview of the indexed simulations 
    * filter(dictionary): this method takes as input a dictionary with parameter
      name/value pairs.   It returns a new Simdex object with those simulations
      that have exactly the same values for those parameters.  
      If a parameter is asked with '*' as value, all simulations that HAVE 
      the parameter are selected.
    * getidentical(sim_id):  this method takes as input a simulation number and 
      returns a new simdex object all simulations with identical parameter set 
      (so all variants of the model with only changed parameter values).  
    * exist(re): similar to Simulation.exist method. Give a search string as re 
      (regular expression) and you get a list of all parameters and variables 
      that satisfy the re. 
    * plot(var): directly create a plot showing the given var for each of the 
      simulations in self.
    - get_values(var or par): get an array with the values of the variable or
      parameter for each of the simulations in the simdex
    * get_parameter(par): to be merged in get_values!!
    * save(filename): saves the simdex by pickling (cPickle) to filename.  The 
      simdex can be loaded later on with the function load_simdex(filename). 
      This is no method of the class Simdex, so to be imported separately from
      this module

@author: Roel De Coninck
"""

import numpy as np
import os
#import scipy.io
import re
import copy
import matplotlib.pyplot as plt
#from matplotlib.dates import date2num
import cPickle as pickle
import bisect
import tables as tbl
#from datetime import datetime, timedelta
#import pandas
import pdb
from .simulation import Simulation
from .result import Result
from .pymosim import analyse_log


class Simdex:
    """
    A Simdex object is an index of all the parameters, variables and 
    post-processing results in several .mat files.  
    It is used to treat different result files with a single or different 
    Process objects and keep the full processed results.
    Methods for easy filtering of the set of simulations are provided.
    
    The variables are stored in an h5 file.  Parameters are only stored locally
    in the attributes of the simdex object. When a filter is applied to a
    simdex, the resulting subset of simulations uses the same h5 file, so 
    there is no unnecessary copying of big files on the hard drive.
        
    Overview of most important attributes :
        - self.simulations: a list with unique ID's (SID) of the simulations
        - self.parameters: a list of all parameters in the set of simulations.
          The parameters are listed here with their FULL names
        - self.parametermap: numpy array mapping which simulation has 
          which parameters     
        - self.parametervalues: numpy array containing the values of the 
          mapped parameters
        - self.variables: a list of all variables in the set of simulations. 
          The variables are listed here with their FULL names.
        - self.variablemap: numpy array mapping which simulation has 
          which variables
        - self.filterset: dictionary tracking all executed filter options on 
          the current set
        - self.vardic (optional): a mapping of shortname:longname pairs for 
          variables
        - self.pardic (optional) a mapping of shortname:longname pairs for 
          parameters
        - self.simulations: list of the SID
        - self.files: dictionary with SID:filename pairs
        - self.identifiers: optional dictionary with short meaningful identifiers
          for the simulations.  Is used as legend in plots by default (unless 
          empty).
          
        
    Most important methods (* = implemented):

    * __init__(folder): create simdex based from all .mat simulation files 
      in the current folder
    - update(folderlist): update simdex with all .mat files found in a list of 
      folders. Checks if the found files are already in simulations, else they 
      are added and their parameters and ALL attributes are updated with the 
      info found in the new files
    * remove(sim_id): remove a simulation from the simdex
    * print(): gives a nice overview of the indexed simulations 
    * filter(dictionary): this method takes as input a dictionary with parameter 
      name/value pairs.   It returns a new Simdex object with those simulations
      that have exactly the same values for those parameters.  
      If a parameter is asked with '*' as value, all simulations that HAVE 
      the parameter are selected.
    * getidentical(sim_id):  this method takes as input a simulation number and 
      returns a new simdex object all simulations with identical parameter set 
      (so all variants of the model with only changed parameter values).  
    * exist(re): similar to Simulation.exist method. Give a search string as re 
      (regular expression) and you get a list of all parameters and variables 
      that satisfy the re. 
    * plot(var): directly create a plot showing the given var for each of the 
      simulations in self.
    - get_values(var or par): get an array with the values of the variable or
      parameter for each of the simulations in the simdex
    * get_parameter(par): to be merged in get_values!!
    - save(filename): saves the simdex by pickling (cPickle) to filename.  The 
      simdex can be loaded later on with the function load_simdex(filename)
      
    FEATURES:
    - all found .mat files are tested.  If they to not have the structure of
      a dymola simulation result, they are not indexed

    Important CONVENTIONS:
        - we make a distinction between parameters (one single value)
          and variables (a timeseries of values)
        - simulation identity (SID) is a string with format SIDxxxx
            
    Possible IMPROVEMENTS:
    - multi folder file search
    - support multiple filter (well, it's possible, but the filtersets will 
      become a bit messed up (last filterset overwrites existing filterset if 
      it has the same keys)

    """
    
    def __init__(self, folder='', h5='simdex.h5', process=None, verbose = False):
        '''
        Create a Simdex object.  
        
        Folder is a single folder to be sought for .mat files
        If folder = '', no work directory is indexed.
        
        h5 is the hdf5 file to which this simdex will be linked.
        '''
        # First, initialise some  attributes
        if verbose == True:
            self.verbose = True
        else:
            self.verbose = False
        
        # List with sim ID's, fixed id's for simulation filenames
        # Format: SID1, SID2, SID3, ...        
        self.simulations = []
        # dictionary with SIDx:path pairs (path are full pathnames)        
        self.files = {}
        self.identifiers = {}
        # used for plotting
        self.year = 2010
        self.time4plots = {}
        
        self.process = process
        
        # The pytables file to which this simdex is linked
        # Will be created in the current work directory, check first if it exists
        self.h5_path = os.path.join(os.getcwd(), h5)
        overwrite = 'y'        
        if os.path.exists(self.h5_path):
            print 'This file already exists: %s' % self.h5_path
            overwrite = raw_input('Overwrite (Y/N)? : \n')
            print '\n'
        
        if overwrite == 'y' or overwrite == 'Y':
            self.h5 = tbl.openFile(self.h5_path, 'w', title='Simdex file')
            self.h5.close()
        else:
            raise NotImplementedError("Remove the file first !")
        # dictionary with the filters previously applied on this simdex
        self.filterset = dict()

        if folder == '' :
            # an empty simdex is created
            self.parameters = []
            self.variables = []
            self.parametermap = np.ndarray((0, 1))
            self.parametervalues = copy.copy(self.parametermap)
            self.variablemap = np.ndarray((0, 1))
            
        
        # here we get a list with all files in 'folder' that end with .mat
        
        elif os.path.exists(folder):
            self.scan(folder, process=process)
        else:
            raise IOError('folder does not exist')
        
        if self.h5.isopen:
            self.h5.close()

    def openh5(self):
        """Open the h5 file in append mode"""
        if not self.h5.isopen:
            self.h5 = tbl.openFile(self.h5_path, 'a')
        
    
    def _last_key(self):
        """Return the last key from the pytables file, or '' if empty file"""
        
        self.openh5()
        
        try:        
            meta = self.h5.getNode(self.h5.root.Metadata)
            return meta.cols.SID[-1]
        except(tbl.NoSuchNodeError):
            return ''
            

    def _gen_key(self):
        '''Generate a new key, based on last key'''
        
        # get last key.  
        last_key = self._last_key()
        
        if last_key == '':
            return 'SID' + str(format(0, '04d'))
        else:            
            last_int = int(last_key.split('SID')[-1])
            return 'SID' + str(format(last_int+1, '04d'))

        
    def __str__(self):
        '''
        Prints the Simdex object as a list of the indexed simulations
        with their simID's
        '''
        
        s= '\nSID     Filename\n'
        for k in self.simulations:
            s = ''.join([s, k, ' ', self.files[k], '\n'])
            
        return s
                    
    def scan(self, folder='', process=None, timecheck=True):
        """
        Scan a folder for .mat files and add them to the simdex
        
        Parameters
        ----------
        - folder: a single folder to be sought for .mat files
          If folder == '', the current work directory is indexed
        - process: a post-processing to be applied to each mat file
        - timecheck: if True, verify that all indexed simulations have the 
          same start and stop times.  
        
        """
        
        #pdb.set_trace()
        
        if process is None:
            process = self.process
        
        if folder == '' :
            folder = os.getcwd()
        
        try:
            filenames = self.__get_files(folder, '.mat')
        except IOError:
            raise IOError('folder %s does not exist' % (folder))
        
        if len(filenames) == 0:
            raise ValueError("No .mat files found in %s" % (folder))
        
        # Convert to full path filenames to avoid confusion
        full_path_filenames = []
        for i in range(len(filenames)):
            full_path_filenames.append(os.path.join(folder,filenames[i]))
        
        # run the following loop only when this is the first time files are
        # being indexed
        if self.simulations == []:
            # index is the pointer to the current file in full_path_filenames
            index = -1
            ########################################################################
            # Now we take the first .mat file and use this as a basis for our index
            first_file_indexed = False
            while first_file_indexed == False:
                index += 1
                simulation_file = False
                try:
                    # We try the .mat files one by one until 
                    # we find a first Dymola file
                    sim = Simulation(full_path_filenames[index])
                    simulation_file = True
                except MemoryError:
                    print 'WARNING: %s could not be indexed because of a MemoryError.\nThe file is probably too big.  It could help to try in a fresh python instance' % (full_path_filenames[index])
                except:
                    print '%s is no Dymola file.  It is not indexed' % \
                        (full_path_filenames[index])
    
                if simulation_file:                
                    # Now, check the simulation runtime and confirm with user that it is 
                    # correct.  For the next simulation files, the runtime will be compared
                    # to this one to decide if the file is ok or not. 
                    time = sim.get_value('Time')
                    
                    print 'The first found simulation, %s, runs from %d s till %d s' % \
                        (sim.filename, time[0],time[-1])
                    timeOK = raw_input('Is this correct? y/n : ')
                    print '\n'
                    if timeOK == 'y' or timeOK == 'Y':
                        self.simulationstart = time[0]
                        self.simulationstop = time[-1]
                        first_file_indexed = True
                    else:
                        print '%s is NOT indexed' % (sim.filename)
                

            # The first simulation file is indexed and the attributes are 
            # initialised.  
            self.index_one_sim(sim, process=process)
            print '%s indexed' % (sim.filename)
            ########################################################################
            # The next step is to index all remaining files
            
            index += 1
            while index < len(full_path_filenames):
                # We try to index the remaining .mat files one by one 
                try:
                    sim = Simulation(full_path_filenames[index])
                except :
                    pass
                else:
                    if timecheck:                    
                        # Now, check the simulation runtime against previously 
                        # confirmed start and stop times
                        time = sim.get_value('Time')
                        
                        if self.simulationstart == time[0] and \
                            self.simulationstop == time[-1]:
                            # index this new simulation 
                            self.index_one_sim(sim, process=process)
                            print '%s indexed' % (sim.filename)
                                                
                        else:
                            print '%s, runs from %d s till %d s, therefore, it \
                               is NOT indexed' % (sim.filename, time[0],time[-1])
                    else:
                        # index this new simulation 
                        self.index_one_sim(sim, process=process)
                        print '%s indexed' % (sim.filename)
                index += 1
        
        else:
            # there are already files in the simdex
            index = 0
            while index < len(full_path_filenames):
                # We try to index the remaining .mat files one by one 
                try:
                    sim = Simulation(full_path_filenames[index])
                except :
                    pass
                else:
                    # Now, check the simulation runtime against previously confirmed
                    # start and stop times
                    time = sim.get_value('Time')
                    
                    if self.simulationstart == time[0] and \
                        self.simulationstop == time[-1]:
                        # index this new simulation 
                        self.index_one_sim(sim, process=process)
                        print '%s indexed' % (sim.filename)
                            
                    else:
                        print '%s, runs from %d s till %d s, therefore, it is NOT \
                             indexed' % (sim.filename, time[0],time[-1])
                
                index += 1            
        
        self.h5.close()
                


    def get_filenames(self, form='filename'):
        """
        Return a list of the filenames
        
        form = 'filename' (default): only the filenames
        form = 'path' : full path name
        """
        
        if form == 'path':
            result = self.files.values()
        elif form == 'filename':
            result = [os.path.split(x)[1] for x in self.files.values()]
        else:
            print 'form is not recognised'
            raise ValueError
            
        return result
    
    
    
    def exist(self, regex, tp = 'all'): 
        '''
        exist(regex, tp='all') 
        
        regex = regular expression
        tp = 'all' (default), 'par' or 'var'
        
        This function checks if a variable name exists in the index.
        In all cases, it returns a list.
        If tp = 'all' this list contains 2 lists of strings: 
            the first with all parameters that satisfy the regex, 
            the second with variables that satisfy regex
        If tp = 'par' or 'var' the return list only contains 
        the corresponding list.
        
        Attention: if you want to check if eg.  c[3].T exists, you have to 
        escape the [ and ] with a backslash, like this:
        self.exist('c\[3\].T). Otherwise c3.T is sought for. This is 
        because in regex syntax, [] is used to indicate a set of characters.
        '''
        
        
        p = re.compile(regex, re.IGNORECASE)
        if tp == 'all' or tp == 'par':
            # we search for parameters in the fullnames
            matchespar = []
            for par in self.parameters:
                m = p.search(par)
                if m:
                    matchespar.append(par)
            # we also search in the shortnames, if present
            if self.__dict__.has_key('pardic'):            
                for par in self.pardic:
                    m = p.search(par)
                    if m:
                        matchespar.append(par)
            
            
        if tp == 'all' or tp == 'var':
            # we search for variables in the fullnames
            matchesvar = []
            for var in self.variables:
                m = p.search(var)
                if m:
                    matchesvar.append(var)
                    
            # we also search in the shortnames
            if self.__dict__.has_key('vardic'):
                for var in self.vardic:
                    m = p.search(var)
                    if m:
                        matchesvar.append(var)                    
          
        if tp == 'all':
            result = [matchespar, matchesvar]
        elif tp == 'par':
            result = matchespar
        elif tp == 'var':
            result = matchesvar
        else:
            print 'wrong input for tp'
            raise ValueError

        return result
    
    
    def __get_files(self, directory, non_wildcard_pattern):
        '''
        This function returns a list of filenames as strings, satisfying the 
        nonWildCardPattern
        '''
        
        # This function was found here on 26/11/2010:
        # http://codecomments.wordpress.com/2008/07/10/find-files-in-directory-using-python/
    
        fileList = os.listdir(directory)
        return [f for f in fileList if f.find(non_wildcard_pattern) > -1]



        
    def index_one_sim(self, simulation, process=None):
        '''
        Add a Simulation instanct to a Simdex instance
        
        This method indexes a single simulation into the simdex.  All simdex
        attributes are updated, and the h5 file is completed with the variables
        defined in process.variables AND with the results of the postprocessing
        as defined in process.pp
        
        simulation has to be a Simulation object
        process is a simman.Process object
        
        Convention: in the h5 file, the short names are used, but with any '.'
        replaced by '_dot_'.  The h5 file only contains variables, no parameters.
        
        '''
        
        # internal function to enhance readibility
        def index_one_var(variables, varmap, var, index):
            """
            Updates the variables and varmap with var, keeping everything sorted
            Important: index is the position of the last found var in variables
            
            Returns the new variables, varmap and index
            """
            
            
            if var==variables[index]:
                # var is the first element in variables, update varmap
                varmap[index,-1] = 1
                pos = index+1
            else:
                try:
                    # search for it in variables, but only in the part AFTER index
                    pos=variables[index:].index(var)+index
                    varmap[pos,-1] = 1
                except(ValueError):
                    # this variable was not found.  Add it in the right position
                    # keeping the list in sorted order
                    pos = bisect.bisect_left(variables, var, lo=index)
                    variables.insert(pos,var)
                    # make new row in variablemap and add '1' in the last column
                    varmap = np.insert(varmap, pos, 0, axis=0)
                    varmap[pos,-1] = 1
                    pos+=1
            return variables, varmap, pos
                    
        # internal function to enhance readibility
        def index_one_par(parameters, parmap, parvalues, par, index, parvalue):
            """
            Updates the parameters, parvalues and parmap with par, 
            keeping everything sorted.
            Important: 
                - index is the position of the last found par in parameters
                - parvalue is the value of par 
            
            Returns the new parameters, parmap, parvalues and index
            """
          
            if par==parameters[index]:
                # par is the first element in parameters, update parmap, parvalues
                parmap[index,-1] = 1
                parvalues[index, -1] = parvalue
                pos = index+1
            else:
                try:
                    # search for it in parameters, but only in the part AFTER index
                    pos=parameters[index:].index(par)+index
                    parmap[pos,-1] = 1
                    parvalues[pos, -1] = parvalue
                except(ValueError):
                    # this parameter was not found.  Add it in the right position
                    # keeping the list in sorted order
                    pos = bisect.bisect_left(parameters, par, lo=index)
                    parameters.insert(pos,par)
                    # make new row in parametermap and add '1' in the last column
                    parmap = np.insert(parmap, pos, 0, axis=0)
                    parmap[pos,-1] = 1
                    parvalues = np.insert(parvalues, pos, 0, axis=0)
                    parvalues[pos,-1] = parvalue
                    pos+=1
            return parameters, parmap, parvalues, pos
        
        
        def add_meta(simulation, key):
            """Create a node for the simulation and add data to /Metadata"""
            
         
            
            class Meta(tbl.IsDescription):
                SID = tbl.StringCol(itemsize=16)
                path = tbl.StringCol(itemsize=160)
                log_analysed = tbl.BoolCol()
                successful = tbl.BoolCol()
                algorithm = tbl.StringCol(itemsize=16)
                cpu_time = tbl.Float32Col()
                successful_steps = tbl.Int32Col()
                steps_nok = tbl.Int32Col()
                timed_out = tbl.BoolCol()
                perc_wrong = tbl.Float32Col()
                time_events_model = tbl.Int32Col()
                time_events_U = tbl.Int32Col()
                state_events = tbl.Int32Col()
                step_events = tbl.Int32Col()
                step_size_min = tbl.Float32Col()
                step_size_max = tbl.Float32Col()
                int_order_max = tbl.Int32Col()
                
                
                
            self.openh5()
            
            # if it's the first simulation, we need to create the Metadata tbl
            try:
                meta = self.h5.getNode(self.h5.root.Metadata)
            except(tbl.NoSuchNodeError):
                meta = self.h5.createTable('/', 'Metadata', Meta, 
                            title='All metadata for the simulations')
            
            # check if there's a log file 
            logfilename = simulation.filename.replace('result_','dslog_')\
                                             .replace('.mat','.txt')
            try:
                log = analyse_log(logfilename)
                loganalysis=True
            except(IOError):
                print 'No %s found, log-analysis not possible' % logfilename
                loganalysis=False

            # create all values for a new row            
            row = meta.row
            row['SID'] = key
            row['path'] = simulation.filename
            row['log_analysed'] = loganalysis
               
            if loganalysis:
                for k in log:
                    row[k] = log[k]
                
            row.append()
            meta.flush
            
            # Create the node for all variable arrays
            var_grp = self.h5.createGroup('/', key, title='All variables, as arrays')
            
            self.h5.flush()
        
        
        def update_h5(simulation, key):
            """
            Update the h5 file with the variables defined by the process.
            Return a dictionary with shortname/longname pairs of everything
            that has been added to the h5.
            

            Still to add: extraction of metadata from the log
            """
            
            #pdb.set_trace()
            var_grp = self.h5.getNode('/', key)
           
            if process is None:
                # add all variables to the h5, with full names
                vardic = dict(zip(simulation.variables, simulation.variables))
                extracted = simulation.extract(var=vardic, arrays = 'each')           
                for shortname, arr in extracted.iteritems():
                    name = shortname.replace('.', '_dot_')
                    self.h5.createArray(var_grp, name, arr)
                
            else:
                extracted = simulation.postprocess(process)
                vardic = {}
                for shortname, arr in extracted.iteritems():
                    name = shortname.replace('.', '_dot_')
                    ispar = process.parameters.has_key(shortname) or \
                            process.parameters.has_key(name)
                                
                    if not ispar:
                        self.h5.createArray(var_grp, name, arr)
                        try:
                            longname = process.variables[shortname]
                        except(KeyError):
                            longname = shortname
                        vardic[shortname] = longname
                
            self.h5.flush()
            return vardic
        
        # separate parameters from variables for simulation 
        simulation.separate()
        
        if self.simulations == []:
            # this is the first simulation to be added to self
           key = self._gen_key()            
           self.simulations.append(key)
           self.files[key] = simulation.filename
           add_meta(simulation, key)
           vardic = update_h5(simulation, key)
               
           if self.verbose:
               print "key = %s, filename = %s" % (key, self.files[key])
           self.parameters = simulation.parameters # a LIST
           self.variables = simulation.variables  # a LIST           
           self.parametermap = np.ndarray((len(self.parameters), 1))
           self.parametermap[:, 0] = 1
           self.parametervalues = copy.copy(self.parametermap)
           self.parametervalues[:, 0] = np.array(simulation.parametervalues)
           self.variablemap = np.ndarray((len(self.variables), 1))
           self.variablemap[:, 0] = 1
           self.vardic = vardic
           
           self.h5.close()
        
        else:
            # new simulation to be added to existing ones            
            # First, add the simulation key to self.simulations            
            key = self._gen_key()            
            self.simulations.append(key)
            self.files[key] = simulation.filename
            add_meta(simulation, key)
            vardic = update_h5(simulation, key)
                      
            if self.verbose:
                print "Added simulation to set with at least one other simulation"
                print "key = %s, filename = %s" % (key, self.files[key])
            
            # Second, make new columns for parametermap, parametervalues and 
            # variablemap
            self.parametermap = np.append(self.parametermap,
                                          np.zeros((len(self.parameters), 1)),
                                          axis=1)
            self.parametervalues = np.append(self.parametervalues,
                                          np.zeros((len(self.parameters), 1)),
                                          axis=1)
            self.variablemap = np.append(self.variablemap,
                                          np.zeros((len(self.variables), 1)),
                                          axis=1)                                          
            
            position = 0            
            for var in simulation.variables:
                self.variables, self.variablemap, position = index_one_var(self.variables, self.variablemap, var, position)
            
            position = 0            
            for par, parvalue in zip(simulation.parameters, simulation.parametervalues):
                self.parameters, self.parametermap, self.parametervalues, position = index_one_par(self.parameters, self.parametermap, self.parametervalues, par, position, parvalue)
                
                  
            # finally, create or update self.vardic and self.pardic
            if process is not None:
                if not self.__dict__.has_key('vardic'):
                    # it's the first time we do a postprocessing on this simdex,
                    # probably because this was the first simulation to be indexed
                    self.vardic = vardic
                else:
                    self.vardic.update(vardic)
                               
                if process.parameters is not None:
                    try:
                        self.pardic.update(process.parameters)
                    except(AttributeError):
                        self.pardic=process.parameters

            # this method can be called on itself: close the h5 file afterwards
            self.h5.close()
            
    def filter_similar(self, SID):
        '''
        Return a new simdex with similar simulations as SID (SIDxxxx)        
        
        Create a new Simdex object from self with only those simulations 
        that have identical parameter and variable lists as SID.  This means
        that we suppose the model is identical, but some parameters may have 
        a different value
        '''
        
        # Approach: copy self and remove the unneeded columns from 
        # parametermap, parametervalues and variablemap by slicing
        
        # Make sure the h5 file is closed (for the deepcopy to work)
        self.h5.close()
        
        try:
            seqnb = self.simulations.index(SID)
        except(ValueError):
            print "This SID is not present in the simdex: %s" % SID
            raise
        
        newsimdex = copy.deepcopy(self)
        newsimdex.simulations = []
        
        parmap = self.parametermap[:, seqnb]
        varmap = self.variablemap[:, seqnb]
        
        sims_to_keep = []
        
        for i in range(len(self.simulations)):
            if np.all(self.parametermap[:, i] == parmap) and \
                np.all(self.variablemap[:, i] == varmap):
                # we have catched an identical simulation
                sims_to_keep.append(True)
                newsimdex.simulations.append(self.simulations[i])
            else:
                sims_to_keep.append(False)
            
        # slicing only works with an array
        s = np.array(sims_to_keep)
        newsimdex.parametermap = newsimdex.parametermap[ : , s]
        newsimdex.parametervalues = newsimdex.parametervalues[ : , s]
        newsimdex.variablemap = newsimdex.variablemap[ : , s]
        
        # remove all empty rows and corresponding parameters/variables
        newsimdex.cleanup()
        
        return newsimdex

    def filter_selection(self, selection):
        """
        Return a new simdex containing only the SID's in selection.
        
        selection is a list with 'SIDxxxx' strings.
        It DOES NOT have to be in the right order.
        """
        
        # Make sure the h5 file is closed (for the deepcopy to work)
        self.h5.close()
        
        newsimdex = copy.deepcopy(self)
        cols_to_remove = []
        
        for col, sid in enumerate(self.simulations):
            try:
                selection.index(sid)
            except(ValueError):
                # sid not in selection: remove it
                cols_to_remove.append(col)
                newsimdex.simulations.remove(sid)
                
        newsimdex.parametermap = np.delete(newsimdex.parametermap, 
                                           cols_to_remove, 1)
        newsimdex.parametervalues = np.delete(newsimdex.parametervalues, 
                                              cols_to_remove, 1)
        newsimdex.variablemap = np.delete(newsimdex.variablemap, 
                                          cols_to_remove, 1)
        
        newsimdex.cleanup()
        return newsimdex

    def filter_remove(self, selection):
        """
        Return a new simdex without the SID's in selection.
        
        selection is a list with 'SIDxxxx' strings.
        It DOES NOT have to be in the right order.
        """
        
        # Make sure the h5 file is closed (for the deepcopy to work)
        self.h5.close()
        
        newsimdex = copy.deepcopy(self)
        cols_to_remove = []
        
        for col, sid in enumerate(self.simulations):
            try:
                selection.index(sid)
            except(ValueError):
                # this sim is to be kept           
                pass
            else:
                # this simulation is found in selection ==> to be removed
                cols_to_remove.append(col)
                newsimdex.simulations.remove(sid)
                
        newsimdex.parametermap = np.delete(newsimdex.parametermap, 
                                           cols_to_remove, 1)
        newsimdex.parametervalues = np.delete(newsimdex.parametervalues, 
                                              cols_to_remove, 1)
        newsimdex.variablemap = np.delete(newsimdex.variablemap, 
                                          cols_to_remove, 1)
        
        newsimdex.cleanup()
        return newsimdex    
    
    
    def filter(self, pardic):
        '''
        Return a new simdex, filtered with the criteria as in pardic
        
        pardic is a dictionary of parameter:value pairs
        If a value is omitted (empty string), all simulations that have any 
        value for this parameter are fine
        
        Get all simulations that satisfy pardic (AND relation) and 
        return them as a new Simdex object
        
        Attention: a tolerance is defined internally in this method to enable
        filtering with floats.  The tolerance is currently set to 0.5% 
        '''
        
        # Approach: first find the parameters from pardic in self.parameters
        # and remove all rows from parametermap and parametervalues that aren't 
        # playing the game.
        # Then separate the parameters with values from the ones without
        
        # I select the rows by creating another array with the row numbers
        # and slice self.parametermap and self.parametervalues with that array
        
        tolerance = 0.005        
        # list of rownumbers from self.parameters with concerned parameters
        rows = []
        values = []    
        for i in pardic:
            rows.append(self.parameters.index(i))
            values.append(pardic[i])
        
        values = np.array(values)
        arows = np.array(rows)
        reduced_par_map = self.parametermap[arows]
        reduced_par_val = self.parametervalues[arows]
        
        
            
        # next step: remove all simulations that do NOT have the parameters
        # with empty strings in pardic
        
        # get row numbers of rows in reduced_par_map where the value doesn't 
        # matter
        parmaprows = np.array([x for (x, y) in zip(range(len(pardic)), \
            pardic.values()) if y==''])
        
        # if there are no parmaprows, we don't need to filter on empty strings
        if len(parmaprows)>0:
            selmap = reduced_par_map[parmaprows]
            # selmap contains rows with 0 and 1's for each of the simulations
            # (columns).  We need to get the simulation numbers (column numbers)
            # that are FINE, meaning that the columns are unit columns
            satisfyingmap = selmap.all(axis = 0)
            # satisfying is a boolean array, true if corresponding simulation 
            # is still in the run for selection
        else:
            # satisfyingmap has to be known.  
            # In this case a boolean array with only True values
            satisfyingmap = np.array(range(len(self.simulations)))>-1
        
        # now we need to get only the rows for which the values matter
        # and compare those values for each simulation with 'values'
        parvalrows = np.array([x for (x, y) in zip(range(len(pardic)), \
            pardic.values()) if y != ''])
        if len(parvalrows) > 0:
            selval = reduced_par_val[parvalrows]
            values = values[parvalrows]
            # we do not compare the values directly, but allow deviations 
            # smaller than tolerance.  
            abs_diff_values = np.abs(selval - values)
            satisfyingval = np.all(abs_diff_values < (tolerance * values), 
                                   axis = 0)
            # again a boolean array
        else:
            # satisfyingval has to be known.  
            # In this case a boolean array with only True values
            satisfyingval = np.array(range(len(self.simulations))) > -1
            
        
        # only simulations satisfying both the requirements are selected
        # first (dummy) row of self.simulations has also to be kept
        satisfying = satisfyingmap & satisfyingval
                
        # we create a new simdex object, with identical properties as self
        # but containing only the simulations we have selected
        
        newsimdex = copy.deepcopy(self)
        newsimdex.simulations = \
            [x for (x, y) in zip(self.simulations, satisfying) if y == True]
        newsimdex.parametermap = self.parametermap[:, satisfying]
        newsimdex.parametervalues = self.parametervalues[:, satisfying]
        newsimdex.variablemap = self.variablemap[:, satisfying]
            
        # we want to keep track of the parameters we have filtered on
        # this should be improved: if two identical keys occur, take the key
        # with associated value (instead of '')
        newsimdex.filterset.update(pardic)
        
        # Removing unused parameters and variables from the filtered simdex
        newsimdex.cleanup()
        if newsimdex.get_filenames() == []:
            raise ValueError("No single simulation could satisfy this filter")
        
        print newsimdex


        return newsimdex
    
    def cleanup(self):
        '''
        Removes unused parameters, variables and filenamesfrom a simdex
               
        '''
        # First, remove all columns from parametermap, parametervalues and 
        # variablemap that are not corresponding to self.simulations anymore
        
        
        
        new_files = {}
        for col, sid in enumerate(self.simulations):
            new_files[sid] = self.files[sid]
        self.files = new_files
        

        # next, remove all parameters/variables that are not in the maps anymore
        pars_to_keep = np.any(self.parametermap, 1)
        self.parametermap = self.parametermap[pars_to_keep]
        self.parametervalues = self.parametervalues[pars_to_keep]
        self.parameters = [x for (x, y) in \
            zip(self.parameters, pars_to_keep) if y == True]
        vars_to_keep = np.any(self.variablemap, 1)
        self.variables = [x for (x, y) in \
            zip(self.variables, vars_to_keep) if y == True]
        self.variablemap = self.variablemap[vars_to_keep]
        

    def get(self, name, aggregate=None):
        """
        Return a Result instance with SID:value pairs for par or var name
        
        If the name is a sub-variable, the corresponding variable for all
        mothers will be extracted and put in a single array according to aggregate.  
        
        aggregate = None, 'sum' or 'mean' : if None, the values in the Result
        object will contain the trajectories for all variables.  
        If aggregate is 'sum', all trajectories are summed, if it is 'mean', 
        the mean value of all trajectories is computed.
        
        If name is a parameter and a simulation does NOT have the parameter, 
        the value in the result object is None.
        
        """
        
        # There are many different options for name
        found_name = False        
        # 1. name is a short parameter name 
        # ==> attribute the longname to name
        try:
            if self.pardic.has_key(name):
                resdic = self._get_par(self.pardic[name])
                time = None
                found_name = True
        except(AttributeError):
            pass
        
        # 2. it is a short variable name
        # keep the short name cause that is the array name in the h5 file
        if not found_name:
            try:
                if self.vardic.has_key(name):
                    resdic = self._get_var_h5(name, selection=self.simulations)
                    time = self._get_var_h5('Time', selection=self.simulations)
                    found_name = True
            except(AttributeError):
                pass

        # 3. it is a long parameter name
        # keep it            
        if not found_name:
            try:
                parindex = self.parameters.index(name)
                resdic = self._get_par(name)
                time = None
                found_name = True
            except:
                pass
            
        # 4. Last regular option, it is a long variable name
        # ==> two options: if it is in a vardic, use the short name, 
        # else use it as it comes
        if not found_name:
            try:
                varindex = self.variables.index(name)
                found_name = True
            except:
                pass
            else:
                #it is a long variable name
                try:
                    for shortname, longname in self.vardic.iteritems():
                        if name == longname:
                            print 'shortname found', shortname
                            resdic = self._get_var_h5(shortname, 
                                                    selection=self.simulations)
                            time = self._get_var_h5('Time', 
                                                    selection=self.simulations)
                except(AttributeError):
                    pass
                else:
                    # it is a long variable name that is not in the h5
                    raise NotImplementedError('This variable name was not yet in the h5 file. \
                    \nAdapt the process to get it in there')
                    
                
        # 5. aggregation option: the name is a sub_var
        if not found_name and self.process.sub_vars.has_key(name):
            # we loop over the mothers and put all the arrays together
            for m in self.process.mothers:
                single_array = self._get_var_h5(var=m+'_'+name, 
                                                selection=self.simulations)
                if m == self.process.mothers[0]:
                    #initiate the resulting dictionary
                    resdic = copy.deepcopy(single_array)
                else:
                    for k in resdic.keys():
                        resdic[k] = np.column_stack((resdic[k], single_array[k]))
            
            if aggregate == 'sum':
                resdic = {k:np.sum(v, axis=1) for k,v in resdic.items()}
            elif aggregate == 'mean':
                resdic = {k:np.mean(v, axis=1) for k,v in resdic.items()}
            
            # reshape the array if the second dimension is larger than the first
            shape = resdic.values()[0].shape      
            try:
                if shape[0]==1 and shape[1] > 1:
                    resdic = {k:v.reshape((shape[1])) for k,v in resdic.items()}
            except:
                # the shape probably is of length 1 or even 0
                pass
            
            time = self._get_var_h5('Time', selection=self.simulations)
            found_name = True
  
        if not found_name:
            print "%s was not found in this simdex" % name
            print 'maybe you want to use any of these parameters/variables?'
            return self.exist(name)
        else:
            return Result(resdic, time=time, 
                          identifiers = self.identifiers, year=self.year)
        

        
    def _get_var_h5(self, var, selection=[]):
        """Get values of variables that are stored in the h5 file"""
        
        self.openh5()
        
        values = {}
        var_replaced = var.replace('.', '_dot_')
        if selection == []:
            selection = [n._v_name for n in self.h5.listNodes('/')]
            
        for node in self.h5.iterNodes('/'):
            try:
                # move on if this node is NOT in the selection
                selection.index(node._v_name)
                try:
                    # look up the variable in this node
                    array = self.h5.getNode(node, name=var_replaced)
                    values[node._v_name] = array.read()
                except(tbl.NoSuchNodeError):
                    # either the node is Metadata, or this variable does not
                    # exist in this node (perfectly possible and normal)
                    pass                    
                    #raise tbl.NoSuchNodeError(var + " not found in node " + node._v_pathname)
            except(ValueError):
                # it's a node that was not in the selection
                if self.verbose:
                    print " node not selected: ", node._v_name
                pass
            
            
       
        self.h5.close()
        return values        
    
    
    
    def _get_par(self, parameter):
        '''
        Return a dictionary with SID:parametervalue pairs
        
        This is a private method.  You should use get()
        parameter = string with exact (long) parameter name. 
        
        If a simulation does NOT have the parameter, the value in 
        the dictionary is None
        '''
        
        parindex = self.parameters.index(parameter)
            # row number of the parameter to be returned
        
#        result = [range(1, len(self.simulations)),\
#            self.parametervalues[parindex, 1:], self.simulations[1:]]
        
        presence = self.parametermap[parindex,:]        
        value = self.parametervalues[parindex, :]
        # take care, first element is dummy value (zero)
        
        result = {}
        for i, sid in enumerate(self.simulations):
            if presence[i] == 0:
                result[sid] = None
            elif presence[i] == 1:
                result[sid] = value[i]
            else:
                raise NotImplementedError('It is not possible to have something\
                  else than 0 or 1 in parametermap')
        
        return result

    def plot(self, variable):
        '''
        plot(variable) - variable = string with variable name (short or long)
        
        Creates a matplotlib figure with a simple plot of the timeseries for 
        each of the simulations in self
        '''
        
        result = self.get(variable)
        [fig, lines, leg] = result.plot(variable)

        return [fig, lines, leg]
        

    def scatterplot(self, X, Y):
        '''
        Creates a matplotlib figure with a simple plot of Y versus X for 
        each of the simulations in self 
        '''
        
           
        # 1. and 2.
        toplot_X = self.get(X).val
        toplot_Y = self.get(Y).val
        
        
        # 3. and 4.
        plotstring = ''
        plotlegend = ''
        
        for sid in toplot_X:
            if toplot_Y.has_key(sid):
                plotstring += ''.join(['toplot_X["', sid, '"], toplot_Y["', sid, '"], "D", '])
                try:
                    label = self.identifiers[sid]
                except KeyError:
                    label=sid
                plotlegend += ''.join(['"', label + '", '])        
        
        # remove last semicolon
        plotstring = plotstring[:-1]
        plotlegend = plotlegend[:-1]
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        lines = eval("ax.plot(" + plotstring + ")")
        leg = eval("ax.legend((" + plotlegend + "))")
        
        ax.set_xlabel(X)
        ax.set_ylabel(Y)
        
        return [fig, lines, leg]

    
    def get_SID(self, regex):
        '''
        Get a list with SID's for a given search expression (regex)
        
        regex = regular expression, not case sensitive
        
        Return a list of simulations (SID) of which the filenames match the regex
        '''
        p = re.compile(regex, re.IGNORECASE)
        matches = []
        sids = []
        # important: iterate always with self.simulations to keep the right order!!        
        for k in self.simulations:
            m = p.search(self.files[k])
            if m:
                matches.append(self.files[k])
                sids.append(k)
        
        print 'SID     ', 'Filename\n'
        for i, sim in zip(sids, matches):
            print i, '   ', sim
        return sids
        
    
    def save(self, filename):
        """
        save(filename)
        
        Save the Simdex object by pickling it with cPickle
        
        
        To unpickle (= load) use the following command:
            objectname = pickle.load(open(filename,'rb'))
            # 'rb' stands for 'read, binary'
            
        """
        
        # the h5 file raises errors, so let's remove it.  Anyway, we keep the 
        # reference to it via the h5_path string
        
        del self.h5
        print 'self.h5 removed'

        f = file(filename,'wb')
        # wb stands for 'write, binary'
        pickle.dump(self, f)
        f.close()
        
        return filename + ' created'
        
    def postproc(self):
        """Run the post-processing"""
        pass
    
    def apply(self, function_call):
        """
        Apply the function_call to each variable in each simulation.
        E.g, if the simdex contains a variable called QHeat as a timeseries,
        you can use this method like this to compute the total integrated QHeat
        for each of the simulations:
            
            simdex.apply(np.trapz(QHeat, time))
            
        Returns a dictionary with SID/result pairs
        
        """
        result={}
        for SID in self.simulations:
            result[SID] = eval(function_call)
                
        return result
        
    def apply2(self, function, variable):
        result={}
        for k,v in self.get(variable).items():
            result[k] = function(v)
        return result
        
    def apply3(self, function, variable, *args, **kwargs):
        result={}
        for k,v in self.get(variable).items():
            result[k] = function(v, *args, **kwargs)
        return result
        
    def apply4(self, function, **kwargs):
        
        arg_dict = {fun_arg: self.get(sim_args) for fun_arg, sim_args in kwargs.items()}
        result = {}
        for SID in arg_dict[kwargs.keys()[0]]:
            fun_kwargs = {fun_arg: sim_dict[SID] for fun_arg, sim_dict in arg_dict.items()} 
            result[SID] = function(**fun_kwargs)
        return result
 

        
def apply(function, results):
    """
    Apply the function on each of the values in results.
    Results is a dictionary with SID/value pairs, typically as a result 
    from a simdex.get() call. 
    
    Returns a SID/function(value) dictionary    
    """
        
    result = {}
    for k,v in results.items():
        result[k] = function(v)
    return result


def load_simdex(filename):
    """load and return a previously saved Simdex object"""
    
    result = pickle.load(open(filename,'rb'))
    result.h5 = tbl.openFile(result.h5_path, 'a')
    result.h5.close()
    return result