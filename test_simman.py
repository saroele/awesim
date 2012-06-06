# -*- coding: utf-8 -*-
"""
Test module for the simman module

Created on Thu Jan 06 20:46:26 2011

@author: RDC
"""

import numpy as np
import unittest
from os import getcwd, path, remove
from cStringIO import StringIO
import sys
import matplotlib
from simman import Simulation, Simdex, Result, Process, load_simdex


class ProcessTest(unittest.TestCase):
    """
    Class for testing the class simman.Process
    """

    def setUp(self):
        self.mothers=['c1', 'c2']
        self.parameters={'cap1':'c1.C', 'res':'r.R'}
        self.sub_pars={'cap':'C'}
        self.variables={}
        self.sub_vars={'Qflow':'heatPort.Q_flow'}
        self.pp=['Qflow10 = 10 * Qflow']
        

    def test_init1(self):
        """init process without attributes"""
        
        p = Process()
        self.assertEqual(p.mothers, [])
        
    def test_init2(self):
        """init Process with only variables should work"""
        process = Process(mothers=['c1', 'c2'], sub_vars={'T':'T'})
        self.assertEqual(process.variables, {'Time':'Time',
                                             'c1_T':'c1.T', 'c2_T':'c2.T'})
                        

    def test_init3(self):
        """init process with variables and parameters"""
        
        p = Process(mothers=self.mothers, parameters=self.parameters, 
                    sub_pars=self.sub_pars, variables=self.variables, 
                    sub_vars=self.sub_vars, pp=self.pp)
        self.assertEqual(p.parameters, {'cap1':'c1.C', 'res':'r.R',
                                        'c1_cap':'c1.C', 'c2_cap':'c2.C'})        
        self.assertEqual(p.variables, {'Time':'Time', 
                                       'c1_Qflow':'c1.heatPort.Q_flow',
                                       'c2_Qflow':'c2.heatPort.Q_flow'}) 
        print p


    def test_init_4(self):
        """Process initialisation should not change original variables"""
        
        variables = {}  
        parameters={'cap1':'c1.C', 'res':'r.R'}
        sub_pars={'cap':'C'} 
        process = Process(mothers=['c1', 'c2'], variables = variables,
                          parameters=parameters, sub_pars=sub_pars,
                          sub_vars={'T':'T'})
        self.assertEqual(process.variables, {'Time':'Time',
                                             'c1_T':'c1.T', 'c2_T':'c2.T'})
        self.assertEqual(variables, {})
        self.assertEqual(parameters, {'cap1':'c1.C', 'res':'r.R'})


    def test_init_integration(self):
        """init process with variables to integrate"""
        
        J2kWh = 1e-6/3.6
        vars_to_integrate = {'Qflow':J2kWh}
                             
        p = Process(mothers=self.mothers, parameters=self.parameters, 
                    sub_pars=self.sub_pars, variables=self.variables, 
                    sub_vars=self.sub_vars, pp=self.pp, integrate=vars_to_integrate)
        self.assertEqual(p.parameters, {'cap1':'c1.C', 'res':'r.R',
                                        'c1_cap':'c1.C', 'c2_cap':'c2.C'})        
        self.assertEqual(p.variables, {'Time':'Time', 
                                       'c1_Qflow':'c1.heatPort.Q_flow',
                                       'c2_Qflow':'c2.heatPort.Q_flow'}) 
        self.assertEqual(p.pp, ['Qflow_Int = np.trapz( Qflow , Time ,axis=0)*'+str(J2kWh)+' if Qflow .shape[0]== Time .shape[0] else np.array([0.0])',
                                'Qflow10 = 10 * Qflow'])
                                       


class ResultTest(unittest.TestCase):
    """Class for testing Result"""
    
    def setUp(self):
        self.values = {'c':np.arange(4), 
                       'a':np.array([1,3,5,7]), 
                       'b':np.array([2,8])}
        self.time = {'c':np.arange(4), 
                       'a':np.arange(4), 
                       'b':np.array([0,1])}
        self.identifiers = {'c':'sim C', 
                       'a':'sim A', 
                       'b':'sim B'}
    
    
    def test___init__(self):
        """Initiation of a Result with only values"""
        res = Result(self.values)
        
        self.assertEqual(res.val, self.values)
        self.assertFalse(hasattr(res, 'time'))
        self.assertFalse(hasattr(res, 'identifiers'))
        
    def test___init__2(self):
        """Initiation of a Result with time and identifiers"""
        res = Result(self.values, self.time, self.identifiers)
        
        self.assertEqual(res.val, self.values)
        self.assertTrue(hasattr(res, 'time'))
        self.assertTrue(hasattr(res, 'identifiers'))
        
    def test___init__3(self):
        """Initiation of a Result with additional **kwargs"""
        res = Result(self.values, self.time, self.identifiers, year=2010)
        
        self.assertEqual(res.val, self.values)
        self.assertTrue(hasattr(res, 'time'))
        self.assertTrue(hasattr(res, 'identifiers'))
        
        self.assertEqual(res.year, 2010)

    def test_values(self):
        """Get values correctly sorted"""
        res = Result(self.values, self.time, self.identifiers)
        v=res.values()
        self.assertTrue(np.all(v[0]==np.array([1,3,5,7])))
        self.assertTrue(np.all(v[1]==np.array([2,8])))
        self.assertTrue(np.all(v[2]==np.arange(4)))
        
    def test_trapz(self):
        """Integrate the values according to provide time"""
        res = Result(self.values, self.time, self.identifiers)
        v=res.trapz()
        
        self.assertEqual(v[0], np.trapz(np.array([1,3,5,7]), x=np.arange(4)))
        self.assertEqual(v[1], np.trapz(np.array([2,8]), x=np.array([0,1])))
        self.assertEqual(v[2], np.trapz(np.arange(4), x=np.arange(4)))
        

class SimulationTest(unittest.TestCase):
    """
    Class for testing the class simman.Simulation
    """
    
    def test___init__attributes(self):
        """
        Tests creation of Simulation object from well known file.
        Checks if attributes are present 
        Make sure the file LinkedCapacities.mat is in the current work directory

        """
        sim = Simulation('LinkedCapacities')
        for attr in ['dataInfo', 'data_1', 'data_2', 'filename', 'names']:
            self.assertTrue(sim.__dict__.has_key(attr), 
                            'Simulation.__init__() dit not create \
                            attribute %s' % attr)
                            
    def test___init__subfolder(self):
        """
        Tests creation of Simulation object from well known file with full
        pathname.
        Checks if attributes are present 
        Make sure the file LinkedCapacities.mat is in the folder 'Subfolder' 
        of the current work directory

        """
        cwd = getcwd()
        filename = path.join(cwd, 'Subfolder', 'LinkedCapacities')
        sim = Simulation(filename)
        for attr in ['dataInfo', 'data_1', 'data_2', 'filename', 'names']:
            self.assertTrue(sim.__dict__.has_key(attr), 
                            'Simulation.__init__() dit not create \
                            attribute %s' % attr)

    def test___init__subfolder_with_spaces(self):
        """
        Tests creation of Simulation object from well known file with full
        pathname including spaces.
        Checks if attributes are present 
        Make sure the file LinkedCapacities.mat is in the folder 
        'Subfolder with spaces' of the current work directory

        """
        cwd = getcwd()
        filename = path.join(cwd, 'A Subfolder with Spaces', 'LinkedCapacities')
        sim = Simulation(filename)
        for attr in ['dataInfo', 'data_1', 'data_2', 'filename', 'names']:
            self.assertTrue(sim.__dict__.has_key(attr), 
                            'Simulation.__init__() dit not create \
                            attribute %s' % attr)

                            
    def test___init__attributes_array(self):
        """
        Tests creation of Simulation object from well known file with arrays.
        Checks if attributes are present 
        Make sure the file Array.mat is in the current work directory

        """
        sim = Simulation('Array')
        for attr in ['dataInfo', 'data_1', 'data_2', 'filename', 'names']:
            self.assertTrue(sim.__dict__.has_key(attr), 
                            'Simulation.__init__() dit not create \
                            attribute %s' % attr)
                            
    def test___init__wrongmatfile(self):
        """
        Test if an IOError is raised when a wrong type of .mat file is used
        
        """
        # important not to forget the syntax of assertRaises. You have to 
        # give a callable object, so don't call it.  Arguments can be added
        # after the next , as shown below
        self.assertRaises(IOError, Simulation, 'EmptyMatFile.mat')
            
    
    def test_exist(self):
        """
        Test if :
        1) running Simulation.exist() gives a match when it should
           (different cases are tested)
        2) running Simulation.exist() gives a zero list for non-present regex
        
        """
        
        sim = Simulation('LinkedCapacities')
        self.assertEqual(sim.names[1], sim.exist(sim.names[1])[0], 
                         'Simulation.exist() does NOT return a present name')
        self.assertEqual(5, len(sim.exist('c1')), "sim.exist('c1') should give \
                         5 matches")
        self.assertEqual(u'Time', sim.exist('time')[0], 
                         'Simulation.exist() should ignore case' )
        self.assertEqual([], sim.exist('something that does not exist'))

    def test_exist_array(self):
        """
        For a simulation file which uses arrays, Test if :
        1) running Simulation.exist() gives a match when it should
           (different cases are tested)
        2) running Simulation.exist() gives a zero list for non-pressent regex
        
        """
        
        sim = Simulation('Array')
        self.assertEqual(sim.names[1], sim.exist(sim.names[1])[0], 
                         'Simulation.exist() does NOT return a present name')
        self.assertEqual(4, len(sim.exist('heatport.q_flow')), 
                         "sim.exist('heatport.q_flow') should give 4 matches")
        self.assertEqual(12, len(sim.exist('[1]')), 
                         "sim.exist('[1]') should give 12 matches")
        
                         
    def test_get_value_wrongname(self):
        """ Tests get_value with a wrong or empty name """
        
        sim = Simulation('LinkedCapacities')
        self.assertRaises(ValueError, sim.get_value, 'wrongname')
        self.assertRaises(ValueError, sim.get_value, '')

    def test_get_value_par(self):
        """ Tests get_value on a parameter  """
        
        sim = Simulation('LinkedCapacities')
        self.assertEqual(600.0, sim.get_value('c1.C'), 'c1.C should be 600.0')
        
    def test_get_value_time1(self):
        """ 
        Tests get_value on the variable 'Time'
        The testes file is the result of a simulation over 10000s with fixed 
        interval of 200s 
        """
        
        sim = Simulation('LinkedCapacities')
        time = np.array(np.arange(0., 10200., 200.))
        sim_time = sim.get_value('Time')
        
        self.assertTrue((time == sim_time).all(), 
                        'Time should contain exactly 51 values, \
                        from 0 to 10000 (incl.) in steps of 200')
                        
    def test_get_value_time2(self):
        """ 
        Tests get_value on the variable 'Time'
        The testes file is the result of a simulation over 10000s with fixed 
        amount of intervals (50) 
        """
        
        sim = Simulation('Array')
        time = np.array(np.arange(0., 10200., 200.))
        Time = sim.get_value('Time')
        
        self.assertTrue((time == Time).all(), 
                        'Time should contain exactly 51 values, \
                        from 0 to 10000 (incl.) in steps of 200')
        
    def test_separate_attributes_present(self):
        """ Tests if the right attributes are created """
        
        sim = Simulation('LinkedCapacities')    
        sim.separate()
        
        self.assertTrue(isinstance(sim.parameters, list), 
                        'Simulation.separate() should create attribute \
                        parameters as list')
        self.assertTrue(isinstance(sim.variables, list), 
                        'Simulation.separate() should create attribute \
                        variables as list')
        self.assertTrue(isinstance(sim.parametervalues, np.ndarray), 
                        'Simulation.separate() should create attribute \
                        parametervalues as numpy array')                        

    def test_separate_attributes_correct(self):
        """ Tests if the created attributes are correct """
        
        sim = Simulation('LinkedCapacities')    
        sim.separate()
        
        parameters = np.array([600.0, 1000.0, 3])
        

        self.assertEqual([u'c1.C', u'c2.C', u'r.R'], sim.parameters)
        self.assertTrue((parameters == sim.parametervalues).all())
        variables = sorted([u'Time', u'c1.heatPort.T', u'c1.heatPort.Q_flow',
                          u'c1.T', u'c1.der(T)', u'c2.heatPort.T',
                          u'c2.heatPort.Q_flow', u'c2.T', u'c2.der(T)',
                          u'r.heatPort_a.T',  u'r.heatPort_a.Q_flow',
                          u'r.heatPort_b.T', u'r.heatPort_b.Q_flow'])       
        self.assertEqual(variables, 
                          sim.variables) 
                          
    def test_separate_twice(self):
        """ Tests if all goes fine if Simulation.separate() is called twice """
        
        sim = Simulation('LinkedCapacities')    
        sim.separate()
        sim.separate()
        
        parameters = np.array([600.0, 1000.0, 3])
        

        self.assertEqual([u'c1.C', u'c2.C', u'r.R'], sim.parameters)
        self.assertTrue((parameters == sim.parametervalues).all())
        variables = sorted([u'Time', u'c1.heatPort.T', u'c1.heatPort.Q_flow',
                          u'c1.T', u'c1.der(T)', u'c2.heatPort.T',
                          u'c2.heatPort.Q_flow', u'c2.T', u'c2.der(T)',
                          u'r.heatPort_a.T',  u'r.heatPort_a.Q_flow',
                          u'r.heatPort_b.T', u'r.heatPort_b.Q_flow'])        
        
        
        self.assertEqual(variables, sim.variables)
        
    def test_extract(self):
        """to be implemented"""
        pass

    def test_get_objects(self):
        """Test if get_objects works with empty mother model"""
        
        sim = Simulation('LinkedCapacities') 
        obj = sim.get_objects()
        obj_sorted = sorted(obj)
        self.assertEqual(obj_sorted, sorted([u'c1', u'c2','r']))

    def test_get_objects_with_mother(self):
        """Test if get_objects works with empty mother model"""
        
        sim = Simulation('LinkedCapacities') 
        obj = sim.get_objects(mother = 'c1')
        obj_sorted = sorted(obj)
        self.assertEqual(obj_sorted, sorted(['heatPort', 'C', u'T', u'der(T)']))                

    def test_postprocess_nopp(self):
        """Try a basic extraction via a post processing on a simulation"""
        
        sim = Simulation('LinkedCapacities')
        process = Process(mothers=['c1', 'c2'], sub_vars={'T':'T'})
        result_pp = sim.postprocess(process)
        self.assertTrue(result_pp.has_key('c1_T'))
        
    def test_postprocess_mothers(self):
        """Try a basic postprocessing with mothers on a simulation"""
        
        sim = Simulation('LinkedCapacities')
        process = Process(mothers=['c1', 'c2'], sub_vars={'T':'T'},
                          sub_pars={'cap':'C'},
                          pp = ['T_degC = T + 273.15'])
        result_pp = sim.postprocess(process)
        self.assertTrue(np.all(result_pp['c1_T_degC'] == result_pp['c1_T']+273.15))
        self.assertEqual(result_pp['c1_cap'], 600.0)
        
    def test_postprocess_nomothers(self):
        """Try a basic postprocessing without mothers on a simulation"""
        
        sim = Simulation('LinkedCapacities')
        process = Process(mothers=['c1', 'c2'], sub_vars={'T':'T'},
                          pp = ['timehours = Time / 3600'])
        result_pp = sim.postprocess(process)
        self.assertTrue(np.all(result_pp['timehours']==result_pp['Time']/3600))

    def test_postprocess_recurring(self):
        """Some advanced postprocessing using previous results with mothers"""
        
        sim = Simulation('LinkedCapacities')
        process = Process(mothers=['c1', 'c2'], sub_vars={'T':'T'},
                          sub_pars={'cap':'C'},
                          pp = ['T_degC = T + 273.15', 'T_degC2 =  T_degC '])
        result_pp = sim.postprocess(process)
        self.assertTrue(np.all(result_pp['c1_T_degC2']==result_pp['c1_T_degC']))

    def test_postprocess_advanced(self):
        """Some more advanced postprocessing with mothers"""
        
        sim = Simulation('LinkedCapacities')
        process = Process(mothers=['c1', 'c2'], sub_vars={'T':'T', 'Q':'heatPort.Q_flow'},
                          sub_pars={'cap':'C'},
                          pp = ['T_degC = T + 273.15', 
                                'T_max =  np.amax( T_degC )',
                                'Thigh = np.nonzero( T_degC > 640)[0]',
                                'Qsel = Q [ Thigh ]'])
        result_pp = sim.postprocess(process)
        self.assertAlmostEqual(result_pp['c1_T_max'], 673.15, 2)
        self.assertEqual(len(result_pp['c1_Qsel']), len(result_pp['c1_Thigh']))

    def test_postprocess_integration(self):
        """Postprocessing with integration on standard variables and mothers"""
        
        J2kWh = 1e-6/3.6
        vars_to_integrate = {'Q':J2kWh, 'Time':1}
        
        sim = Simulation('LinkedCapacities')
        process = Process(mothers=['c1', 'c2'], sub_vars={'T':'T', 'Q':'heatPort.Q_flow'},
                          sub_pars={'cap':'C'},
                          pp = ['T_degC = T + 273.15', 
                                'T_max =  np.amax( T_degC )',
                                'Thigh = np.nonzero( T_degC > 640)[0]',
                                'Qsel = Q [ Thigh ]'],
                          integrate = vars_to_integrate)
        result_pp = sim.postprocess(process)
        self.assertAlmostEqual(result_pp['c1_T_max'], 673.15, 2)
        self.assertEqual(len(result_pp['c1_Qsel']), len(result_pp['c1_Thigh']))
        self.assertAlmostEqual(result_pp['c1_Q_Int'], -result_pp['c2_Q_Int'], 10)
        self.assertIsNotNone(result_pp['Time_Int'])     
        

        
class SimdexTest(unittest.TestCase):
    """
    Class for testing the class simman.Simdex
    """

    def setUp(self):
      
        # The next N times, that a question is being asked interactively,
        # we will automatically answer y.
        # 
        # It doesn't matter if N exceeds the numer of questions
        # NOTE: this may affect interactive debugging tools
        N = 1000
        f = StringIO("y\n" * N)
        sys.stdin = f
        
        self.cwd = getcwd()
        # sims contains the simualations we expect to be in self.simdex in the 
        # current work directory        
        self.sims = []
        self.filenames = ['Array.mat', 'LinkedCapacities.mat', \
                   'LinkedCapacities_A.mat',  'LinkedCapacities_B.mat', \
                   'LinkedCapacities_C.mat', 'LinkedCapacities_D.mat', \
                   'LinkedCapacities_E.mat', 'LinkedCapacities_F.mat']
        self.filenames.sort()
        for fn in self.filenames:
            self.sims.append(path.join(self.cwd, fn))
            
        self.simdex = Simdex()
        self.simdex.scan()    
                       
    def tearDown(self):
        """ Restore back interactivity of keyboard and close h5 file """
        sys.stdin = sys.__stdin__
        try:
            self.simdex.h5.close()
        except:
            pass
        if path.exists(path.join(getcwd(), 'simdex.h5')):
            remove(path.join(getcwd(), 'simdex.h5'))
        
    
    def test_init(self):
        """
        Tests if a Simdex object is created correctly based on the following
        files in the current work directory, containing:
            - 'Array.mat',
            - 'LinkedCapacities.mat',
            - 'LinkedCapacities_A.mat',
            - 'LinkedCapacities_B.mat',
            - 'LinkedCapacities_C.mat',
            - 'LinkedCapacities_D.mat',
            - 'LinkedCapacities_E.mat'
            - 'LinkedCapacities_F.mat']
      
        """
        filenames = self.simdex.get_filenames('path')
        filenames.sort()
        self.assertEqual(self.sims, filenames)
        self.simdex.h5.close()

    def test_init_with_process(self):
        """
        Tests if a Simdex object is created correctly when a process is passed
      
        """
        
        process=Process(variables={'T2':'c2.T', 'dt2':'c[2].der(T)'})
        self.simdex=Simdex(folder = getcwd(), process=process)        
        
        filenames = self.simdex.get_filenames('path')
        filenames.sort()
        self.assertEqual(self.sims, filenames)
        self.assertEqual(self.simdex.vardic, 
                         {'Time':'Time','T2':'c2.T', 'dt2':'c[2].der(T)'})

        self.simdex.h5.close()        
        
    def test_init_with_process_pp(self):
        """Tests if vardic is updated correctly after postprocessing"""
        
        J2kWh = 1e-6/3.6
        vars_to_integrate = {'Q':J2kWh, 'Time':1}        
        process = Process(mothers=['c1', 'c2'], 
                          sub_vars={'T':'T', 'Q':'heatPort.Q_flow'},
                          sub_pars={'cap':'C'},
                          pp = ['T_degC = T + 273.15', 
                                'T_max =  np.amax( T_degC )',
                                'Thigh = np.nonzero( T_degC > 640)[0]',
                                'Qsel = Q [ Thigh ]'],
                          integrate = vars_to_integrate)
        
        self.simdex=Simdex(folder = getcwd(), process=process)        
        
        filenames = self.simdex.get_filenames('path')
        filenames.sort()
        self.assertEqual(self.sims, filenames)
        for v in ['c1_Thigh', 'c2_Q', 'c2_Q_Int']:
            self.assertTrue(self.simdex.vardic.has_key(v))
        

        self.simdex.h5.close()


    def test_init_subfolder_with_crappy_files(self):
       """ Test initiation from a folder including wrong .mat files"""
       
       folder = path.join(self.cwd, 'SubfolderWithCrappyFiles')
       self.simdex = Simdex(folder)
       filenames = self.simdex.get_filenames()
       filenames.sort()
       self.assertEqual(self.filenames, filenames)
       self.simdex.h5.close()
       
    def test_scan_additional_folder(self):
        """Test adding files from another folder to existing simdex"""
        
        J2kWh = 1e-6/3.6
        vars_to_integrate = {'Q':J2kWh, 'Time':1}        
        process = Process(mothers=['c1', 'c2'], 
                          sub_vars={'T':'T', 'Q':'heatPort.Q_flow'},
                          sub_pars={'cap':'C'},
                          pp = ['T_degC = T + 273.15', 
                                'T_max =  np.amax( T_degC )',
                                'Thigh = np.nonzero( T_degC > 640)[0]',
                                'Qsel = Q [ Thigh ]'],
                          integrate = vars_to_integrate)
        
        self.simdex=Simdex(folder = getcwd(), process=process)
        
        # now add another folder (with exactly the same files)
        folder = path.join(self.cwd, 'SubfolderWithCrappyFiles')
        self.simdex.scan(folder = folder, process=process)
        self.assertEqual(len(self.simdex.simulations), 16)


    def test_exist(self):
        """
        Test if :
        1) running Simdex.exist() gives a match when it should
           (different cases are tested)
        2) running Simdex.exist() gives a zero list for non-present regex
        
        """
        
        self.assertEqual(u'c1.C', 
                         self.simdex.exist('C1.c')[0][0], 
                         'Simdex.exist() does NOT return a present name')
        self.assertEqual(1, len(self.simdex.exist('c1')[0]), 
                         "sim.exist('c1') should return 1 parameter")
        self.assertEqual(1, len(self.simdex.exist('c1', 'par')), 
                         "sim.exist('c1', 'par') should return 1 parameter")                 
        self.assertEqual(4, len(self.simdex.exist('c1')[1]), 
                         "sim.exist('c1') should return 4 variables")
        self.assertEqual(4, len(self.simdex.exist('c1', 'var')), 
                         "sim.exist('c1', 'var') should return 4 variables")                 
        self.assertEqual(u'Time', self.simdex.exist('time')[1][0], 
                         'Simulation.exist() should ignore case' )
        self.assertEqual([[],[]], self.simdex.exist('this does not exist'))
        self.simdex.h5.close()
        
    def test_cleanup_1(self):
        """Check if parameters are removed if not present"""
                
        
        n = 1
        # set 1 value in the parametermap to 0
        par = self.simdex.parameters[n]
        npars = len(self.simdex.parameters)
        nvars = len(self.simdex.variables)
        self.simdex.parametermap[n, 0] = 0
        self.simdex.cleanup()
        self.assertEqual(npars-1, len(self.simdex.parameters))
        self.assertEqual(npars-1, self.simdex.parametermap.shape[0])
        self.assertEqual(npars-1, self.simdex.parametervalues.shape[0])
        self.assertEqual(nvars, len(self.simdex.variables))
        
        self.simdex.variablemap[44, 0] = 0
        self.simdex.cleanup()
        self.assertEqual(nvars-1, len(self.simdex.variables))
        self.assertEqual(nvars-1, self.simdex.variablemap.shape[0])
        self.assertEqual(npars-1, self.simdex.parametervalues.shape[0])
        self.simdex.h5.close()


    def test_filter_selection(self):
        """Simdex.filter_selection() should return only selected SIDs"""
        
        selection = ['SID0001', 'SID0003', 'SID0004']
        c1_C = self.simdex.get('c1.C')
        self.simdex_filtered = self.simdex.filter_selection(selection)
        exp_results = ['LinkedCapacities.mat', \
                   'LinkedCapacities_B.mat', 'LinkedCapacities_C.mat']
        exp_results.sort()
        self.simdex_filtered_fn = self.simdex_filtered.get_filenames()
        self.simdex_filtered_fn.sort()
        self.assertEqual(exp_results, self.simdex_filtered_fn, 
                         'get_identical on LinkedCapacities_C.mat should \
                         return the files mentioned in exp_results')
        c1_C_f = self.simdex_filtered.get('c1.C')                 
        for sid in c1_C_f.val:
            self.assertEqual(c1_C.val[sid], c1_C_f.val[sid])
                 
    def test_filter_remove(self):
        """Simdex.filter_remove() should remove only selected SIDs"""
        
        selection = ['SID0001', 'SID0003', 'SID0004', 'SID0006']
        self.simdex_filtered = self.simdex.filter_remove(selection)
        exp_results = ['Array.mat', 
                   'LinkedCapacities_A.mat',  'LinkedCapacities_D.mat', \
                   'LinkedCapacities_F.mat']
        exp_results.sort()
        self.simdex_filtered_fn = self.simdex_filtered.get_filenames()
        self.simdex_filtered_fn.sort()
        self.assertEqual(exp_results, self.simdex_filtered_fn, 
                         'get_identical on LinkedCapacities_C.mat should \
                         return the files mentioned in exp_results')
        
    
    def test_filter_remove_2(self):
        """Check if colums are removed from parametermap if needed"""
                
        # we test the shape and the values
        shape = self.simdex.parametermap.shape
        c1_C = self.simdex.get('c1.C')
        filtered = self.simdex.filter_remove(['SID0002'])
        shape_filtered = filtered.parametermap.shape
        self.assertEqual(shape[1]-1, shape_filtered[1])
        c1_C_f = filtered.get('c1.C')
        for sid in c1_C_f.val:
            self.assertEqual(c1_C.val[sid], c1_C_f.val[sid])


        self.simdex.h5.close()  
        filtered.h5.close()


    def test_filter_similar_single_result(self):
        """ Simdex.get_identical() for Array.mat should return only Array.mat"""
        

        SID_array = self.simdex.get_SID('array')
        self.simdex_array = self.simdex.filter_similar(SID_array[0])
        self.assertEqual(['Array.mat'], self.simdex_array.get_filenames(), 
                         'get_identical on array should only return Array.mat')
        
        # check if the self.simdex and sim objects have the same parameters
        sim_array = Simulation('Array.mat')
        sim_array.separate()
        sim_array_pars = sim_array.parameters
        sim_array_pars.sort()
        self.simdex_array_pars = self.simdex_array.parameters
        self.simdex_array_pars.sort()
        self.assertEqual(self.simdex_array_pars, sim_array_pars,
                         'the self.simdex and sim objects should have the same \
                         parameters')
        
        # check if the self.simdex and sim objects have the same variables
        sim_array_vars = sim_array.variables
        sim_array_vars.sort()
        self.simdex_array_vars = self.simdex_array.variables
        self.simdex_array_vars.sort()
        self.assertEqual(self.simdex_array_vars, sim_array_vars,
                         'the self.simdex and sim objects should have the same \
                         variables')
        self.simdex.h5.close()

    def test_filter_similar_multiple_results(self):
        """ 
        Simdex.get_identical() for LinkedCapacities_C.mat should return 
        all LinkedCapacities* except LinkedCapacities_F 
        
        Attention: currently, the String parameter that is defined in the 
        Modelica files is NOT present in the .mat files.  If this would be 
        changed (desirable) than the result of this test should also 
        exclude LinkedCapacities.mat
                
        """
        

        SID_lc = self.simdex.get_SID('_C')
        
        self.simdex_lc = self.simdex.filter_similar(SID_lc[0])
        exp_results = ['LinkedCapacities.mat', \
                   'LinkedCapacities_A.mat',  'LinkedCapacities_B.mat', \
                   'LinkedCapacities_C.mat', 'LinkedCapacities_D.mat', \
                   'LinkedCapacities_E.mat']
        exp_results.sort()
        self.simdex_lc_fn = self.simdex_lc.get_filenames()
        self.simdex_lc_fn.sort()
        self.assertEqual(exp_results, self.simdex_lc_fn, 
                         'get_identical on LinkedCapacities_C.mat should \
                         return the files mentioned in exp_results')
        
        # check if the self.simdex and sim objects have the same parameters
        sim_lc = Simulation('LinkedCapacities_A.mat')
        sim_lc.separate()
        sim_lc_pars = sim_lc.parameters
        sim_lc_pars.sort()
        self.simdex_lc_pars = self.simdex_lc.parameters
        self.simdex_lc_pars.sort()
        self.assertEqual(self.simdex_lc_pars, sim_lc_pars,
                         'the self.simdex and sim objects should have the same \
                         parameters')
        
        # check if the self.simdex and sim objects have the same variables
        sim_lc_vars = sim_lc.variables
        sim_lc_vars.sort()
        self.simdex_lc_vars = self.simdex_lc.variables
        self.simdex_lc_vars.sort()
        self.assertEqual(self.simdex_lc_vars, sim_lc_vars,
                         'the self.simdex and sim objects should have the same \
                         variables')
                         
        self.simdex.h5.close()
    
  
    
    def test_filter_intvalues(self):
        """Simdex.filter() with integer values should work well"""
        

        filt_dic = {'c1.C': 800}
        self.simdex_filtered = self.simdex.filter(filt_dic)
        self.simdex_filtered_fn = self.simdex_filtered.get_filenames()
        self.simdex_filtered_fn.sort()
        exp_result = ['LinkedCapacities_A.mat',  'LinkedCapacities_D.mat', \
                   'LinkedCapacities_E.mat', 'LinkedCapacities_F.mat']
        exp_result.sort()
        self.assertEqual(exp_result, self.simdex_filtered_fn,
                         'filtering self.simdex with c1.C = 800 should return\
                         exp_result')
        self.assertEqual(filt_dic, self.simdex_filtered.filterset,
                         'After filtering, filterset has to be updated')
        self.simdex.h5.close()
        
       
    def test_filter_nonexist_values(self):
        """Simdex.filter() with non-existent values should return False"""
        pass
        
        filt_dic = {'c1.C': 850} 
        filt_dic2 = {'c1.C': 800, 'r.R': 800}
        self.assertRaises(ValueError, self.simdex.filter, filt_dic)
        self.assertRaises(ValueError, self.simdex.filter, filt_dic2)
        self.simdex.h5.close()
        
    def test_filter_wildcard(self):
        """Simdex.filter() with '' values should return any sim having \
        that parameter"""
        
        filt_dic = {'c1.C': ''}
        self.simdex_filtered = self.simdex.filter(filt_dic)
        self.simdex_filtered_fn = self.simdex_filtered.get_filenames()
        self.simdex_filtered_fn.sort()
        exp_result = ['LinkedCapacities.mat', \
                   'LinkedCapacities_A.mat',  'LinkedCapacities_B.mat', \
                   'LinkedCapacities_C.mat', 'LinkedCapacities_D.mat', \
                   'LinkedCapacities_E.mat', 'LinkedCapacities_F.mat']
        exp_result.sort()
        self.assertEqual(exp_result, self.simdex_filtered_fn,
                         'filtering self.simdex with c1.C = '' should return\
                         all LinkedCapacities*')
        self.assertEqual(filt_dic, self.simdex_filtered.filterset,
                         'After filtering, filterset has to be updated')
        self.simdex.h5.close()
    
    def test_filter_twice(self):
        """ filtering twice should give correct end results and filterset"""

        filt_dic = {'c1.C': 800}
        filt_dic2 = {'r.R': 3}
        self.simdex_filtered1 = self.simdex.filter(filt_dic)
        self.simdex_filtered = self.simdex_filtered1.filter(filt_dic2)
        self.simdex_filtered_fn = self.simdex_filtered.get_filenames()
        self.simdex_filtered_fn.sort()
        exp_result = ['LinkedCapacities_A.mat',  'LinkedCapacities_F.mat']
        exp_result.sort()
        self.assertEqual(exp_result, self.simdex_filtered_fn,
                         "filtering self.simdex with c1.C = 800 then with r.R = 3 \
                         should return LinkedCapacities_A and _F")
        self.assertEqual({'c1.C': 800, 'r.R': 3}, self.simdex_filtered.filterset,
                         'After filtering, filterset has to be updated')
        self.simdex.h5.close()
        
    def test_filter_unchanged_original(self):
        """self.simdex.filter() should not change self.simdex (bug and issue on github)"""

        filt_dic = {'c1.C': 800}
        self.simdex_filtered = self.simdex.filter(filt_dic)
        
        self.assertEqual({}, self.simdex.filterset,
                        'Simdex.filter() should NOT change filterset of Simdex')
        self.simdex.h5.close()

    def test_filter_floatvalues(self):
        """Simdex.filter() with float values should work well"""
        

        filt_dic = {'r.R': 8.15}
        self.simdex_filtered = self.simdex.filter(filt_dic)
        self.simdex_filtered_fn = self.simdex_filtered.get_filenames()
        self.simdex_filtered_fn.sort()
        exp_result = ['LinkedCapacities_E.mat']
        exp_result.sort()
        self.assertEqual(exp_result, self.simdex_filtered_fn,
                         'filtering self.simdex with r.R = 8.15 should return\
                         exp_result')
        self.assertEqual(filt_dic, self.simdex_filtered.filterset,
                         'After filtering, filterset has to be updated')
        self.simdex.h5.close()
    
    def test_plot(self):
        """Simdex.plot() should return [fig, lines, leg]"""
        
        self.simdex.h5.close()
        vardic = {'T2': 'c2.T'}
        pardic = {'parc': 'c1.C'}
        process=Process(parameters=pardic, variables=vardic)
        self.simdex = Simdex(folder=getcwd(), process=process)        
        
        self.simdex_filtered = self.simdex.filter({'c1.C': ''})
        [fig, lines, leg] = self.simdex_filtered.plot('T2')
        self.assertTrue(isinstance(fig, matplotlib.figure.Figure))
        self.assertEqual(7, len(lines))
        for line in lines:
            self.assertTrue(isinstance(line, matplotlib.lines.Line2D))
        self.assertTrue(isinstance(leg, matplotlib.legend.Legend))
        self.simdex.h5.close()
        
        
    def test_save_and_load(self):
        """Saving and loading a self.simdex object should return exactly the same object"""
        
        self.simdex.save('Test_save.dat')
        loaded = load_simdex('Test_save.dat')
        for attr in self.simdex.__dict__:
            exec("s = self.simdex." + attr)
            exec("l = loaded." + attr)
                   
            if isinstance(s, np.ndarray):
                self.assertTrue((l == s).all())
            else: 
                self.assertEqual(s, l)

    
    def test_scatterplot(self):
        """Simdex.scatterplot() should return [fig, lines, leg]"""
        
        self.simdex.h5.close()
        vardic = {'T1':'c1.T', 'T2': 'c2.T'}
        pardic = {'parc': 'c1.C'}
        process=Process(parameters=pardic, variables=vardic)
        self.simdex = Simdex(folder=getcwd(), process=process) 
        
        self.simdex_filtered = self.simdex.filter({'c1.C': ''})
        [fig, lines, leg] = self.simdex_filtered.scatterplot('T1', 'T2')
        self.assertTrue(isinstance(fig, matplotlib.figure.Figure))
        self.assertEqual(7, len(lines))
        for line in lines:
            self.assertTrue(isinstance(line, matplotlib.lines.Line2D))
        self.assertTrue(isinstance(leg, matplotlib.legend.Legend))  
        self.simdex.h5.close()
        
    def test_get_fullnames(self):
        """Simdex.get() should return correctly for variables and parameters"""
        
#        result = self.simdex.get('c2.T')
#        expectedsids = self.simdex.get_SID('linkedcapacities')
#        self.assertEqual(sorted(result.keys()), sorted(expectedsids))
#        sim = Simulation('LinkedCapacities_B.mat')
#        sid = self.simdex.get_SID('_B')[0]
#        self.assertTrue(np.all(result[sid]==sim.get_value('c2.T')))

        # now a parameter
        c1_C = self.simdex.get('c1.C').values()
        c1_C.sort()
        print '\n'*10, 'c1_C = ', c1_C
        exp_result_sorted = np.array([600.,   600.,   800.,  
                                      800.,   800.,   800., 1000., np.NaN])
        print c1_C
        print exp_result_sorted
        self.assertTrue(((c1_C == exp_result_sorted) | (np.isnan(c1_C) & np.isnan(exp_result_sorted))).all())
        self.simdex.h5.close()
        
        
    def test_get_shortnames(self):
        """Simdex.get() should return correctly for shortnames too"""
        
        self.simdex.h5.close()
        vardic = {'T2': 'c2.T'}
        pardic = {'parc': 'c1.C'}
        process=Process(parameters=pardic, variables=vardic)
        self.simdex = Simdex(folder=getcwd(), process=process)        
        result = self.simdex.get('T2')
        expectedsids = self.simdex.get_SID('linkedcapacities')
        self.assertEqual(sorted(result.val.keys()), sorted(expectedsids))
        sim = Simulation('LinkedCapacities_B.mat')
        sid = self.simdex.get_SID('LinkedCapacities_B')[0]
        self.assertTrue(np.all(result.val[sid]==sim.get_value('c2.T')))
        
        c1_C = self.simdex.get('parc').values()
        c1_C.sort()
        exp_result_sorted = np.array([600.,   600.,   800.,  
                                      800.,   800.,   800., 1000., np.NaN])
        #self.assertTrue(((c1_C == exp_result_sorted) | (np.isnan(c1_C) & np.isnan(exp_result_sorted))).all())
        np.testing.assert_equal(c1_C , exp_result_sorted)
        self.simdex.h5.close()        
        
    def test_get_sub_var(self):
        """Simdex.get() should return correctly for aggregated subvariables"""
        
        self.simdex.h5.close()
        
        mothers=['c1', 'c2']
        parameters={'cap1':'c1.C', 'res':'r.R'}
        sub_pars={'cap':'C'}
        sub_vars={'Qflow':'heatPort.Q_flow'}        
        process=Process(parameters=parameters, sub_vars = sub_vars, 
                        mothers=mothers, sub_pars=sub_pars)

        self.simdex = Simdex(folder=getcwd(), process=process)        

        Q1 = self.simdex.get('c1_Qflow').val
        Q2 = self.simdex.get('c2_Qflow').val
        Q = self.simdex.get('Qflow').val
        
        for sid in Q1.keys():        
            self.assertTrue(np.all(Q1[sid] == Q[sid][:,0]))
            self.assertTrue(np.all(Q2[sid] == Q[sid][:,1]))

    def test_get_sub_var2(self):
        """Simdex.get() should return correctly for sum and mean-aggregated subvariables"""
        
        self.simdex.h5.close()
        
        mothers=['c1', 'c2']
        parameters={'cap1':'c1.C', 'res':'r.R'}
        sub_pars={'cap':'C'}
        sub_vars={'Qflow':'heatPort.Q_flow'}        
        process=Process(parameters=parameters, sub_vars = sub_vars, 
                        mothers=mothers, sub_pars=sub_pars)

        self.simdex = Simdex(folder=getcwd(), process=process)        

        Q1 = self.simdex.get('c1_Qflow').val
        Q2 = self.simdex.get('c2_Qflow').val
        QSum = self.simdex.get('Qflow', aggregate='sum').val
        QMean = self.simdex.get('Qflow', aggregate='mean').val
        
        for sid in Q1.keys():        
            self.assertTrue(np.all((Q1[sid]+Q2[sid])/2 == QMean[sid]))
            self.assertTrue(np.all((Q1[sid]+Q2[sid]) == QSum[sid]))

        self.simdex.h5.close()                
#if __name__ == '__main__':
#    unittest.main()


suite1 = unittest.TestLoader().loadTestsFromTestCase(ProcessTest)
suite2 = unittest.TestLoader().loadTestsFromTestCase(ResultTest)
suite3 = unittest.TestLoader().loadTestsFromTestCase(SimulationTest)
suite4 = unittest.TestLoader().loadTestsFromTestCase(SimdexTest)


alltests = unittest.TestSuite([suite1, suite2, suite3, suite4])

unittest.TextTestRunner(verbosity=0, failfast=True).run(alltests)
#unittest.TextTestRunner(verbosity=1).run(suite2)

