# -*- coding: utf-8 -*-
"""
Test module for the simman module

Created on Thu Jan 06 20:46:26 2011

@author: RDC
"""

import numpy as np
import unittest
from os import getcwd, path
from cStringIO import StringIO
import sys
from simman import Simulation, Simdex

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
            
    def test___init__bigfile(self):
        """
        Tests creation of Simulation object from well known file.
        Checks if attributes are present 
        Make sure the file Array_Big.mat is in the current work directory

        """
        cwd = getcwd()
        filename = path.join(cwd, 'SubfolderForArray_Big', 'Array_Big')        
        try:
            sim = Simulation(filename)
            for attr in ['dataInfo', 'data_1', 'data_2', 'filename', 'names']:
                self.assertTrue(sim.__dict__.has_key(attr), 
                            'Simulation.__init__() dit not create \
                            attribute %s' % attr)
        except IOError:
            print """
            
        !!!!!!!!!!!!!!!!!!!!  PAY ATTENTION !!!!!!!!!!!!!!!!!!!!
        
        Check if Array_Big.mat' is present in folder 
        'SubfolderForArray_Big' of your work directory"
        
        If yes, this test did NOT succeed, even if you get OK!
        """
        
    
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
        self.assertEqual([u'Time', u'c1.heatPort.T', u'c1.heatPort.Q_flow',
                          u'c1.T', u'c1.der(T)', u'c2.heatPort.T',
                          u'c2.heatPort.Q_flow', u'c2.T', u'c2.der(T)',
                          u'r.heatPort_a.T',  u'r.heatPort_a.Q_flow',
                          u'r.heatPort_b.T', u'r.heatPort_b.Q_flow'], 
                          sim.variables) 
                          
    def test_separate_twice(self):
        """ Tests if all goes fine if Simulation.separate() is called twice """
        
        sim = Simulation('LinkedCapacities')    
        sim.separate()
        sim.separate()
        
        parameters = np.array([600.0, 1000.0, 3])
        

        self.assertEqual([u'c1.C', u'c2.C', u'r.R'], sim.parameters)
        self.assertTrue((parameters == sim.parametervalues).all())
        self.assertEqual([u'Time', u'c1.heatPort.T', u'c1.heatPort.Q_flow',
                          u'c1.T', u'c1.der(T)', u'c2.heatPort.T',
                          u'c2.heatPort.Q_flow', u'c2.T', u'c2.der(T)',
                          u'r.heatPort_a.T',  u'r.heatPort_a.Q_flow',
                          u'r.heatPort_b.T', u'r.heatPort_b.Q_flow'], 
                          sim.variables)

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
        f=StringIO("y\n" * N)
        sys.stdin = f
        
        self.cwd = getcwd()
        # sims contains the simualations we expect to be in simdex in the 
        # current work directory        
        self.sims = ['']
        self.filenames = ['Array.mat', 'LinkedCapacities.mat', \
                   'LinkedCapacities_A.mat',  'LinkedCapacities_B.mat', \
                   'LinkedCapacities_C.mat', 'LinkedCapacities_D.mat', \
                   'LinkedCapacities_E.mat', 'LinkedCapacities_F.mat']
        self.filenames.sort()
        for fn in self.filenames:
            self.sims.append(path.join(self.cwd, fn))
                       
    
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
        simdex = Simdex()
        filenames = simdex.get_filenames('path')
        filenames.sort()
        self.assertEqual(self.sims[1:], filenames)
        
        
    def test_init_subfolder(self):
       """ Test initiation from a different folder with absolute pathname"""
       
       folder = path.join(self.cwd, 'Subfolder')
       simdex = Simdex(folder)
       filenames = simdex.get_filenames()
       filenames.sort()
       self.assertEqual(self.filenames, filenames)
       
    def test_init_subfolder_with_spaces(self):
       """ Test initiation from a folder with spaces with absolute pathname"""
       
       folder = path.join(self.cwd, 'A Subfolder with Spaces')
       simdex = Simdex(folder)
       filenames = simdex.get_filenames()
       filenames.sort()
       self.assertEqual(self.filenames, filenames)
       
    def test_init_array_big(self):
       """ Test initiation with a large .mat file"""
       
       folder = path.join(self.cwd, 'SubfolderForArray_Big')
       
       try:
            simdex = Simdex(folder)
            self.assertEqual(path.join(folder, 'Array_Big.mat'), 
                             simdex.get_filenames('path')[0])
       except IOError:
            print """
            
        !!!!!!!!!!!!!!!!!!!!  PAY ATTENTION !!!!!!!!!!!!!!!!!!!!
        
        Check if Array_Big.mat' is present in folder 
        'SubfolderForArray_Big' of your work directory"
        
        If yes, this test did NOT succeed, even if you get OK!
        """
            
    def test_init_subfolder_with_crappy_files(self):
       """ Test initiation from a folder including wrong .mat files"""
       
       folder = path.join(self.cwd, 'SubfolderWithCrappyFiles')
       simdex = Simdex(folder)
       filenames = simdex.get_filenames()
       filenames.sort()
       self.assertEqual(self.filenames, filenames)
       
    def test_exist(self):
        """
        Test if :
        1) running Simdex.exist() gives a match when it should
           (different cases are tested)
        2) running Simdex.exist() gives a zero list for non-present regex
        
        """
        
        simdex = Simdex()
        self.assertEqual(simdex.parameters[-1], 
                         simdex.exist(simdex.parameters[-1])[0][0], 
                         'Simdex.exist() does NOT return a present name')
        self.assertEqual(1, len(simdex.exist('c1')[0]), 
                         "sim.exist('c1') should return 1 parameter")
        self.assertEqual(1, len(simdex.exist('c1', 'par')), 
                         "sim.exist('c1', 'par') should return 1 parameter")                 
        self.assertEqual(4, len(simdex.exist('c1')[1]), 
                         "sim.exist('c1') should return 4 variables")
        self.assertEqual(4, len(simdex.exist('c1', 'var')), 
                         "sim.exist('c1', 'var') should return 4 variables")                 
        self.assertEqual(u'Time', simdex.exist('time')[1][0], 
                         'Simulation.exist() should ignore case' )
        self.assertEqual([[],[]], simdex.exist('this does not exist'))
        
    def test_cleanup(self):
        """
        Test cleanup for pars and vars separately
        I had trouble with checking if the right par and right var is removed
        therefore I only control the sizes of attributes and suppose the 
        right one is gone
                
        """
        simdex = Simdex()
        n = 4
        # set 1 value in the parametermap to 0
        par = simdex.parameters[n]
        npars = len(simdex.parameters)
        nvars = len(simdex.variables)
        simdex.parametermap[n, 1] = 0
        simdex.cleanup()
        self.assertEqual(npars-1, len(simdex.parameters))
        self.assertEqual(npars-1, simdex.parametermap.shape[0])
        self.assertEqual(npars-1, simdex.parametervalues.shape[0])
        self.assertEqual(nvars, len(simdex.variables))
        simdex.variablemap[n, 1] = 0
        simdex.cleanup()
        self.assertEqual(nvars-1, len(simdex.variables))
        self.assertEqual(nvars-1, simdex.variablemap.shape[0])
        self.assertEqual(npars-1, simdex.parametervalues.shape[0])
    
    def test_get_simID(self):
        """
        check if Simdex.get_simID() returns correct simID
        maybe this is not really an important check
        
        """
        

    
    
    def test_filter(self):
        """
        dummy docstring
        """
        pass
       
       
#if __name__ == '__main__':
#    unittest.main()

suite1 = unittest.TestLoader().loadTestsFromTestCase(SimulationTest)
suite2 = unittest.TestLoader().loadTestsFromTestCase(SimdexTest)
alltests = unittest.TestSuite([suite1, suite2])

unittest.TextTestRunner(verbosity=1).run(alltests)
#unittest.TextTestRunner(verbosity=1).run(suite2)

# Restore back interactivity of keyboard
sys.stdin = sys.__stdin__