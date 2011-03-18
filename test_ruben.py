# -*- coding: utf-8 -*-
"""
Created on Tue Mar 08 10:49:08 2011

@author: RBA
"""

import buiuser
import numpy as np
import unittest
# from os import getcwd, path
# from cStringIO import StringIO
# import sys
import matplotlib.pyplot as plt
# from simman import Simulation, Simdex, load_simdex

class BuiuserTest(unittest.TestCase):
    """
    Class for testing buiuser.py
    """
    def test___init__smoothendata(self):
        """
        Tests the function for smoothening data from well known file.
        """
        data = np.array([0, 10, 25, 30, 48, 21, 81, 79, 80, 91, 100])
        time = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
        to_time, to_data = buiuser.get_smoothened_data(data, time, 4, 'tryout')
        self.assertEqual(len(to_time),len(to_data),'smoothening failed')

suite1 = unittest.TestLoader().loadTestsFromTestCase(BuiuserTest)
# suite2 = unittest.TestLoader().loadTestsFromTestCase(SimdexTest)
alltests = unittest.TestSuite([suite1])
# alltests = unittest.TestSuite([suite1, suite2])

unittest.TextTestRunner(verbosity=1).run(alltests)
            