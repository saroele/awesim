# -*- coding: utf-8 -*-
"""
Created on Fri Jan 07 11:10:01 2011

@author: RDC
"""

import unittest

class MyClass:
    def my_function(self):
        raise IOError
    
class my_function_tester(unittest.TestCase):
    def test_my_function(self):
        self.assertRaises(IOError, MyClass.my_function(), 'should raise IOError')
        
if __name__ == '__main__':
    unittest.main()
    