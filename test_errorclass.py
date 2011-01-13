# -*- coding: utf-8 -*-
"""
Created on Fri Jan 07 11:17:20 2011

@author: RDC
"""
import unittest


class MyClass:
    def __init__(self, value):
        self.value = value
    
    def my_method(self):
        if self.value < 0 :
            result = self.value
        else:
            result = False
            raise myError
        return result

class myError(Exception):
    "this is the error text"

class my_method_tester(unittest.TestCase):
    def test_my_function_neg(self):
        self.myclass = MyClass(-3)
        self.assertEqual(-3, self.myclass.my_method(), 
                         'with negative values myclass.my_method() \
                          should give input value as result')
    
    def test_my_method_pos_workaround(self):
        self.myclass = MyClass(4)
        error_occurred = False
        try:
            self.myclass.my_method()
        except myError:
            error_occurred = True
        self.assertTrue(error_occurred,
                        'with positive values myclass.my_method() \
                        should raise myError')
    
    def test_my_method_pos(self):
        self.myclass = MyClass(4)
        self.assertRaises(myError, self.myclass.my_method)
        
if __name__ == '__main__':
    unittest.main()
    