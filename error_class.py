# -*- coding: utf-8 -*-
"""
Created on Fri Jan 07 11:10:01 2011

@author: RDC
"""


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