# -*- coding: utf-8 -*-
"""
Created on Tue Mar 08 10:49:08 2011

@author: RBA
"""

from buiuser import *
import matplotlib.pyplot as plt

to_time=[]
to_data=[]
a=np.array([0,10,25,30,48,21,81,79,80,91,100])
b=np.array([0,2,4,6,8,10,12,14,16,18,20])
to_time, to_data = get_smoothened_data(a,b,4,'tryout')
print(to_time,to_data)

plt.plot(b,a,to_time,to_data)
plt.show()
