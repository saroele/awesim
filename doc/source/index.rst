.. SimulationManagement documentation master file, created by
   sphinx-quickstart on Wed Aug 29 15:35:15 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Simulation Management
=====================

*Optimize your simulation workflow.*  

Simulation Management is written in **Python** and contains tools for pre- and post-processing of **Modelica** models.      

With Simulation Management you can:

    1. compile models, set parameters and solver options and start simulations 
    2. run simulations in parallel on the same computer or different machines
    3. open result files, access simulation parameters and trajectories, make plots
    4. keep an overview of your cases by creating an index (called simdex) of different results
    5. filter your results based on filenames or parameter sets 
      
Most of the tools provided are currently based on Dymola and the .mat result file. 

Simulation management is licensed under the xxx license.

See the documentation for more information.



Simulation management should be renamed to simman I think, and I should then have a separate file for Simdex, Simulation, Process and Result, all in the folder simman with an __init__.py in there. new test for sourcecode

Hier is een eerste test::

    def myfunc(foo, bar=True):
        print foo

Je kan ook gewoon code schrijven, maar er moet een lege lijn tussen::

    a = range(10)

Ofwel met sourcecode ipython, tweede test

.. ipython:: python

   x = 2
   x**3
   x**200
   

.. literalinclude:: DemoProcess.py   
   
   
Contents:

.. toctree::
   :maxdepth: 2

Class Simulation
================    
.. autoclass:: simman.Simulation
   :members:
   :private-members:
   :special-members:
   
Class Process
================    
.. autoclass:: simman.Process
   :members:
   :private-members:
   :special-members:
   
Class Simdex
================    
.. autoclass:: simman.Simdex
   :members:
   :private-members:
   :special-members:
   
Class Result
================    
.. autoclass:: simman.Result
   :members:
   :private-members:
   :special-members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

