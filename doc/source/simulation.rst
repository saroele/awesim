.. currentmodule:: awesim
.. _simulation:

Working with single simulations
===============================

To get started, we'll always use the following import statements

.. ipython:: python
    :suppress:
    
    import os
    print os.getcwd()
    os.chdir('../tests')
    print os.getcwd()
    
.. ipython:: python
    
    import numpy as np
    from awesim import Simulation
    
.. ipython:: python
    :suppress:
    
    np.set_printoptions(precision=2, suppress=True)
    
Suppose we have a simulation result file called 'LinkedCapacities.mat'.
We can instantiate a :class:'~awesim.Simulation' object from this file.

.. ipython:: python

    sim = Simulation('LinkedCapacities') #with our without .mat extension
    print sim

When instantiating a Simulation object, the .mat file is read in memory. 
The filename is stored as an attribute.
    
.. ipython:: python
    
    sim.filename

A simulation object can be introspected in different ways.  A list of all
known parameters and variables can be obtained.  For large simulations, the method ``get_objects()`` will be more practical.  It lists all the sub-objects from a given parent.  The simulation root will be taken if no parent is given.
.. ipython:: python

    sim.names
    sim.get_objects()
    sim.get_objects('C1')
    
The list of known names can be split in parameters and variables with the method ``separate()``.
This will create three attributes:

1. ``Simulation.variables``
2. ``Simulation.parameters``
3. ``Simulation.parametervalues``

.. ipython:: python

    sim.separate()
    sim.variables
    for p,v in zip(sim.parameters, sim.parametervalues):
        print p, ' = ', str(v)

A search method is foreseen: ``Simulation.exist()``.  It will return a list of all names that satisfy the search criterium.  You can even use regular expressions in your search.

.. ipython:: python

    sim.exist('q_flow')

There are two methods to obtain the values.  The first one is ``Simulation.get_value()``.  This method will return a numpy.float for paramters and a numpy.array for trajectories (variables). 

    Note: the simulation time is always accessible through 'Time'

..  ipython:: python
    
    time = sim.get_value('Time')
    Q = sim.get_value(u'r.heatPort_a.Q_flow')
    Q_sum = np.trapz(Q, time)
    print "The total energy that flowed through the resistance is %.1f J" %Q_sum

The second method is ``Simulation.extract()``::
    
    def extract(self, var, arrays='sum'):
        """
        Return dictionary with values of the variables/parameters to extract.
        """
        
extract() takes a dictionary as input argument, and will return a dictionary of the same length with the same keys.  The variable or parameter names will be replaced with their values.  

.. ipython:: python

    sim.extract({'c1': 'c1.C', 'T1': 'c1.T', 'Q1':'c1.heatPort.Q_flow'})
    
If you want to extract an array of variables, just replace the index (between the []) by 'x'.
There are three possible options for the processing of arrays:

#. arrays='each': return each of the arrays
#. arrays='sum': return the sum (default)
#. arrays='mean': return the mean

Here's an example of a simulation result file that contains arrays.
  
.. ipython:: python

    sim = Simulation('Array.mat')
    sim.get_objects()
    sim.extract({'T_array': 'c[x].T'}, arrays='each')
    sim.extract({'T_array': 'c[x].T'}, arrays='mean')
    
There is one more method in the Simulation class: ``postprocess()``.  This method is related to the Result class, it will be explained in the documentation of Result.
   
At this moment, there is no automatic plotting method in Simulation.  The reason is that is is more convenient to work with the Simdex class, which has methods for plotting.  See the documentation of Simdex.

