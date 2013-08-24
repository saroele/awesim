.. currentmodule:: awesim
.. _process:

Defining a post-process
=======================

Post-processing is the data treatment that is unavoidable after running your simulations.  There are different things to get done:

* extract all useful information from the simulation file(s)
* treat the results: number crunching like scaling, integration, aggregation, or simple arithmetic operations
* compare results of different runs with plots, tables, etc. 
* store your results for later use

Post processing is often a time-intensive step in the simulation set-up.  To ease often used post-processing operations, awsim defines a Process class.  A Process is a set of operations that can be executed on a Simulation instance.  Let's start with a simple example.

.. ipython:: python
    :suppress:
    
    import os
    print os.getcwd()
    os.chdir('../tests')
    print os.getcwd()
    
.. ipython:: python
    
    import numpy as np
    from awesim import Simulation, Process
    
.. ipython:: python
    :suppress:
    
    np.set_printoptions(precision=2, suppress=True)
    
.. ipython:: python
    
    process = Process(variables={'T': 'c1.T'})
    print process

This process has a dictionary ``variables`` that contains a mapping of short to full variable names. You'll notice that the variable 'Time' is automatically added to the variables.  When we apply this simple process to a Simulation, the result will be the same as invoking ``simulation.extract(process.variables)``
    
.. ipython:: python

    sim = Simulation('LinkedCapacities')
    sim.postprocess(process)
    
If we need parameter values in the postprocessing, we can add (or edit) the attribute ``parameters``.

.. ipython:: python
    
    process.parameters = {'c1': 'c1.C', 'c2': 'c2.C'}
    print process
    
One of the main uses of the process class is the definition of post-processing actions.  These are defined as strings.  A simple example could be the conversion of the temperature of capacity c1 from Kelvin to degree Celsius. 

.. ipython:: python
    
    post_proc_string = 'T_degC = T - 273.15'
    process.pp.append(post_proc_string)
    sim.postprocess(process)
    
You'll notice that we have created a new variable, T_degC which is added to the result of the post-processing.  Note that in the post-processing string, we were able to use the shortname 'T' that was introduced as key in the ``variables`` of the process.  
