PyDesign Library for designing RNA molecules
============================================

Installation
------------

Simple use setuptools to install the module and the few sample scripts::

    python setup.py install

I highly recommend to install the module locally into a direcory specified in your
:code:`PYTHONPATH` variable. This can be achieved either by using the :code:`--prefix=$HOME/local_path` 
argument, or in case :code:`PYTHONUSERBASE` is set, just use the :code:`--user` flag.

Build the documentation
-----------------------

Documentation available online: http://ribonets.github.io/PyDesign

Call these commands to build the html documentation::

    cd doc
    make html

Provided Example Scripts
------------------------

Design a Multistate Riboswitch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    echo -e '(((((....)))))....\n....(((((....)))))' | multistate-riboswitch.py -i -m random -e 500


Design a sRNA mediated translational regulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cofold_design.py  

Design a simple Thermoswitch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    echo -e "((((((((((....)))))))))) 24.0\n((((....))))((((....)))) 37.0\n((((....))))............ 46.0" | multistate-thermoswitch.py -m random -e 1000

This results e.g in a sequence like `GGGUUGAUACCCGAGUGUUGAUUC` which has the given MFE structures at the specified temperatures.
Folding it at all Temeratures from 10 to 100 degree celsius shows, that the first structural change happens at 30.0 degree
celsius and the second one at 45 degrees. After 85 degrees, the sequence occurs only in the open chain conformation.

`RNAheat` further confirms that the designed sequence is indeed a three-stable thermoswitch:

.. image:: doc/data/RNAheat.png?raw=true 

RNAheat Plot, ViennaRNA v2.2.9, GGGUUGAUACCCGAGUGUUGAUUC

Display the Dependency Graph
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_graphml.py  print_graphml.py
