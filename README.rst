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

Design a Multistate Design
~~~~~~~~~~~~~~~~~~~~~~~~~~

This simple script generates a multistate design, which is a RNA molecule that is able to fold into
all the given structural conformation. In case of a bistable molecule just call:

.. code:: bash

    echo -e '(((((....)))))....\n....(((((....)))))' | multistate-design.py -i -m random -e 500

The program `barriers` can the be used to visualize the energy landscape to confirm the design goals:

.. code:: text

      GUGACCGCGGUCACGUGG
    1 (((((....)))))....  -7.00    0  10.00
    2 ....(((((....)))))  -7.00    1   9.50
    3 ..................   0.00    2   1.60
    4 ....(((......)).).   0.80    1   2.00
    5 ((....))..........   1.10    1   1.50
    6 .......((......)).   1.50    1   1.10
    7 ...((.........))..   2.40    1   0.40
    8 .((....)).........   2.60    1   0.20

.. figure:: data/barriers.png
    :width: 350px

    Barriers Tree showing the two desired states as deep minima (1, 2) and the open chain (3) as neighbouring
    minimum.

Design a sRNA mediated translational regulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:code:`TODO`: `cofold_design.py`

Design a simple Thermoswitch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    echo -e "((((((((((....)))))))))) 24.0\n((((....))))((((....)))) 37.0\n((((....))))............ 46.0" | multistate-thermoswitch.py -m random -e 1000

This results e.g in a sequence like `GGGUUGAUACCCGAGUGUUGAUUC` which has the given MFE structures at the specified temperatures.
Folding it at all Temeratures from 10 to 100 degree celsius shows, that the first structural change happens at 30.0 degree
celsius and the second one at 45 degrees. After 85 degrees, the sequence occurs only in the open chain conformation.

`RNAheat` further confirms that the designed sequence is indeed a three-stable thermoswitch:

.. figure:: data/RNAheat.png
    :width: 350px
    
    RNAheat Plot, ViennaRNA v2.2.9, GGGUUGAUACCCGAGUGUUGAUUC

Display the Dependency Graph
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We generated two example scripts which can dump the Dependency Graph in the common GraphML format and,
by using the :code:`igraph` python library, render these files as images.
Following example input is possible:

.. code:: bash
   
    echo -e '(((((....)))))....\n....(((((....)))))' | generate_graphml.py -i > dependency_graph.gml
    print_graphml.py -g dependency_graph.gml -o dependency_graph.png

Or use the second script directly:

.. code:: bash
    
    echo -e '(((((....)))))....\n....(((((....)))))\n(((((((....)))))))' | print_graphml.py -i

This results in a nice representation of the dependency graph:

.. figure:: data/graph.png
    :width: 350px
    
    Very simple dependency graph visualized using igraph.
