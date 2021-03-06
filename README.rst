RNAsketch Library for designing RNA molecules
=============================================

Installation
------------

Download the lastes RNAsketch package from https://github.com/ViennaRNA/RNAsketch/releases/latest, exctract and simply use setuptools to install the module and the few sample scripts::

    python setup.py install

I highly recommend to install the module locally into a directory specified in your
:code:`PYTHONPATH` variable. This can be achieved either by using the :code:`--prefix=$HOME/local_path` 
argument, or in case :code:`PYTHONUSERBASE` is set, just use the :code:`--user` flag.

Build the documentation
-----------------------

Documentation available online: http://viennarna.github.io/RNAsketch

Call these commands to build the html documentation::

    cd doc
    make html

Cite
----

Stefan Hammer, Birgit Tschiatschek, Christoph Flamm, Ivo L. Hofacker, and Sven Findeiß. “RNAblueprint: Flexible Multiple Target Nucleic Acid Sequence Design.” Bioinformatics, 2017. https://doi.org/10.1093/bioinformatics/btx263.

If you used the multi-dimensional boltzmann sampling scripts `design-redprint-multistate.py` or `design-energyshift.py` please cite:

Stefan Hammer, Yann Ponty, Wei Wang, and Sebastian Will. “Fixed-Parameter Tractable Sampling for RNA Design with Multiple Target Structures.” In RECOMB 2018. Paris, France, 2018. https://arxiv.org/abs/1804.00841.


Provided Example Scripts
------------------------

Design a Multistate Design
~~~~~~~~~~~~~~~~~~~~~~~~~~

This simple script generates a multistate design, which is a RNA molecule that is able to fold into
all the given structural conformation. In case of a bistable molecule just call:

.. code:: bash

    echo -e '(((((....)))))....\n....(((((....)))))' | design-multistate.py -i -m random -s 500

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
.. figure:: doc/data/barriers.png
    :width: 350px

    Barriers Tree showing the two desired states as deep minima (1, 2) and the open chain (3) as neighbouring
    minimum.

Design a sRNA mediated translational regulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
    
    echo '...(((((((((((((((((&)))))))))))))))))((((....))));...(((((((((((((((((&)))))))))))))))))............
    ....................&.................((((....))));...xxxxxxxxxxxxxxxxx&xxxxxxxxxxxxxxxxx............
    NNNNAAGGAGNNNNNNNAUG&NNNNNNNNNNNNNNNNNNNNNNNNNNNNN' | design-cofold.py -n 1 -s 1000

This small example will design a simple device consisting of a 5'UTR region which can
be translationally controlled by a sRNA molecule. In this case the sRNA will shut down
translation by directly binding the RBS (Ribosome Binding Site) and the AUG start codon.

`RNAcofold -a -p -d2` calculates three dot-plots showing the base pair probabilities in the ensemble of states which
confirms the design objective:

.. figure:: data/cofold.png
    :width: 350px
.. figure:: doc/data/cofold.png
    :width: 350px
    
    RNAcofold Dot-Plots, ViennaRNA v2.2.9, AAAUAAGGAGUAAAUGAAUG&CAUUCAUUUACUCCUUACCGCACUCGCGG
    Plots were assembled in a single picture for better comparison. Only base pair probabilities
    are shown in the plots.
    
    Score: 0.89; complex concentration: 1.00; P(5UTR unpaired): 0.97; P(sRNA unpaired): 0.99; P(mRNA context): 0.18

Design a multistate Thermoswitch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    echo -e "(((((((((((((....))))))))))))) 5.0\n(((((.....)))))(((((.....))))) 10.0\n(((((.....)))))............... 37.0" | design-thermoswitch.py -m random -s 1000

This results e.g in a sequence like `GAUCUGUGUGGGGUCGAUUUUGUGUGGGUU` which has the given MFE structures at the specified temperatures (lower plot).
Folding it at all Temeratures from 10 to 100 degree Celsius shows, that the first structural change happens at ~7.0 degree
Celsius and the second one at ~26 degrees. After _72 degrees, the sequence occurs only in the open chain conformation.

`RNAheat` further confirms that the designed sequence is indeed a three-stable thermoswitch:

.. figure:: data/thermoswitch.png
    :width: 350px
.. figure:: doc/data/thermoswitch.png
    :width: 350px
    
    RNAheat Plot, ViennaRNA v2.2.9, GAUCUGUGUGGGGUCGAUUUUGUGUGGGUU

Design a ligand triggered switch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
    
    echo -e "(((((...((((((((.....)))))...)))...)))))........................\n.........................(((((((((((......)))))))))))...........\nAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCANNNNNNNNNNNNNNNNNNNNNN" | design-ligandswitch.py -r 70:30 --ligand "GAUACCAG&CCCUUGGCAGC;(...((((&)...)))...);-9.22"

This designs a simple theophylline triggered switch. It adapts to a certain ratio of 
the aptamer structure (ligand competent state) and an alternative state as specified
with the `--ratio` option. To model ligand binding we use the soft-constraints framework
of the ViennaRNA package, similar to the `--motif` option of `RNAfold`.
The specified objective function calculates the probabilities of the structural features with
and without ligand. Thus, we optimize towards the given ratio without the ligand and maximize the
ligand binding competent state in the presence of the ligand.

.. figure:: data/ligandswitch.png
    :width: 350px
.. figure:: doc/data/ligandswitch.png
    :width: 350px
    
Design a multistate riboswitch with equal target energies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
    
    export REDPRINT=</path/to/redprint/>
    echo -e "(((((((((((((....)))))))))))))\n(((((.....)))))(((((.....)))))\n(((((.....)))))..............." | design-redprint-multistate.py -i -n 10

Using this script you can design a multistate riboswitch similar to the `design-multistate.py` script with fundamental differences in the sampling technology.
Here, we use `RNARedPrint` to gain Boltzmann-weighted sampling of sequences given the structural inputs. This script then implements a strategy to achieve multi-dimensional boltzmann sampling with the target energy being the mean achieveable energy of the target structures.

The resulting sequences exhibit equal energies for the given target structurs. This does not mean that the target structures are highly populated in the ensemble. To additionally achieve this, specify the `-s 200` option to subsequently do a optimization towards the multi-target objective function.

Either put the complete `RNARedPrint` project folder next to the scripts, or set a `REDPRINT` environmental variable before executing this script!

Design a multistate riboswitch with specific target energies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
    
    export REDPRINT=</path/to/redprint/>
    echo -e ".((((((......)))))).((((...((((((...((((...(((.......)))..........(((....))).))))..))))))...))))....\n.((((((......)))))).((((...((((((...((((((.(((.......((........))..)))....)).))))..))))))...))))....\n......((.((((.(((((((.((.((...((((..((.....(((.......))).....))..)))).)).)))))))))..)).)).))........" | design-energyshift.py -i -e 40,40,20

Design a multistate riboswitch with specific target energies and GC-content. This script uses RNARedPrint for Boltzmann-weighted sampling and implements a strategy to achieve multi-dimensional boltzmann sampling.

Either put the complete `RNARedPrint` project folder next to the scripts, or set a `REDPRINT` environmental variable before executing this script!


The resulting sequences are nicely distributed around the specified target energies:

.. figure:: data/energyshift.png
    :width: 350px
.. figure:: doc/data/energyshift.png
    :width: 350px


Display the Dependency Graph
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We generated two example scripts which can dump the Dependency Graph in the common GraphML format and,
by using the :code:`igraph` python library, render these files as images.
Following example input is possible:

.. code:: bash
   
    echo -e '(((((....)))))....\n....(((((....)))))' | design-generategraphml.py -i > dependency-graph.gml
    design-printgraphml.py -g dependency-graph.gml -o dependency-graph.png

Or use the second script directly:

.. code:: bash
    
    echo -e '(((((....)))))....\n....(((((....)))))\n(((((((....)))))))' | design-printgraphml.py -i

This results in a nice representation of the dependency graph:

.. figure:: data/graph.png
    :width: 350px
.. figure:: doc/data/graph.png
    :width: 350px
    
    Very simple dependency graph visualized using igraph.
