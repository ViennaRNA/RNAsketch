#### PyDesign Library for designing RNA molecules ####

Full documentation: [http://ribonets.github.io/PyDesign](http://ribonets.github.io/PyDesign)

##### Build the documentation #####

sphinx-build -b html doc ~/html

##### Installation #####

python setup.py install

##### Examples #####

* Classic optimization script 

> echo -e '(((((....)))))....\n....(((((....)))))' | classic_optimization.py -i -m random -e 500

