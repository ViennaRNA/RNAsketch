# PyDesign Library for designing RNA molecules #


## Installation ##

```
python setup.py install
```

## Build the documentation ##

Documentation available online: [http://ribonets.github.io/PyDesign](http://ribonets.github.io/PyDesign)

Call these commands to build the html documentation:
```
cd doc
make html
```

## Provided Example Scripts ##

### Design a Multistate Riboswitch ###

```
echo -e '(((((....)))))....\n....(((((....)))))' | multistate-riboswitch.py -i -m random -e 500
```

### Design a sRNA mediated translational regulation ###
cofold_design.py  

### Design a simple Thermoswitch ###

multistate-thermoswitch.py

### Display the Dependency Graph ###

generate_graphml.py  print_graphml.py
