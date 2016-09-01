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

### Classic optimization script ###

```
echo -e '(((((....)))))....\n....(((((....)))))' | classic_optimization.py -i -m random -e 500
```

