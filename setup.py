import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "PyDesign",
    version = "1.0",
    author = "Stefan Hammer",
    author_email = "s.hammer@univie.ac.at",
    description = ("A wrapper framework to design RNA molecules."),
    license = "GPLv3",
    keywords = "RNA design synthetic biology",
    url = "http://github.com/ribonets/PyDesign",
    packages=['PyDesign'],
    test_suite="tests",
    install_requires=[
                    ],
    scripts=['bin/cofold_design.py',  'bin/generate_graphml.py', 'bin/multistate-design.py', 'bin/print_graphml.py'],
    long_description=read('README.rst'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
