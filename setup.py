import os
from os.path import join
from distutils.core import setup

def read(fname):
    return open(join(os.path.dirname(__file__), fname)).read()

setup(
    name = "RNAsketch",
    version = "1.4",
    author = "Stefan Hammer",
    author_email = "s.hammer@univie.ac.at",
    description = ("A wrapper framework to design RNA molecules."),
    license = "GPLv3",
    keywords = "RNA design synthetic biology",
    url = "http://github.com/ribonets/RNAsketch",
    packages=['RNAsketch'],
    test_suite="tests",
    install_requires=[
                    ],
    scripts=[ join('bin', 'design-cofold.py'), join('bin', 'design-generategraphml.py'), join('bin', 'design-multistate.py'), join('bin', 'design-printgraphml.py'), join('bin', 'design-thermoswitch.py'), join('bin', 'design-ligandswitch.py')],
    long_description=read('README.rst'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
