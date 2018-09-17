#!/usr/bin/env python
from __future__ import print_function

'''
    HotKnots.py: Wrapper for the HotKnots software.
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2018"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

import math
import subprocess as sub
import os
#import tempfile
import re
#import shutil

def _get_exec_path(exec_name):
  """ If the HOTKNOTS environment variable is set, use that as the directory
  of the hotknots executables. Otherwise, have Python search the PATH directly. """
  if 'HOTKNOTS' in os.environ:
    return os.environ['HOTKNOTS'] + '/bin/' + exec_name
  else:
    return exec_name

def _setup_args(**kargs):
    '''
    Takes the arguments from the various calls and generates a string array for the command-line call.

    :param kargs: List of the function arguments
    :return: List of strings with the command line arguments
    :return: String with command line input
    '''

    # Check for structure input
    evalEnergy = False
    try:
        if kargs['structure']:
            evalEnergy = True
    except:
        pass

    binary = _get_exec_path('HotKnots')
    if evalEnergy:
        binary = _get_exec_path('computeEnergy')

    args =  [binary, '-m', kargs['model'], '-s', kargs['sequence']]

    #'-c' # sequence file includes structural constraints
    if evalEnergy:
        args += [kargs['structure']]
    else:
        args += ['-noPS']

    if not kargs['allowGU']:
        args += ['-noGU']

    # command line arguments must be strings
    args = list(map(str, args))

    return args

def _call_with_pipe(args, cmd_input, outputdir=False, raiseOnError = True):
    '''
    Performs the command line call.

    :param args: List of strings containing the program arguments
    :param cmd_input: String containing the program stdin input
    :return: String with output of the program
    '''

    # create the output directory
    if outputdir and not os.path.exists(_get_exec_path('output')):
        os.makedirs(_get_exec_path('output'))

    #print(" ".join(['# '] + args))
    p=sub.Popen(args, stdin = sub.PIPE, stdout = sub.PIPE, stderr = sub.PIPE, cwd = _get_exec_path(''))
    output, error = p.communicate(cmd_input)

    e = error.decode().rstrip('\n')
    if e and raiseOnError:
        raise IOError('HotKnots returned an error:\n' + e)

    # delete output directory
    #if outputdir:
    #    shutil.rmtree(_get_exec_path('output'), ignore_errors=True)

    return output

def _parse_output(output):
    '''
    Regex that parses the output of the programs.

    :param output: Output string of the pKiss program
    :param multistate: Boolean stating if it is a multistate output
    :param window: Boolean stating if it is window mode output
    :return: Either a tuple of structure/shape and energy values, or in window mode a dictionary with the start position as key and the tuple as value.
    '''

    eval_pattern = re.compile(r'^(?P<model>[A-Za-z&]+)\s+(?P<energy>-?\d+\.?\d*(e-?\d+)?)\s+(?P<nodangling>-?\d+\.?\d*(e-?\d+)?)', flags = re.M)

    mfe_pattern = re.compile(r'S(?P<state>\d+):\s+(?P<structure>[\{\}\<\>\.\(\)\[\]AaBb]+)\s+(?P<energy>-?\d+\.?\d*(e-?\d+)?)', flags = re.M)

    m = re.search(eval_pattern, output)
    if m:
        return float(m.group('energy'))
    else:
        values = {}
        for m in re.finditer(mfe_pattern, output):
            values[int(m.group('state'))] = (m.group('structure'), float(m.group('energy')))
    if not values:
        raise IOError("Could not parse HotKnots output:\n" + output)
    return values

def mfe(sequence, allowGU = True, model = 'DP'):
    '''
    MFE prediction using HotKnots with the given parameters

    :param sequence: String containing the RNA sequence [AUGC]
    :param allowGU: Boolean whether to allow GU base-pairs (default: True)
    :param model: String specifying the PK model; RE for Rivas&Eddy, DP for Dirks&Pierce (default), CC for Cao&Chen
    :return: Tuple of string and float containing the mfe structure and energy
    '''

    args = \
    _setup_args(sequence = sequence, allowGU = allowGU, model = model)

    output = _call_with_pipe(args, None, outputdir = True)
    try:
        d = _parse_output(output)
        returnvalue = min(d.values(), key=lambda x: x[1])
    except:
        raise IOError("Could not parse HotKnots output:\n" + output)
    return returnvalue

def subopt(sequence, allowGU = True, model = 'DP'):
    '''
    Subpotimal structure prediction using HotKnots with the given parameters

    :param sequence: String containing the RNA sequence [AUGC]
    :param allowGU: Boolean whether to allow GU base-pairs (default: True)
    :param model: String specifying the PK model; RE for Rivas&Eddy, DP for Dirks&Pierce (default), CC for Cao&Chen
    :return: Tuple of string and float containing the mfe structure and energy
    '''

    args = \
    _setup_args(sequence = sequence, allowGU = allowGU, model = model)

    output = _call_with_pipe(args, None, outputdir = True)

    return _parse_output(output)

def eval(sequence, structure, allowGU = True, model = 'DP'):
    '''
    Energy evaluation of a sequence structure pair using HotKnots with the given parameters

    :param sequence: String containing the RNA sequence [AUGC]
    :param structure: Dot-bracket string specifying the secondary structure [.(){}[]<>]
    :param allowGU: Boolean whether to allow GU base-pairs (default: True)
    :param model: String specifying the PK model; RE for Rivas&Eddy, DP for Dirks&Pierce (default), CC for Cao&Chen
    :return: Tuple of string and float containing the mfe structure and energy
    '''

    args = \
    _setup_args(sequence = sequence, structure = structure, allowGU = allowGU, model = model)

    output = _call_with_pipe(args, None)

    return structure, _parse_output(output)

# Test for binary precense
try:
    output = _call_with_pipe([_get_exec_path('HotKnots'), '-h'], None, raiseOnError = False)
    output = _call_with_pipe([_get_exec_path('computeEnergy'), '-h'], None, raiseOnError = False)
except:
    raise ImportError('HotKnots not found. Please install HotKnots and speficy the HotKnots folder path in an environmental variable called HOTKNOTS!')
