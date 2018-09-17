#!/usr/bin/env python
from __future__ import print_function

'''
    Structure.py: Class for RNA secondary structures
    This implements various useful functions to manipulate and analyze secondary structures.
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2018"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

from collections import Counter

class RNAStructure(object):
    '''
    RNAStructure object to handle RNA secondary structures. Holds a pairtable internally to operate on.

    :param dotbracket: Dot-bracket notation of the RNA secondary structure as string
    '''

    def __init__(self, dotbracket):
        self._dotbracket = dotbracket
        self._pairtable = self._parseRNAStructure(dotbracket)

    def __str__(self):
        return self._dotbracket

    def __unicode__(self):
        return unicode(self._dotbracket)

    def _parseRNAStructure(self, structure):
        '''
        Parse RNA structure including pseudoknots

        :param structure: string, RNA secondary structure in dot-bracket notation
        :return: list of positions specifying the pair-table
        '''

        opening = "([{<"
        closing = ")]}>"

        stack = { op:list() for op in opening }
        bps = [-1]*len(structure)

        for i,c in enumerate(structure):
            for (op,cl) in zip(opening,closing):
                if c==op:
                    stack[op].append(i)
                elif c==cl:
                    j = stack[op].pop()
                    bps[i] = j
                    bps[j] = i

        return bps

    def resolvePKs(self):
        '''
        Removes any pseudoknots from the structure. Pseudoknotted base-pairs are missing in the returned structure!

        :return: list of pseudoknot-free structure strings in dot-bracket notation
        '''
        opening = "([{<"
        closing = ")]}>"
        structures = []
        c = Counter(self._dotbracket)
        for (op,cl) in zip(opening,closing):
            if c[op] > 0:
                newstruct = []
                for n in self._dotbracket:
                    if n == op:
                        newstruct.append('(')
                    elif n == cl:
                        newstruct.append(')')
                    else:
                        newstruct.append('.')
                structures.append(''.join(newstruct))
        return structures

    def removeLonelyPairs(self):
        '''
        Returns a dotbracket string without lonely pairs.

        :return: structure string in dot-bracket notation without any lonely pairs
        '''
        newstruct = ['.'] * len(self._dotbracket)
        for (i,j) in enumerate(self._pairtable):
            # this pair is stacking
            if i+1<j-1 and self._pairtable[i+1]==j-1:
                for n in [i, j, i+1, j-1]:
                    newstruct[n] = self._dotbracket[n]

        return ''.join(newstruct)

    def hasLonelyPairs(self):
        '''
        Returns a boolean whether the structure has lonely pairs
        
        :return: boolean
        '''
        return self.removeLonelyPairs() is self._dotbracket
