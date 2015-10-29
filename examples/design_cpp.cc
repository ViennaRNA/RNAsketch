/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 23.01.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "energy.h"


// include standard library parts
#include <sstream>

// include boost components

// typedefs

// functions
float energy_of_structure (std::string& sequence, std::string& structure);

float fold (std::string& sequence, std::string& structure);


extern "C" {
#include "ViennaRNA/fold.h"
}

// typedefs


// global variables


// functions

float energy_of_structure(std::string& sequence, std::string& structure) {
    float energy = energy_of_structure(sequence.c_str(), structure.c_str(), 0);
    return energy;
}

float fold(std::string& sequence, std::string& structure) {
    char* structure_cstr = new char[sequence.length()+1];
    float energy = fold(sequence.c_str(), structure_cstr);
    structure = structure_cstr;
    delete structure_cstr;
    return energy;
}
