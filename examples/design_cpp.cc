/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 23.01.2014
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 * 
 * Compile with: g++ -fopenmp -g -Wall --std=c++11 -I/home/jango/local/include/RNAdesign -L/home/jango/local/lib/ -L/usr/lib64/ design_cpp.cc -lRNA -lRNAdesign -lm -o design_cpp
 */


// include standard library parts
#include <vector>
#include <string>
#include <iostream>

// include RNA header
extern "C" {
	#include "ViennaRNA/fold.h"
	#include "ViennaRNA/part_func.h"
}
// include RNAdesign
#include "RNAdesign.h"

// typedefs

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

float pf_fold(std::string& sequence, std::string& structure) {
    char* structure_cstr = new char[sequence.size()+1];
    float energy = pf_fold(sequence.c_str(), structure_cstr);
    structure = structure_cstr;
    delete structure_cstr;
    return energy;
}

//objective function: eos(1)+eos(2)+eos(3) - 3 * gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)
float objective_function(std::string& sequence, std::vector<std::string>& structures) {
    std::vector<float> eos;
    for (auto s : structures) {
        eos.push_back(energy_of_structure(sequence, s));
    }
    std::string pf_fold_struct;
    float gibbs = pf_fold(sequence, pf_fold_struct);
    
    float objective_difference_part = 0.0;
    for (unsigned int i=0; i < eos.size(); i++) {
        for (unsigned int j=i+1; j < eos.size(); j++) {
            objective_difference_part += pow(eos[i] - eos[j], 2);
        }
    }
    float eos_sum = 0;
    for (int n : eos)
        eos_sum += n;
    
    return eos_sum - eos.size() * gibbs + 1 * objective_difference_part;
}

// main program starts here
int main () {
	
	std::vector<std::string> structures;
	std::cout << "Input structures in dot-bracket (end with empty line): " << std::endl;
	while (true) {
	    std::string structure;
	    std::getline(std::cin, structure);
	    if (structure.empty())
	        break;
        else
	        structures.push_back(structure);
    }
    
    design::DependencyGraph<std::mt19937> dependency_graph(structures);
    std::cout << "Number of Sequences: " << dependency_graph.number_of_sequences() << std::endl;
    //std::cout << "Graph: " << dependency_graph.get_graphml() << std::endl;
    
    
    for (unsigned int n=0; n<10; n++) {
        std::string result_sequence = dependency_graph.get_sequence();
        float score = objective_function(result_sequence, structures);
    
        for (unsigned int i=0; i<10000; i++) {
            dependency_graph.sample_global();
            std::string current_sequence = dependency_graph.get_sequence();
            float this_score = objective_function(current_sequence, structures);
            
            if (this_score < score) {
                score = this_score;
                result_sequence = current_sequence;
            }
        }
        std::cout << result_sequence << "\t" << score << std::endl;
    }
    exit(0);
}
