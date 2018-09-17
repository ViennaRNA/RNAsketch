#!/usr/bin/env python
from __future__ import print_function

try:
    from RNAsketch import *
except ImportError, e:
    print(e.message)
    exit(1)

import RNAblueprint as rbp
import argparse
import sys
import os
import time
from collections import Counter

def main():
    parser = argparse.ArgumentParser(description='Design a multi-stable riboswitch similar using Boltzmann sampling.')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-q", "--package", type=str, default='vrna', help='Chose the calculation package: hotknots, pkiss, nupack, or vrna/ViennaRNA (default: vrna)')
    parser.add_argument("-j", "--objective", type=str, default='1', help='Chose the objective function: 1 for abs differences and 2 for squared (default: 1)')
    parser.add_argument("-T", "--temperature", type=float, default=37.0, help='Temperature of the energy calculations.')
    parser.add_argument("-n", "--number", type=int, default=1000, help='Number of designs to generate')
    parser.add_argument("-m", "--model", type=str, default='stacking', help='Model for getting a new sequence: uniform, nussinov, basepairs, stacking')
    parser.add_argument("-s", "--stop", type=int, default=0, help='Stop optimization run of unpaired bases if no better solution is aquired after (stop) trials. 0 is no unpaired bases optimization.')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-c", "--csv", default=False, action='store_true', help='Write output as semi-colon csv file to stdout')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()
    
    if (args.debug):
        print("# Options: number={0:d}, model={1:}, stop={2:d}, package={3:}, temperature={4:}".format(args.number, args.model, args.stop, args.package, args.temperature))
    rbp.initialize_library(args.debug, args.kill)
    # define structures
    structures = []
    constraint = ''
    start_sequence = ''

    if (args.input):
        data = ''
        for line in sys.stdin:
            data = data + '\n' + line
        (structures, constraint, start_sequence) = read_input(data)
    elif (args.file is not None):
        if (args.debug):
            print("# Input File: {0:}".format(args.file))
        (structures, constraint, start_sequence) = read_inp_file(args.file)
    else:
        structures = ['((((....))))....((((....))))........',
            '........((((....((((....))))....))))',
            '((((((((....))))((((....))))....))))']
        constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    # try to construct dependency graph, catch errors and timeouts
    construction_time = 0.0
    sample_time = 0.0

    # remove lonely pairs
    structures = [RNAStructure(s).removeLonelyPairs() for s in structures]

    # general DG values
    if (args.debug):
        print("# " + "\n# ".join(structures) + "\n# " + constraint)

    design = get_Design(structures, start_sequence, args.package, args.temperature)

    # print header for csv file
    if (args.csv):
        print(";".join(["stop",
                    "model",
                    "score",
                    "num_mutations",
                    "construction_time",
                    "sample_time",
                    design.write_csv_header()]))

    target_energies, offsets, construction_time = getTargetEnergy(structures, args)
    if (args.debug):
        print("# Target Energies are: ", target_energies)
    bs = BalancedSamples(structures, target_energies, offsets, energy_step=0.5, args=args)

    count = 0
    for b in sorted(bs.keys()):
        if count > args.number:
            break
        count += 1

        design = get_Design(structures, bs[b]['seq'], args.package, args.temperature)
        #out = '$;' + str(b) + ';'
        #for i in range(0, design.number_of_structures):
        #    out = ';'.join([out, str(bs[b]['energies'][i]), str(bs[b]['energies'][i]+offsets[i]), str(design.eos[str(i)])])
        #print(out)
        #print('$ simple model: ', bs[b]['energies'], ' viennaRNA: ', design.eos)
        objective = calculate_objective
        if (args.objective == 2):
            objective = squared_objective
        # if args.stop is not 0, we want to optimize unpaired positions in the structure
        if (args.stop):
            score, num, sample_time = local_optimization(design, objective, args)
        else:
            score = objective(design)
            num = 0
        # output sequence
        if (args.csv):
            print(args.stop,
                    "\"" + args.model + "\"",
                    score,
                    num, # number of sequences sampled
                    construction_time,
                    sample_time, # sample time until now
                    design.write_csv(), sep=";")
        else:
            print(design.write_out(score))

def getTargetEnergy(structures, args):
    sampler = RPSampler(structures, model=args.model, weights=[1]*len(structures), temperature=args.temperature, stacksize=1000, StopConstruct=True, debug=args.debug)

    # get new sequecne
    newsample, energies = sampler.dump_new_stack()
    # get energy offsets
    offsets = getEnergyOffsets(structures, newsample, energies, args)
    # find target energy_step
    phi = 99999.9
    target_energies= []
    for i, s in enumerate(newsample):
        current_phi = getPhi(energies[i], offsets)
        if current_phi < phi:
            phi = current_phi
            target_energies = [np.mean(energies[i].values() + offsets.values())]*len(structures)
            if (args.debug):
                print('# Curren Phi and Simple Target Energies: ', current_phi, energies[i])
    # correct target energies with offsets
    for i, o in offsets.iteritems():
        target_energies[i] -= o

    return target_energies, offsets, sampler.construction_time

def getEnergyOffsets(structures, newsample, energies, args):
    nstr = len(structures)
    offsets = np.zeros( (len(newsample), nstr) )
    # iterate over sample
    for i, s in enumerate(newsample):
        # get design object
        design = get_Design(structures, s, args.package, args.temperature)
        # iterate over structures
        for t in range(0, nstr):
            # calculate offset between turner eos and simple model eos
            offsets[i,t] = design.eos[str(t)] - energies[i][t]
    # calculate mean offsets
    mean_offsets = {}
    for t in range(0, nstr):
        mean_offsets[t] = np.mean(offsets[:,t])
    if (args.debug):
        print('# mean offsets are: ', mean_offsets)
    return mean_offsets

def getPhi(energies, offsets):
    mean_eos = np.mean(energies.values() + offsets.values())
    phi = 0
    for i, eos in energies.iteritems():
        phi += abs(eos + offsets[i] - mean_eos)
    return phi

def BalancedSamples(structures, target_energies, offsets, args, energy_step=0.5):
    BalancedSample = {}
    # construct redprint sampler object
    nstr = len(structures)
    wastefactor = 20
    number = 1000
    sampler = RPSampler(structures, model=args.model, weights=([1.0] * nstr), gcweight=1.0, temperature=args.temperature, stacksize=(wastefactor*number), debug=args.debug)

    for shift in np.arange(0, 9999, (energy_step)):
        te = [x-shift for x in target_energies]
        if (args.debug):
            print("# Current Target energies are: ", te)
        AdmissibleSample = Sample(sampler, nstr, te, target_GC=0.5, number=number, args=args)
        for s in AdmissibleSample:
            BalancedSample[getPhi(s['energies'], offsets)] = s
        # Stop criterion
        eos = []
        for s in AdmissibleSample:
            eos_mean = np.mean(s['energies'].values() + offsets.values())
            eos.append(eos_mean)
        tartet_energy = target_energies[0]+offsets[0]
        if (args.debug):
            print('# Stop: ', abs(np.mean(eos)-tartet_energy))
            print("# Already found: ", len(BalancedSample)/float(args.number), "%")
        if (abs(abs(np.mean(eos)-tartet_energy)) > energy_step) and len(BalancedSample) > args.number:
            break
    return BalancedSample

def Sample(sampler, nstr, target_energies, target_GC, args, target_energy_eps = 0.20, target_GC_eps=0.05, maxiterations=10, number=1000):
    # weights = [math.exp(1/((args.temperature + 273.15)*0.00198717))] * nstr

    AdmissibleSample = []
    # count iterations
    count = 0
    while count < maxiterations:
        count += 1
        # get new sequences
        if (args.debug):
            print("# Weights: ", sampler.weights)
            print('# GC weight: ', sampler.gcweight)
        newsample, energies = sampler.dump_new_stack()

        # get average structue energies for newsample
        eos = np.zeros( (len(newsample), nstr) )
        GC_freq = []

        for i, s in enumerate(newsample):
            # count GC content
            c = Counter(s)
            sigma = 2.0 # one per nucleotide, laplace
            GC = (c['G'] + c['C'] + sigma) / (len(s) + 2*sigma)
            GC_freq.append(GC)
            # add if it is eps-admissible
            admissible = True
            if not (1-target_GC_eps <= GC/target_GC <= 1+target_GC_eps):
                admissible = False
            for t in range(0, nstr):
                # add to eos np array in any case
                eos[i,t] = energies[i][t]
                # check if eps admissible
                #print(eos[i,t], eos[i,t]/target_energies[t], target_energies[t])
                if not (1-target_energy_eps <= eos[i,t]/target_energies[t] <= 1+target_energy_eps):
                    admissible = False
            if admissible:
                #print('# is admissible:', eos[i,:], GC)
                AdmissibleSample.append({'seq': s, 'energies': energies[i]})
        # update weights
        for t in range(0, nstr):
            e_mean = np.mean(eos[:,t])
            if (args.debug):
                print('# Energy mean: ', str(t), e_mean)
            # exp version
            sampler.weights[t] = sampler.weights[t] * (1.1**(e_mean-target_energies[t]))
            # Yann old version without positive e_mean
            #weights[t] = weights[t]*target_energies[t]/e_mean
        # update gcweight
        GC_mean = np.mean(GC_freq)
        if (args.debug):
            print('# GC mean: ', GC_mean)
            print('# Found for current Target: ', len(AdmissibleSample)/float(number), '%')
        sampler.gcweight = sampler.gcweight * target_GC/GC_mean
        # return if large enough
        if len(AdmissibleSample) >= number:
            break
    return AdmissibleSample

def local_optimization(design, objective, args):
    # in case we want to optimize unpaired positions we can call this here

    dg = rbp.DependencyGraphMT(design.structures)
    start = time.clock()

    (score, number_of_mutations) = adaptive_walk_fixed(dg, design, objective_function=objective, number=args.stop, mode='sample_clocal', progress=args.progress)

    sample_time = time.clock() - start

    return score, number_of_mutations, sample_time

def squared_objective(design, weight=0.5):
    return calculate_objective_1(design) + weight * calculate_objective_2_squared(design)

if __name__ == "__main__":
    main()
