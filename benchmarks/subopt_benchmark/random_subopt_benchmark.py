from __future__ import print_function
import RNAdesign as rd
import RNA
import argparse
import sys
import random
import signal
import time


class timeout:
    def __init__(self, seconds=10, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise Exception(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

def main():
    parser = argparse.ArgumentParser(description='Generate a random sequence, get suboptimal structures and see if dependencygraph exists.')
    parser.add_argument("-n", "--number", type=int, default=100, help='Number of tests to generate')
    parser.add_argument("-l", "--length", type=int, default=60, help='Length of random sequence')
    parser.add_argument("-s", "--structures", type=int, default=4, help='Number of structures as constraints input')
    parser.add_argument("-m", "--mode", type=str, default='sample', help='Mode for getting a new sequence: sample, sample_local, sample_global')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()
    
    print ("# random_subopt_benchmark.py")
    print ("# Options: number={0:d}, length={1:d}, structures={2:d}, mode={3:}".format(args.number, args.length, args.structures, args.mode))
    print ("# length", "structures", "graph_construction", "num_cc", "max_special_ratio", "mean_special_ratio", "nos", "construction_time", "sample_time", sep=";")
    
    rd.initialize_library(args.debug)
    
    dg1 = rd.DependencyGraphMT(["." * args.length])
    
    # benchmark results
    graph_construction_count = 0
    input_sequence_found_count = 0
    
    for n in range(0, args.number):
        dg1.sample()
        input_sequence = dg1.get_sequence()
        if(args.debug):
            print (input_sequence)
        
        # calculate a RNAsubopt
        subopt_solution = RNA.subopt(dg1.get_sequence(), "N" * args.length, 600, None)
        #for i in range(0, subopt_solution.size()-1):
        #    print (subopt_solution.get(i).structure, subopt_solution.get(i).energy)
        
        # get random structures
        structures = []
        for i in range(0, args.structures):
            rand = random.randint(0, subopt_solution.size()-1)
            structures.append(subopt_solution.get(rand).structure)
        if (args.debug):
            print(*structures, sep="\n")
        
        # try to construct dependency graph, catch errors and timeouts
        dg = None
        graph_construction = 0
        construction_time = 0.0
        sample_time = 0.0
        max_special_ratio = 0
        mean_special_ratio = 0
        num_cc = 0
        nos = 0
        
        try:
            with timeout(seconds=15):
                start = time.clock()
                dg = rd.DependencyGraphMT(structures)
                construction_time = time.clock() - start
        except Exception as e:
            print( "Error: %s" % e )
        
        if (dg is not None):
            graph_construction = 1
            num_cc = dg.number_of_connected_components()
            nos = dg.number_of_sequences()
            
            special_ratios = []
            for cc in range(0, num_cc):
                special_ratios.append(float(len(dg.special_vertices(cc)))/float(len(dg.component_vertices(cc))))
            
            max_special_ratio = max(special_ratios)
            mean_special_ratio = sum(special_ratios)/len(special_ratios)
            
            # do sampling of sequences
            start = time.clock()
            dg.sample()
            for o in range(0, 10000):
                current_sequence = dg.get_sequence()
                # sample sequence completely random
                if args.mode == 'sample':
                    mut_nos = dg.sample()
                elif args.mode == 'sample_global':
                    mut_nos = dg.sample_global()
                elif args.mode == 'sample_local':
                    mut_nos = dg.sample_local()
                else:
                    sys.stdout.write("Wrong sample argument: " + args.mode + "\n")
                    sys.exit(1)
            sample_time = time.clock() - start

        
        print (args.length, args.structures, graph_construction, num_cc, max_special_ratio, mean_special_ratio, nos, construction_time, sample_time, sep=";")
    
if __name__ == "__main__":
    main()
