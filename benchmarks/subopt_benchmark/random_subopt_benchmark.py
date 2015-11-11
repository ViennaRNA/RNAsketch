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
    parser.add_argument("-n", "--number", type=int, default=40, help='Number of subpots to generate')
    parser.add_argument("-l", "--length", type=int, default=60, help='Length of random sequence')
    parser.add_argument("-s", "--structures", type=int, default=4, help='Number of structures as constraints input')
    args = parser.parse_args()
    
    dg1 = rd.DependencyGraphMT(["." * args.length])
    
    # benchmark results
    graph_construction_count = 0
    input_sequence_found_count = 0
    
    for n in range(0, args.number):
        dg1.set_sequence()
        input_sequence = dg1.get_sequence()
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
        print(*structures, sep="\n")
        
        # try to construct dependency graph, catch errors and timeouts
        dg = None
        
        try:
            with timeout(seconds=8):
                dg = rd.DependencyGraphMT(structures)
        except Exception as e:
            print( "Error: %s" % e )
        
        if (dg is not None):
            graph_construction_count += 1
            print ("NOS: %d" % dg.number_of_sequences())
            
            # do optimization to mfe
            dg.set_sequence()
            for o in range(0, 10000):
                current_sequence = dg.get_sequence()
                if (current_sequence == input_sequence):
                    input_sequence_found_count += 1
                    print ("Input sequence found!")
                    break
                else:
                    for pos, c in enumerate(input_sequence):
                        if c != current_sequence[pos]:
                            dg.mutate(pos)
                            break
            
            print (dg.get_sequence())
        
        print ("-" * args.length)
        
    # Print final benchmark result
    print ("%d%% of dependency graphs could be built" % (float(graph_construction_count) / float(args.number) * 100))
    print ("%d%% of input sequences could be found" % (float(input_sequence_found_count) / float(args.number) * 100))
    
if __name__ == "__main__":
    main()

