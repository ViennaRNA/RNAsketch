import RNAdesign as rd

# define structures
structures = ['(((((....)))))', '(((....)))....']
# construct dependency graph with these structures
dg = rd.DependencyGraphMT(structures)

# print this sequence
print dg.get_sequence()

# mutate globally for 1000 times and print
for i in range(0, 1000):
    dg.sample_global()
    print dg.get_sequence()

# print the amount of solutions
print 'Maximal number of solutions: ' + str(dg.number_of_sequences())
# print the amount of connected components
print 'Number of Connected Components: ' + str(dg.number_of_connected_components())
