use RNAdesign;

# define structures
@structures = ['(((((....)))))', '(((....)))....'];
# construct dependency graph with these structures
$dg = new RNAdesign::DependencyGraphMT(@structures);

# print this sequence
print $dg->get_sequence()."\n";

# mutate globally for 1000 times and print
for($i=0; $i<1000; $i++) {
    $dg->sample_global();
    print $dg->get_sequence()."\n";
}
# print the amount of solutions
print 'Maximal number of solutions: '.$dg->number_of_sequences()."\n";
# print the amount of connected components
print 'Number of Connected Components: '.$dg->number_of_connected_components()."\n";
