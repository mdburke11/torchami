#include <iostream>
#include <fstream>
#include <string> 
#include <sstream> 
#include<Eigen/Dense>
#include<Eigen/Core>

#include <complex>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include<experimental/filesystem>

#include "tami_base.hpp"

// Boost headers 
#include <boost/config.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>

// #include <boost/program_options.hpp>
#include "boost/graph/graphviz.hpp"
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
//#include <boost/graph/depth_first_search.hpp>
#include <boost/range/irange.hpp>
//#include <boost/pending/indirect_cmp.hpp>
//#include <boost/graph/undirected_dfs.hpp>
//#include <boost/cstdlib.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

int add(int, int);
