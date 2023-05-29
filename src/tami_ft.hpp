#pragma once

#include <filesystem>
#include <tami_base.hpp>

#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <unistd.h>
#include <random>

#include <chrono>
#include <thread>

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


class FermiTree {
public:



/*/* 
typedef std::vector<int> epsilon_t;
typedef std::vector<int> alpha_t;
struct pole_struct {
    pole_struct(epsilon_t eps, alpha_t alpha) {
      eps_ = eps;
      alpha_ = alpha;
    }

    pole_struct() {}

    epsilon_t eps_;
    alpha_t alpha_;
    // / Index that specifies which frequency it is a pole with respect to.
    int index_;
    // / The multiplicity of the pole, starts at 1 and increments as needed.
    int multiplicity_ = 1;
    int der_ = 0;              //< Counter for derivatives.
    std::vector<int> which_g_; //< Index to identify which `g_struct` a pole originated from.

    // / Experimental component of Spectral evaluation.
    alpha_t x_alpha_;
  };
 */
enum operation{add,mult,end};

// Vertex info  

struct vertex_info {
  
  vertex_info(){
visited=0;
depth_=-1;  // default value so easier to catch. 
operation_=-1; // set to nonsense value for catching bugs.
prefactor_=1;
}

vertex_info(int op, double depth){
  operation_=op;
  depth_=depth;
  prefactor_=1;
}

vertex_info(TamiBase::pole_struct pole,  int op){
  operation_=op;
  pole_=pole;
  prefactor_=1;
}


vertex_info(TamiBase::pole_struct pole, double prefactor,  int op){
  operation_=op;
  pole_=pole;
  prefactor_=prefactor;
}


vertex_info(TamiBase::pole_struct pole, std::complex<double> value, int op, double depth){
  operation_=op;
  depth_=depth;
  pole_=pole;
  value_=value;
  prefactor_=1;
}

vertex_info(TamiBase::pole_struct pole, std::complex<double> value, int op, double depth, double prefactor){
  operation_=op;
  depth_=depth;
  pole_=pole;
  value_=value;
  prefactor_=prefactor;
}

 
  std::complex<double> value_;
  int operation_;

  int visited;
  boost::property< boost::vertex_index_t, int> vertex_index;
  int index_;
  int depth_;
  
  TamiBase::pole_struct pole_;
  double prefactor_;
  
};

// Graph global info 

struct graph_info {
  
  int max_depth_=0;
  // double overall_prefactor=1.0;
     
};

/////////////////////////

////////////////////////
// Edge info structure
////////////////////////

struct edge_info {

int depth_; // say this is the depth of either the target or the source


  edge_info(){ 
depth_=-1;
}

// todo fix this - need to remove loop id
edge_info(int depth){
depth_=depth;
}


}; // end edge_struct bracket




typedef boost::adjacency_list < boost::vecS, boost::listS, boost::bidirectionalS,vertex_info, edge_info, graph_info > fermi_tree_t;
typedef boost::graph_traits < fermi_tree_t >::vertex_descriptor vertex_t;
typedef boost::graph_traits < fermi_tree_t >::edge_descriptor edge_t;

// helper functions
void number_vertices(fermi_tree_t &g);





// Initialize the tree to have one 
void initialize_ft( fermi_tree_t &ft);
void initialize_ft( fermi_tree_t &ft, operation op);
void initialize_ft( fermi_tree_t &ft, TamiBase::pole_struct &pole);
void initialize_ft( fermi_tree_t &ft, TamiBase::pole_struct &pole, double &prefactor);

void plist_to_ft(TamiBase::pole_array_t &plist, FermiTree::fermi_tree_t &ft);
void plist_to_ft(TamiBase::pole_array_t &plist,double sign, FermiTree::fermi_tree_t &ft);

// this changes every prefactor and is not in general the same as an overall prefactor. 
void update_prefactors(fermi_tree_t &ft, double sign);
vertex_t get_root(fermi_tree_t  &ft); 
bool is_empty_ft(fermi_tree_t &ft);

void mult_prefactor(fermi_tree_t &ft, double sign);



void get_roots(fermi_tree_t  &ft, std::vector<vertex_t> &vv);


// TODO: Do I want it to destroy the level info on call?
void get_next_level( fermi_tree_t &ft1, vertex_t &root1,  std::vector<vertex_t > &next_level);



void print_vertex(vertex_t &v, fermi_tree_t &ft);


void print_graph(fermi_tree_t &ft);

std::string pretty_print_ft(fermi_tree_t &ft1, vertex_t &v);
std::string pretty_print_pole(TamiBase::pole_struct &pole);
std::string pretty_print(fermi_tree_t &ft, vertex_t &v);
std::string pretty_print(fermi_tree_t &ft);
std::string pretty_print(fermi_tree_t &ft, std::vector<vertex_t> &vv, int op);
std::string pretty_print_level(fermi_tree_t &ft, std::vector<vertex_t> &vv, int source_op);

// this should copy root1 from ft1 to ft2 below root2
// will fill the pair with the vertex_t of the new vertex on ft1 and ft2
void copy_vertex(vertex_t &root1, vertex_t &root2, fermi_tree_t &ft1, fermi_tree_t &ft2, std::pair<vertex_t, vertex_t> &map);


// copy the level below root1 to below root2 
// will populate a map relating the next level of root1 to the new vertices created on ft2 below root2
void copy_level(vertex_t &root1, vertex_t &root2 ,fermi_tree_t &ft1, fermi_tree_t &ft2,std::vector<std::pair<vertex_t, vertex_t>> &map_vec);


//copy entire tree from ft1 starting at root1 below root2 
// either include root1 or don't.  which is it?
void copy_tree(fermi_tree_t &ft1, fermi_tree_t &ft2,vertex_t &root1, vertex_t &root2);




// could either place 
fermi_tree_t add_ft(fermi_tree_t ft1, fermi_tree_t ft2);



fermi_tree_t mult_ft(fermi_tree_t ft1, fermi_tree_t ft2);



std::complex<double> eval_ft(fermi_tree_t &ft1, vertex_t &v);



FermiTree();

private:
};


