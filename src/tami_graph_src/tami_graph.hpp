#include <iostream>
#include <fstream>
#include <string> 
#include <sstream> 
//#include<Eigen/Dense>
//#include<Eigen/Core>

#include <complex>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
//#include<experimental/filesystem>

#include "../tami_base_src/tami_base.hpp"

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


class TamiGraph{

    private:

    public:

        // structs and typedefs

        enum spin_type{up,dn};
        enum bubble_type{directed, both_in, both_out};
        enum label_type{labelled, unlabelled};
        enum vertex_type{three_leg, four_leg };

        ////////////////////////
        // Vertex info structure
        ////////////////////////

        //typedef boost::property< vertex_index_t, int> vertex_index;

        struct vertex_info {
            vertex_info(vertex_type type){
                type_=type;
                visited=0;

            if (type_==four_leg){  throw std::runtime_error("vertex type is 4 leg but no U_value given.");}

                }

            vertex_info(vertex_type type, double U_value){
                type_=type;
                U_value_=U_value;
                visited=0;
                }
            //  vertex_info(){ throw std::runtime_error("type has to be default constructible but we should never use the default constructor.");}
            vertex_info(){
            type_=three_leg; // assume it is a 3-leg vertex if no type specified
            visited=0;
            }
            
            vertex_type type_;
            double U_value_;
            int index_;
            //  size_t vertex_index;
            //int component;
            int visited;
            boost::property< boost::vertex_index_t, int> vertex_index;
            
        };

        /////////////////////////

        ////////////////////////
        // Vertex info structure
        ////////////////////////

        // This is a mirror of the g_struct from ami
        struct edge_info {

            TamiBase::g_struct g_struct_;
            // g_struct_ includes 

            // epsilon_t eps_;
            // alpha_t alpha_;
            // stat_type stat_;

            // where these are defined as 
            // typedef std::vector<int> epsilon_t;  /// the symbolic epsilon - not numeric - hence integers
            // typedef std::vector<int> alpha_t;
            // typedef enum {Bose,Fermi} stat_type ;



            // DEBUGGING ONLY
            int edge_number_;
            //TamiBase::stat_type edge_stat=
            int fermi_loop_id=-1;
            bool bubble=false;
            bool tadpole=false;
            spin_type spin=up;
            int band=0;
            std::pair<int, int> band_element;
            label_type label;

            std::vector<int> fourindex;

            edge_info(){ 
            label=unlabelled;
            }

            // todo fix this - need to remove loop id
            edge_info(int bnd, TamiBase::stat_type type ){
            band=bnd;
            g_struct_.stat_=type;	
            label=unlabelled;
            }

            edge_info(TamiBase::stat_type type, std::pair<int,int> be){
            g_struct_.stat_=type;
            band_element=be; // this is basically the original source-target labelling	
            label=unlabelled;
            }

            edge_info(TamiBase::stat_type type, std::vector<int> four){

            g_struct_.stat_=type;
            fourindex=four;
            label=unlabelled;	
            }

            edge_info(TamiBase::stat_type type, int loop_id, spin_type spin_val){ 
            g_struct_.stat_=type;
            fermi_loop_id=loop_id;
            spin=spin_val;
            label=unlabelled;

            // species info to g_struct
            // g_struct_.species_=(species_t)spin_val;
            // in theory R0 now knows the list of species 

            }

            edge_info(TamiBase::epsilon_t epsilon, TamiBase::alpha_t alpha, TamiBase::stat_type type, int loop_id, spin_type spin_val){
            g_struct_.eps_=epsilon;
            g_struct_.alpha_=alpha;
            g_struct_.stat_=type;
            fermi_loop_id=loop_id;
            spin=spin_val;
            label=unlabelled;
            }

            edge_info(TamiBase::stat_type type, int loop_id){ 
            g_struct_.stat_=type;
            fermi_loop_id=loop_id;
            label=unlabelled;
            }


            edge_info(TamiBase::stat_type type){ 
            g_struct_.stat_=type;
            label=unlabelled;
            }
            edge_info(TamiBase::epsilon_t epsilon, TamiBase::alpha_t alpha, TamiBase::stat_type type){
            g_struct_.eps_=epsilon;
            g_struct_.alpha_=alpha;
            g_struct_.stat_=type;
            label=unlabelled;
            }

            edge_info(TamiBase::epsilon_t epsilon, TamiBase::alpha_t alpha){
            g_struct_.eps_=epsilon;
            g_struct_.alpha_=alpha;
            g_struct_.stat_=static_cast<TamiBase::stat_type>(1); // 1 here means fermi;
            label=unlabelled;
            }
                


        }; // end edge_struct bracket


        ////////////////////////
        // graph_info structure
        ////////////////////////
        struct graph_info {
            
            int fermi_loop_id;
            int num_loops;
            int n_indep;
            int n_labelled;
            int ext_counts;
            bool is_bose;
            int ct_count=0;
                int sigma_ct_count=0;

                
        };

        //   A ->--- a =====b--->--B .  A and B are external vertices. a and b are the bubble vertices 

        struct bubble{

            std::vector<std::pair<int,int>> bubble_edges;
            std::vector<std::pair<int,int>> leg_edges;	

            int A,a,b,B;
  
         };


        typedef boost::adjacency_list <boost::vecS, boost::listS, boost::bidirectionalS, vertex_info, edge_info, graph_info> graph_t;
        typedef std::vector<graph_t> labels_t;
        typedef boost::graph_traits < graph_t >::vertex_descriptor vertex_t;
        typedef boost::graph_traits < graph_t >::edge_descriptor edge_t;

        typedef std::vector< edge_t> edge_vector_t;
        typedef std::vector< vertex_t> vertex_vector_t;

        typedef std::vector<int> git_perm_t;
        typedef std::vector<git_perm_t> git_perm_set_t;
        typedef std::vector<git_perm_set_t> git_perm_list_t;

        // used inside a member function
        typedef
            boost::adjacency_list<
                boost::vecS            // edge list
            , boost::vecS            // vertex list
            , boost::undirectedS     // directedness
            , float                  // property associated with vertices
            >
        cc_graph_t;





        struct git_pair{
                
            git_pair(graph_t g1, graph_t g2, git_perm_set_t pst){ 

            g1_=g1;
            g2_=g2;
            pst_=pst;

            }

            git_pair(){}		
                
            graph_t g1_;
            graph_t g2_;
            git_perm_set_t pst_;
            // !!!!!!!!!!!!!!! not sure if needed - hopefully not! // NewAmiCalc::solution_set s1_, s2_;	
        };

        struct graph_group {
            
            std::vector<graph_t> graph_vec;
            std::vector<git_pair> gp_vec;
            std::vector<graph_t> ct_vec; // this is not used I think due to order issues. 
            
            std::vector<double> prefactor;
            std::vector<labels_t> labels;
            // !!!!!!!!!!!!!!! not sure if needed - hopefully not! // std::vector<NewAmiCalc::solution_set> ss_vec; // solution set vector 
            
            // !!!!!!!!!!!!!!! not sure if needed - hopefully not! // std::vector<AmiSpec::spec_solution_set> sss_vec;// spectral solution set vector 
            
            // git permutation set 
            // git_perm_list_t perm_vec;
            int order_shift=0;
                
        };

        typedef std::vector<graph_group> gg_vec_t;
        typedef std::vector< gg_vec_t > gg_matrix_t;






        // read ggmp functions
        void read_ggmp(std::string folder, gg_matrix_t &ggm, int max_ord);
        void read_ggmp(std::string folder, gg_matrix_t &ggm, int min_ord, int max_ord);
            void graph_read(std::string filename, graph_t &g);
            void pair_read(std::string filename, git_pair &p);


        // ggm label function
        void ggm_label(gg_matrix_t &ggm, int min);
            void repeated_labelling(graph_t &g, bool &result);
                void fix_epsilons(graph_t &g);
                    void find_internal_fermionic_edges(graph_t &g, edge_vector_t &vector);
                    bool edge_alphas_are_equal(edge_t &one, edge_t &two, graph_t &g);
                    bool edge_alphas_are_negative(edge_t &one, edge_t &two, graph_t &g);
                void number_vertices(graph_t &g);
                void label_half_random(graph_t &g);
                    void reset_g(graph_t &g);
                        void find_bose_fermi_edges(graph_t &g, edge_vector_t &bose, edge_vector_t &fermi);
                    void find_external_vertices(graph_t &g, vertex_vector_t &v, edge_vector_t &edges);
                    void label_extern_legs(edge_vector_t &extern_vect_list,graph_t &g);
                    void label_and_find_tadpoles_ami(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges);
                        void find_tadpoles(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges, edge_vector_t &edge_a, edge_vector_t &edge_b );
                            void find_bosonic_edge_on_three_pointvert(graph_t &g, vertex_t vert, edge_t &bose_edge);
                        void label_indep_edge(edge_t &e, graph_t &g);
                        void assign_cons_label(graph_t &g, vertex_t &vin, edge_t &fermi_one, edge_t &fermi_two, edge_t &bose);
                            void mark_labelled(edge_t &e, graph_t &g);
                            void add_labels(edge_t &one, edge_t &two, graph_t &g);
                            void subtract_labels(edge_t &one, edge_t &two, graph_t &g);
                            void label_indep_epsilon(edge_t &e, graph_t &g);
                    void find_internal_vertices(graph_t &g, vertex_vector_t &vector);
                    void sort_labelled_unlabelled_adjacent(graph_t &g, vertex_t &vin, vertex_vector_t &labelled_adj_vertices,vertex_vector_t &unlabelled_adj_vertices, edge_vector_t &labelled_edges, edge_vector_t &unlabelled_edges);
                    int random_int( int min, int max);  // this is a small function to get a random integer.  For real numbers use auto roll_dice example below.
                    bool label_consv_momentum(graph_t &g, vertex_t &vin, vertex_vector_t &labelled_adj_vertices,vertex_vector_t &unlabelled_adj_vertices, edge_vector_t &labelled_edges, edge_vector_t &unlabelled_edges);
                        void print_all_edge_info(graph_t &g);
                void check_momentum_conservation(graph_t &g, bool &result);
                void find_unlabelled_fermionic_edges(graph_t &g, edge_vector_t &vector);
            void sys_label(graph_t &g, bool &result);
                int graph_order(graph_t &g);
                    void find_internal_edges_stat(graph_t &g, edge_vector_t &vector, TamiBase::stat_type requested);
                void find_force_LR_vertices(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv);
                void delete_legs(graph_t &g, vertex_vector_t &v, edge_vector_t &edges);
                void find_unlabelled_edges(graph_t &g, edge_vector_t &vector);
                void fix_force_labels(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv);
                    void find_path_between_vertices(graph_t &g, vertex_t &v1, vertex_t &v2, bool &success, edge_vector_t &final_ev, vertex_vector_t &final_vv);
                        void get_adjacent_vertex_stat(graph_t &g, vertex_t &vin, vertex_t &next, edge_t &eout,  TamiBase::stat_type stat);
                void put_back_legs(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv);

        // print ggm function
        void print_ggm( gg_matrix_t &ggm);

        // ggm to r0
        void graph_to_R0(graph_t &g, TamiBase::g_prod_t &R0);

        // extract bose alphas for non-hubbard
        void extract_bose_alphas(graph_t g, std::vector<TamiBase::alpha_t> &bose);

        // Renorm PT CT diagrams
        void generate_sigma_ct( graph_t &g_in, std::vector< graph_t> &ct_vec, int maxdots);
            void find_fermionic_edges(graph_t &g, edge_vector_t &vector);
            void combinations_repetition(int n, int r, std::vector< std::vector<int>> &list );
                void combinations_repetition_util(std::vector<int> &chosen, std::vector< std::vector<int>> &list, int index, int r, int start, int end);// allows repeated values in combinations 
            void insert_chain(graph_t &g_in, graph_t &g_out,  edge_t &ei, int length);

        // prefactor from fermi loops and ct insertions
        double get_prefactor(graph_t &g, int order);
            int fermi_connected_components(graph_t &g);
                void find_external_vertices(graph_t &g, vertex_vector_t &v, edge_vector_t &edges);
    


                





        




};