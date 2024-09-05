#include "tami_graph.hpp"

TamiGraph::TamiGraph() : rand_gen(0), rand_dist(0, 1.0) {
  std::cout << "TamiGraph Constructor Called" << std::endl;

  initialize(TamiBase::Sigma);
  // auto roll_dice=std::bind ( rand_dist, rand_gen);
}

TamiGraph::TamiGraph(TamiBase::graph_type type, int seed)
    : rand_gen(seed), rand_dist(0, 1.0) {

  std::cout << "TamiGraph Constructor Called" << std::endl;

  initialize(type);
  // auto roll_dice=std::bind ( rand_dist, rand_gen);
  // std::cout<<"Testing random number generator rand="<< roll_dice()
  // <<std::endl;
}

TamiGraph::TamiGraph(TamiBase::graph_type type, int dim, int seed)
    : rand_gen(seed), rand_dist(0, 1.0) {

  std::cout << "TamiGraph Constructor Called" << std::endl;

  initialize(type, dim);
  // auto roll_dice=std::bind ( rand_dist, rand_gen);
  // std::cout<<"Testing random number generator rand="<< roll_dice()
  // <<std::endl;
}

void TamiGraph::initialize(TamiBase::graph_type type, int dim) {

  graph_type = type;
  // TamiBase::graph_type type=AmiCalc::Sigma;

  // std::cout<< type <<std::endl;
  ami_parameters.TYPE_ = type; // type;

  // create_starter_graph(current_graph);

  // TODO: is this needed? -MB current_state.initialize(
  // graph_order(current_graph) , dim);
}

void TamiGraph::initialize(TamiBase::graph_type type) {

  graph_type = type;
  // TamiBase::graph_type type=AmiCalc::Sigma;

  // std::cout<< type <<std::endl;
  ami_parameters.TYPE_ = type; // type;

  // create_starter_graph(current_graph);

}




// Dev Note 08-22-2024: The starter graph construction was part of the graph generation process, but the necessary functions are not here and so for torchami we are marking these for removal.

// void TamiGraph::create_starter_graph(TamiGraph::graph_t &g) {

//   // std::cout<<"Creating starter graph"<<std::endl;
//   switch (ami_parameters.TYPE_) {

//   // typedef enum {Sigma,Pi_phuu, Pi_phud,Hartree, Bare, Greens, density,
//   // doubleocc} graph_type ;
//   case 0:
//     construct_starter_sigma(g);
//     break;
//   case 1:
//     construct_starter_phuu_bubble(g);
//     break;
//   case 2:
//     construct_starter_phud_bubble(g);
//     break;
//   case 3:
//     construct_starter_hartree(g);
//     break;
//   case 4:
//     construct_starter_bare(g);
//     break;
//   case 5:
//     construct_starter_bare(g);
//     break;
//   case 6:
//     construct_starter_bare(g);
//     break;
//   case 7:
//     construct_starter_phuu_bubble(g);
//     break;
//   case 8:
//     construct_starter_ppuu_bubble(g);
//     break;
//   case 9:
//     construct_starter_ppud_bubble(g);
//     break;
//   case 10:
//     construct_starter_bare(g);
//     break;
//   case 12:
//     construct_starter_force(g, f2, f3);
//     break;
//   default:
//     construct_starter_sigma(g);
//   }
// }

// void TamiGraph::construct_starter_sigma(TamiGraph::graph_t &g) {
//   // Set vertices and edges for graph to be 2nd order self energy

//   std::cout << "Constructing Sigma Starter Graph" << std::endl;
//   // Create 6 vertices in that graph
//   vertex_t in_vert = add_vertex(g);
//   vertex_t out_vert = add_vertex(g);

//   vertex_t first = add_vertex(g);
//   vertex_t last = add_vertex(g);

//   vertex_t bubble_left = add_vertex(g);
//   vertex_t bubble_right = add_vertex(g);

//   // Create an edge connecting those two vertices
//   add_edge(in_vert, first, edge_info(TamiBase::Fermi, 1, up), g);
//   add_edge(first, last, edge_info(TamiBase::Fermi, 1, up), g);
//   add_edge(last, out_vert, edge_info(TamiBase::Fermi, 1, up), g);

//   // make the bubble
//   add_edge(bubble_left, bubble_right, edge_info(TamiBase::Fermi, 2, dn), g);
//   add_edge(bubble_right, bubble_left, edge_info(TamiBase::Fermi, 2, dn), g);

//   // Add interaction lines that are Bose not Fermi
//   add_edge(first, bubble_left, edge_info(TamiBase::Bose), g);
//   add_edge(bubble_right, last, edge_info(TamiBase::Bose), g);

//   g[boost::graph_bundle].fermi_loop_id = 2;
// }

// void TamiGraph::construct_starter_ppuu_bubble(graph_t &g) {
//   // set vertices and edges for graph to be 0th order polarization bubble

//   std::cout << "Constructing Bubble Starter Graph" << std::endl;

//   vertex_t in_vert = add_vertex(g);
//   vertex_t out_vert = add_vertex(g);

//   vertex_t bubble_left = add_vertex(g);
//   vertex_t bubble_right = add_vertex(g);

//   // Add interaction lines that are Bose not Fermi
//   add_edge(in_vert, bubble_left, edge_info(TamiBase::Bose), g);
//   add_edge(bubble_right, out_vert, edge_info(TamiBase::Bose), g);

//   // make the bubble
//   add_edge(bubble_left, bubble_right, edge_info(TamiBase::Fermi, 1, up), g);
//   add_edge(bubble_left, bubble_right, edge_info(TamiBase::Fermi, 1, up), g);

//   g[boost::graph_bundle].fermi_loop_id = 1;
// }

// void TamiGraph::construct_starter_ppud_bubble(graph_t &g) {
//   // set vertices and edges for graph to be 0th order polarization bubble

//   std::cout << "Constructing Bubble Starter Graph" << std::endl;

//   vertex_t in_vert = add_vertex(g);
//   vertex_t out_vert = add_vertex(g);

//   vertex_t bubble_left = add_vertex(g);
//   vertex_t bubble_right = add_vertex(g);

//   // Add interaction lines that are Bose not Fermi
//   add_edge(in_vert, bubble_left, edge_info(TamiBase::Bose), g);
//   add_edge(bubble_right, out_vert, edge_info(TamiBase::Bose), g);

//   // make the bubble
//   add_edge(bubble_left, bubble_right, edge_info(TamiBase::Fermi, 1, up), g);
//   add_edge(bubble_left, bubble_right, edge_info(TamiBase::Fermi, 1, dn), g);

//   g[boost::graph_bundle].fermi_loop_id = 1;
// }

// void TamiGraph::construct_starter_hartree(TamiGraph::graph_t &g) {
//   // Set vertices and edges for graph to be 2nd order self energy

//   std::cout << "Constructing Hartree Starter Graph" << std::endl;
//   // Create 6 vertices in that graph
//   vertex_t in_vert = add_vertex(g);
//   vertex_t out_vert = add_vertex(g);

//   vertex_t bot = add_vertex(g);
//   vertex_t top = add_vertex(g);

//   // Create an edge connecting those two vertices
//   add_edge(in_vert, bot, edge_info(TamiBase::Fermi, 1, up), g);
//   add_edge(bot, out_vert, edge_info(TamiBase::Fermi, 1, up), g);

//   // make the int line
//   add_edge(bot, top, edge_info(TamiBase::Bose), g);

//   // Add loop
//   add_edge(top, top, edge_info(TamiBase::Fermi, 2, dn), g);

//   g[boost::graph_bundle].fermi_loop_id = 2;
// }

// void TamiGraph::construct_starter_bare(graph_t &g) {

//   std::cout << "Constructing Bare Starter Graph" << std::endl;

//   vertex_t in_vert = add_vertex(g);
//   vertex_t out_vert = add_vertex(g);

//   add_edge(in_vert, out_vert, edge_info(TamiBase::Fermi, 1, up), g);

//   g[boost::graph_bundle].fermi_loop_id = 1;
//   std::cout << "Complete" << std::endl;
// }

// void TamiGraph::construct_starter_phuu_bubble(graph_t &g) {
//   // set vertices and edges for graph to be 0th order polarization bubble

//   std::cout << "Constructing Bubble Starter Graph" << std::endl;

//   vertex_t in_vert = add_vertex(g);
//   vertex_t out_vert = add_vertex(g);

//   vertex_t bubble_left = add_vertex(g);
//   vertex_t bubble_right = add_vertex(g);

//   // Add interaction lines that are Bose not Fermi
//   add_edge(in_vert, bubble_left, edge_info(TamiBase::Bose), g);
//   add_edge(bubble_right, out_vert, edge_info(TamiBase::Bose), g);

//   // make the bubble
//   add_edge(bubble_left, bubble_right, edge_info(TamiBase::Fermi, 1, up), g);
//   add_edge(bubble_right, bubble_left, edge_info(TamiBase::Fermi, 1, up), g);

//   g[boost::graph_bundle].fermi_loop_id = 1;
// }

// void TamiGraph::construct_starter_phud_bubble(graph_t &g) {
//   // set vertices and edges for graph to be 0th order polarization bubble

//   std::cout << "Constructing Bubble Starter Graph" << std::endl;

//   vertex_t in_vert = add_vertex(g);
//   vertex_t out_vert = add_vertex(g);

//   vertex_t bubble_left = add_vertex(g);
//   vertex_t bubble_right = add_vertex(g);

//   // Add interaction lines that are Bose not Fermi
//   add_edge(in_vert, bubble_left, edge_info(TamiBase::Bose), g);
//   add_edge(bubble_right, out_vert, edge_info(TamiBase::Bose), g);

//   // make the bubble
//   add_edge(bubble_left, bubble_right, edge_info(TamiBase::Fermi, 1, up), g);
//   add_edge(bubble_right, bubble_left, edge_info(TamiBase::Fermi, 1, dn), g);

//   g[boost::graph_bundle].fermi_loop_id = 1;
// }

// void TamiGraph::construct_starter_force(graph_t &g, graph_t &f, graph_t &fp) {

//   std::cout << "Constructing force Starter Graph 1" << std::endl;

//   // vertex_t in_vert = add_vertex(g);
//   // vertex_t out_vert = add_vertex(g);

//   // vertex_t in_vert2 = add_vertex(g);
//   // vertex_t out_vert2 = add_vertex(g);

//   vertex_t bubble_left = add_vertex(g);
//   vertex_t bubble_right = add_vertex(g);

//   vertex_t bubble_left2 = add_vertex(g);
//   vertex_t bubble_right2 = add_vertex(g);

//   // Add force probing lines
//   // add_edge(in_vert,bubble_left, edge_info(TamiBase::Bose),g);
//   // add_edge(in_vert2,bubble_left2, edge_info(TamiBase::Bose),g);

//   // add_edge(bubble_right, out_vert, edge_info(TamiBase::Bose),g);
//   // add_edge(bubble_right2, out_vert2, edge_info(TamiBase::Bose),g);

//   // connect the bubbles
//   add_edge(bubble_left, bubble_left2, edge_info(TamiBase::Bose), g);
//   add_edge(bubble_right, bubble_right2, edge_info(TamiBase::Bose), g);

//   // make the bubble
//   add_edge(bubble_left, bubble_right, edge_info(TamiBase::Fermi, 1, up), g);
//   add_edge(bubble_right, bubble_left, edge_info(TamiBase::Fermi, 1, up), g);

//   // make the bubble
//   add_edge(bubble_left2, bubble_right2, edge_info(TamiBase::Fermi, 1, dn), g);
//   add_edge(bubble_right2, bubble_left2, edge_info(TamiBase::Fermi, 1, dn), g);

//   g[boost::graph_bundle].fermi_loop_id = 1;

//   number_vertices(g);
//   print_all_edge_info(g);

//   std::cout << "Constructing force Starter Graph 2" << std::endl;

//   // vertex_t fin_vert = add_vertex(f);
//   // vertex_t fout_vert = add_vertex(f);

//   // vertex_t fin_vert2 = add_vertex(f);
//   // vertex_t fout_vert2 = add_vertex(f);

//   vertex_t fbubble_left = add_vertex(f);
//   vertex_t fbubble_right = add_vertex(f);

//   vertex_t ftp_left = add_vertex(f);
//   vertex_t ftp_right = add_vertex(f);

//   // Add force probing lines
//   // add_edge(fin_vert,fbubble_left, edge_info(TamiBase::Bose),f);
//   // add_edge(fin_vert2,fbubble_left, edge_info(TamiBase::Bose),f);

//   // add_edge(fbubble_right, fout_vert, edge_info(TamiBase::Bose),f);
//   // add_edge(fbubble_right,fout_vert2, edge_info(TamiBase::Bose),f);

//   // connect those to the bubbles
//   add_edge(ftp_left, fbubble_left, edge_info(TamiBase::Bose), f);
//   add_edge(fbubble_right, ftp_right, edge_info(TamiBase::Bose), f);

//   // make the bubble
//   add_edge(fbubble_left, fbubble_right, edge_info(TamiBase::Fermi, 1, up), f);
//   add_edge(fbubble_right, fbubble_left, edge_info(TamiBase::Fermi, 1, up), f);

//   // make the tps
//   add_edge(ftp_left, ftp_left, edge_info(TamiBase::Fermi, 1, dn), f);
//   add_edge(ftp_right, ftp_right, edge_info(TamiBase::Fermi, 1, dn), f);

//   f[boost::graph_bundle].fermi_loop_id = 1;

//   number_vertices(f);
//   print_all_edge_info(f);

//   ////////////

//   std::cout << "Constructing force Starter Graph 3" << std::endl;

//   // vertex_t fpin_vert = add_vertex(fp);
//   // vertex_t fpout_vert = add_vertex(fp);

//   // vertex_t fpin_vert2 = add_vertex(fp);
//   // vertex_t fpout_vert2 = add_vertex(fp);

//   vertex_t fpbubble_left = add_vertex(fp);
//   vertex_t fpbubble_right = add_vertex(fp);

//   vertex_t fptp_left = add_vertex(fp);
//   vertex_t fptp_right = add_vertex(fp);

//   // Add force probing lines
//   // add_edge(fin_vert,fbubble_left, edge_info(TamiBase::Bose),f);
//   // add_edge(fin_vert2,fbubble_left, edge_info(TamiBase::Bose),f);

//   // add_edge(fbubble_right, fout_vert, edge_info(TamiBase::Bose),f);
//   // add_edge(fbubble_right,fout_vert2, edge_info(TamiBase::Bose),f);

//   // connect those to the bubbles
//   add_edge(fptp_left, fpbubble_left, edge_info(TamiBase::Bose), fp);
//   add_edge(fpbubble_right, fptp_right, edge_info(TamiBase::Bose), fp);

//   // make the bubble
//   add_edge(fpbubble_left, fpbubble_right, edge_info(TamiBase::Fermi, 1, dn),
//            fp);
//   add_edge(fpbubble_right, fpbubble_left, edge_info(TamiBase::Fermi, 1, dn),
//            fp);

//   // make the tps
//   add_edge(fptp_left, fptp_left, edge_info(TamiBase::Fermi, 1, up), fp);
//   add_edge(fptp_right, fptp_right, edge_info(TamiBase::Fermi, 1, up), fp);

//   f[boost::graph_bundle].fermi_loop_id = 1;

//   number_vertices(fp);
//   print_all_edge_info(fp);
// }
