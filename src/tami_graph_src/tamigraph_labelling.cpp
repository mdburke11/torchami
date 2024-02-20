#include "tami_graph.hpp"

using namespace boost;

void TamiGraph::repeated_labelling(graph_t &g, bool &result) {
  int iter = 0;
  result = false;

  number_vertices(g);
  // std::cout<<"Num vertices is "<<num_vertices(g)<<std::endl;
  // std::cout<<"Num edges is "<<num_edges(g)<<std::endl;

  /*
  for(int i=0; i<10; i++){
          i++;
  label_systematic(g);
  check_momentum_conservation(g, result);
  if(result==true){
          iter+=i;
          break;
          }
  } */

  if (result == false) {
    // std::cout<<"Starting half-random labelling"<<std::endl;
    // for(int i=0; i<2000;i++){
    do {
      // label_systematic_redone(g);

      // std::cout<<"Next attempt"<<std::endl;

      label_half_random(g);
      check_momentum_conservation(g, result);
      edge_vector_t unlabelled;
      find_unlabelled_fermionic_edges(g, unlabelled);
      if (unlabelled.size() != 0) {
        result = false;
      }
    } while (result == false);
  }
  /*
  if(result==false){
          // std::cout<<"Starting half-random labelling"<<std::endl;
  for(int i=0; i<20000;i++){

  //label_systematic_redone(g);
  random_labelling(g,result);
  // check_momentum_conservation(g, result);
  if(result==true){
          iter+=i;
          break;
          }
  }
  } */

  fix_epsilons(g);

  if (result == false) {
    std::cerr << "Warning: Labelling Failed." << std::endl;
  }
  // std::cout<<"Labelling took iter="<<iter<<std::endl;
}

void TamiGraph::fix_epsilons(graph_t &g) {

  edge_vector_t internal_fermi;

  find_internal_fermionic_edges(g, internal_fermi);

  for (int i = 0; i < internal_fermi.size(); i++) {
    for (int j = i; j < internal_fermi.size(); j++) {

      if (i != j) {

        bool check =
            edge_alphas_are_equal(internal_fermi[i], internal_fermi[j], g);
        bool check2 =
            edge_alphas_are_negative(internal_fermi[i], internal_fermi[j], g);
        if (check || check2) {
          g[internal_fermi[j]].g_struct_.eps_ =
              g[internal_fermi[i]].g_struct_.eps_;
        }
      }
    }
  }
}

void TamiGraph::find_unlabelled_edges(graph_t &g, edge_vector_t &vector) {

  vector.clear();

  ////std::cout<< "Finding Fermionic edges" <<std::endl;
  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  // Looking through all edges to find source and targets : this could be
  // improved?
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

    if (g[*ei].label == unlabelled) {
      vector.push_back(*ei);
    }
    ////std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_
    ///<< std::endl;

    ////std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<<
    ///total<<std::endl;
  }
}

/* void TamiGraph::find_entry_vertices(graph_t &g, vertex_vector_t &v){
v.clear();

boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
for (boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi){

        if (degree(*vi,g)==1 && in_degree(*vi,g)==0){ v=*vi;}

        }

}*/

void TamiGraph::label_half_random(graph_t &g) {

  reset_g(g);
  // initialize some counters on the interior of the graph
  g[boost::graph_bundle].n_indep = 0;
  g[boost::graph_bundle].n_labelled = 0;

  // number_vertices(g);  // assume it is already numbered. - moved to repeated
  // labelling at the start

  // label the external legs
  vertex_vector_t extern_vect_list;
  edge_vector_t extern_edge_list;
  find_external_vertices(g, extern_vect_list, extern_edge_list);

  // if(ami_parameters.TYPE_!=TamiBase::density && ami_parameters.TYPE_!=
  // TamiBase::doubleocc){
  if (extern_vect_list.size() < 1) {
    throw std::runtime_error("Can't find external leg.");
  }
  // }

  label_extern_legs(extern_edge_list, g);

  // Need function to find and label tadpoles
  vertex_vector_t tp_vert, tp_conn_vert;
  edge_vector_t tp_bose_edges;
  label_and_find_tadpoles_ami(g, tp_vert, tp_conn_vert, tp_bose_edges);

  // now get all the internal vertices

  vertex_vector_t internal;
  find_internal_vertices(g, internal);

  // fermionic loop stuff vs label along line stuff
  vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
  edge_vector_t labelled_edges, unlabelled_edges;

  do {
    int rand_choice = random_int(0, internal.size() - 1);
    sort_labelled_unlabelled_adjacent(
        g, internal[rand_choice], labelled_adj_vertices,
        unlabelled_adj_vertices, labelled_edges, unlabelled_edges);
    //
    // std::cout<<"Currently at "<< g[internal[rand_choice]].index_<<std::endl;
    // std::cout<<"Of list with length "<< internal.size()<<" chose random entry
    // "<< rand_choice<<std::endl; std::cout<<"Labelled are "<<
    // labelled_edges.size()<<"unlabelled are "<<
    // unlabelled_edges.size()<<std::endl;
    bool junk = label_consv_momentum(
        g, internal[rand_choice], labelled_adj_vertices,
        unlabelled_adj_vertices, labelled_edges, unlabelled_edges);

    internal.erase(internal.begin() + rand_choice);

  } while (internal.size() > 0);
}

void TamiGraph::find_unlabelled_fermionic_edges(graph_t &g,
                                                edge_vector_t &vector) {

  vector.clear();

  ////std::cout<< "Finding Fermionic edges" <<std::endl;
  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  // Looking through all edges to find source and targets : this could be
  // improved?
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

    if (g[*ei].g_struct_.stat_ == TamiBase::Fermi &&
        g[*ei].label == unlabelled) {
      vector.push_back(*ei);
    }
    ////std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_
    ///<< std::endl;

    ////std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<<
    ///total<<std::endl;
  }
}

void TamiGraph::reset_g(graph_t &g) {

  // reset labelled counters
  g[boost::graph_bundle].n_indep = 0;
  g[boost::graph_bundle].n_labelled = 0;

  edge_vector_t bose_edges;
  edge_vector_t fermi_edges;

  // edge_t bose_edge, fermi_edge;
  // find_one_fermi_bose_edge(g, bose_edge, fermi_edge);
  find_bose_fermi_edges(g, bose_edges, fermi_edges);
  // question is this internal only?

  // find_bosonic_edges(g, bose_edges);
  // find_fermionic_edges(g, fermi_edges);

  // std::vector<int> test1;

  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  // Looking through all edges to find source and targets : this could be
  // improved?
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

    // if( g[*ei].g_struct_.stat_==TamiBase::Fermi){  /// DONT need THIS LOGIC?
    // Bosonic lines get same reset of labels.
    //  std::vector<int> v{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    // std::fill(v.begin(), v.end(), -1);
    // test1.resize(bose_edges.size()+1);
    // std::fill(test1.begin, test1.end,0);

    if (graph_type == TamiBase::Pi_phuu || graph_type == TamiBase::Pi_phud ||
        (graph_type == TamiBase::doubleocc) ||
        ami_parameters.TYPE_ == TamiBase::Pi_ppud ||
        ami_parameters.TYPE_ == TamiBase::Pi_ppuu) {

      // std::cout<<"Bose edges are "<<bose_edges.size()<<std::endl;
      // std::cout<<"Fermi edges are "<<fermi_edges.size()<<std::endl;

      g[*ei].g_struct_.alpha_.resize(bose_edges.size());
      std::fill(g[*ei].g_struct_.alpha_.begin(), g[*ei].g_struct_.alpha_.end(),
                0);
      // TODO: figure out ACTUALLY how big eps vectors are because this is not
      // correct
      g[*ei].g_struct_.eps_.resize(fermi_edges.size());
      std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

      g[*ei].label = unlabelled;

    } else if (graph_type == TamiBase::density ||
               graph_type == TamiBase::Greens || graph_type == TamiBase::DOS ||
               graph_type == TamiBase::ENERGY) {

      if (num_edges(g) == 1) {

        g[*ei].g_struct_.alpha_.resize(1);
        g[*ei].g_struct_.eps_.resize(1);
      } else {

        g[*ei].g_struct_.alpha_.resize(bose_edges.size() + 1);
        g[*ei].g_struct_.eps_.resize(fermi_edges.size() - 1);
      }

      std::fill(g[*ei].g_struct_.alpha_.begin(), g[*ei].g_struct_.alpha_.end(),
                0);
      std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

      g[*ei].label = unlabelled;

    } else if (graph_type == TamiBase::FORCE) {

      int this_ord = graph_order(g) + 2;
      g[*ei].g_struct_.alpha_.resize(this_ord);
      std::fill(g[*ei].g_struct_.alpha_.begin(), g[*ei].g_struct_.alpha_.end(),
                0);
      // TODO: figure out ACTUALLY how big eps vectors are because this is not
      // correct
      g[*ei].g_struct_.eps_.resize(fermi_edges.size());
      std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

      g[*ei].label = unlabelled;

    }
    // else
    // if(graph_type== AmiCalc::density){

    // g[*ei].g_struct_.alpha_.resize(bose_edges.size()+1);
    // std::fill(g[*ei].g_struct_.alpha_.begin(),g[*ei].g_struct_.alpha_.end(),
    // 0);
    ////TODO: figure out ACTUALLY how big eps vectors are because this is not
    ///correct
    // g[*ei].g_struct_.eps_.resize(fermi_edges.size()-3);
    // std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

    // g[*ei].label=unlabelled;

    // }
    else {

      g[*ei].g_struct_.alpha_.resize(bose_edges.size() + 1);
      std::fill(g[*ei].g_struct_.alpha_.begin(), g[*ei].g_struct_.alpha_.end(),
                0);
      // TODO: figure out ACTUALLY how big eps vectors are because this is not
      // correct
      g[*ei].g_struct_.eps_.resize(fermi_edges.size() - 2);
      std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

      g[*ei].label = unlabelled;
    }
  }
}

// TODO: Issue. what if you give this a graph where the alphas and epsilons are
// empty? This won't work then.
void TamiGraph::label_extern_legs(edge_vector_t &extern_vect_list, graph_t &g) {

  // if(graph_type != TamiBase::Bose){
  edge_t e1;

  // std::cout<<"Supposed to be labelling legs here"<<std::endl;

  for (int i = 0; i < extern_vect_list.size(); i++) {

    // std::cout<<"Supposed to be labelling legs here at i="<<i<<std::endl;
    e1 = extern_vect_list[i];
    g[e1].g_struct_.alpha_[g[e1].g_struct_.alpha_.size() - 1] = 1;

    // TODO these lines are hardcoded for single external lines
    if (graph_type == TamiBase::density || graph_type == TamiBase::DOS ||
        graph_type == TamiBase::Greens || graph_type == TamiBase::ENERGY) {
      g[e1].g_struct_.eps_[g[e1].g_struct_.eps_.size() - 1] = 1;
    }

    // external legs have no epsilon value??? we set the last epsilon for the
    // external legs

    g[e1].label = labelled;
  }

  // }
  // else{

  // edge_t e1;

  // std::cout<<"Supposed to be labelling bose legs here"<<std::endl;

  // e1=extern_vect_list[i];
  // g[e1].g_struct_.alpha_[g[e1].g_struct_.alpha_.size()-1]=1;

  // }

  // g[boost::graph_bundle].n_labelled++;// we don't increment for external legs

  // TODO this is not general - since they all get the same label if there are
  // not two external legs this won't be right.
}

void TamiGraph::label_and_find_tadpoles_ami(graph_t &g, vertex_vector_t &tp_vec,
                                            vertex_vector_t &tp_conn_vec,
                                            edge_vector_t &tp_bose_edges) {

  //	vertex_vector_t tp_vert_vec;
  //	edge_vector_t tp_bose_line_vec;
  //	vertex_vector_t tp_connect_vert_vec;

  edge_vector_t edges_a, edges_b;

  find_tadpoles(g, tp_vec, tp_conn_vec, tp_bose_edges, edges_a, edges_b);

  // std::cout<<"Tadpoles and edges found "<<tp_conn_vec.size()<<"
  // "<<edges_a.size()<<" "<<edges_b.size()<<std::endl;

  if (edges_a.size() != tp_conn_vec.size() ||
      edges_b.size() != tp_conn_vec.size()) {
    throw std::runtime_error("Did not find edges_a and edges_b correctly");
  }

  edge_t tp_edge;
  for (int i = 0; i < tp_vec.size(); i++) {

    tp_edge = edge(tp_vec[i], tp_vec[i], g).first;

    label_indep_edge(tp_edge, g);

    g[tp_bose_edges[i]].label = labelled;

    // logic to label edges connected to tadpoles

    // if(g[edges_a[i]].label==unlabelled && g[edges_b[i]].label==unlabelled){
    // std::cout<<"Logic a triggered"<<std::endl;
    // label_indep_edge(edges_a[i],g);
    // assign_cons_label( g, tp_conn_vec[i], tp_bose_edges[i], edges_a[i],
    // edges_b[i]);

    // }

    if (g[edges_a[i]].label == unlabelled && g[edges_b[i]].label == labelled) {
      // std::cout<<"Logic b triggered"<<std::endl;
      // label_indep_edge(edges_a[i],g);
      assign_cons_label(g, tp_conn_vec[i], tp_bose_edges[i], edges_b[i],
                        edges_a[i]);
    }

    if (g[edges_a[i]].label == labelled && g[edges_b[i]].label == unlabelled) {
      // std::cout<<"Logic c triggered"<<std::endl;
      // label_indep_edge(edges_a[i],g);
      assign_cons_label(g, tp_conn_vec[i], tp_bose_edges[i], edges_a[i],
                        edges_b[i]);
    }
  }
}

void TamiGraph::label_indep_edge(edge_t &e, graph_t &g) {

  // std::cout<<"Attempting to add to alpha, and independent ends in element "<<
  // g[boost::graph_bundle].n_indep<<std::endl; std::cout<<"Alpha size is
  // "<<g[e].g_struct_.alpha_.size()<<std::endl;

  // std::cout<<"Attempting to add to eps, in element "<<
  // g[boost::graph_bundle].n_labelled<<std::endl;
  // std::cout<<"eps size is "<<g[e].g_struct_.eps_.size()<<std::endl;

  if (g[e].g_struct_.alpha_.size() <= g[boost::graph_bundle].n_indep ||
      g[e].g_struct_.eps_.size() <= g[boost::graph_bundle].n_labelled) {
    return;
  }
  g[e].g_struct_.alpha_[g[boost::graph_bundle].n_indep] = 1;
  g[e].g_struct_.eps_[g[boost::graph_bundle].n_labelled] = 1;

  // update tracking integers
  g[e].label = labelled;
  g[boost::graph_bundle].n_indep++;
  g[boost::graph_bundle].n_labelled++;
}

// This function takes two labelled fermi lines at a three point vertex and
// assigns a label to the bosonic line this has assigned only the momentum
// conserving alpha. the epsilon is just like the others.
void TamiGraph::assign_cons_label(graph_t &g, vertex_t &vin, edge_t &fermi_one,
                                  edge_t &fermi_two, edge_t &bose) {

  // std::cout<<print_edge(fermi_one,g)<<" "<<print_edge(fermi_two,g)<<"
  // "<<print_edge(bose,g)<<" "<<std::endl;

  if (fermi_one == fermi_two) {

    mark_labelled(bose, g);
    return;
  }

  if (source(bose, g) == vin) {

    if (target(fermi_one, g) == vin) {
      // std::cout<<"Adding"<<std::endl;
      add_labels(bose, fermi_one, g);
    } else {
      subtract_labels(bose, fermi_one, g);
    }

    if (target(fermi_two, g) == vin) {
      add_labels(bose, fermi_two, g);
      // std::cout<<"Adding"<<std::endl;
    } else {
      subtract_labels(bose, fermi_two, g);
    }

  } else // don't need the following if but for debugging i'm going to put it.
    if (target(bose, g) == vin) {

      if (target(fermi_one, g) == vin) {
        subtract_labels(bose, fermi_one, g);
      } else {
        add_labels(bose, fermi_one, g);
        // std::cout<<"Adding"<<std::endl;
      }

      if (target(fermi_two, g) == vin) {
        subtract_labels(bose, fermi_two, g);
      } else {
        add_labels(bose, fermi_two, g);
        // std::cout<<"Adding"<<std::endl;
      }
    }

  if (g[bose].g_struct_.stat_ == TamiBase::Fermi) {
    label_indep_epsilon(bose, g);
  } else {
    mark_labelled(bose, g);
  }
}

void TamiGraph::check_momentum_conservation(graph_t &g, bool &result) {

  bool cons = true;
  bool vert_cons = true;

  boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

  boost::graph_traits<graph_t>::in_edge_iterator iei, iei_end;
  boost::graph_traits<graph_t>::out_edge_iterator oei, oei_end;

  std::vector<int> sum_vec(g[random_edge(g, rand_gen)].g_struct_.alpha_.size(),
                           0);

  // std::cout<<std::endl;
  // std::cout<<"Checking momentum conservation"<<std::endl;

  for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
    vert_cons = true;
    // std::cout<<"On vertex "<<g[*vi].index_ <<std::endl;
    std::fill(sum_vec.begin(), sum_vec.end(), 0);

    if (degree(*vi, g) != 1) {
      // std::cout<<"On vertex "<<g[*vi].index_ << " checking
      // conserv"<<std::endl;

      for (boost::tie(iei, iei_end) = in_edges(*vi, g); iei != iei_end; ++iei) {
        for (int i = 0; i < g[*iei].g_struct_.alpha_.size(); i++) {
          //	//std::cout<<"in edges i= "<<i<<std::endl;
          sum_vec[i] += g[*iei].g_struct_.alpha_[i];
        }
      }

      for (boost::tie(oei, oei_end) = out_edges(*vi, g); oei != oei_end;
           ++oei) {

        for (int i = 0; i < g[*oei].g_struct_.alpha_.size(); i++) {
          //	//std::cout<<"in edges i= "<<i<<std::endl;
          sum_vec[i] -= g[*oei].g_struct_.alpha_[i];
        }
      }
    }

    for (int i = 0; i < sum_vec.size(); i++) {

      if (sum_vec[i] != 0) {
        // std::cout<<"Momentum not conserved at vertex "<< g[*vi].index_
        // <<std::endl;

        cons = false;
        vert_cons = false;
      }
    }
    // if(vert_cons==false){
    // print_edgeson_vert(g[*vi].index_);
    // }
  }

  // if(cons==true){std::cout<<"Momentum was conserved everywhere "<<std::endl;}
  // if(cons==false){std::cout<<"Momentum was NOT conserved  "<<std::endl;}

  result = cons;
}

void TamiGraph::sys_label(graph_t &g, bool &result) {

  // std::cout<<"Systematic bubble labelling for graph"<<std::endl;
  // print_all_edge_info(g);

  graph_t gc;
  gc = g;
  int ord;

  result = false;

  ord = graph_order(gc);
  // std::cout<<"Ord is "<<ord<<std::endl;
  // std::cout<<"Graph type is "<<graph_type<<std::endl;
  // ord is one bigger

  if ((graph_type == TamiBase::Pi_phuu) || (graph_type == TamiBase::Pi_phud) ||
      (graph_type == TamiBase::doubleocc) || graph_type == TamiBase::Pi_ppud ||
      graph_type == TamiBase::Pi_ppuu || graph_type == TamiBase::FORCE) {
    // std::cout<<"Increasing order by 1"<<std::endl;
    ord = graph_order(gc) + 1;
  }

  bool has_external_legs = true;
  // if a force graph remove the external legs completely.

  vertex_vector_t in_vv, out_vv;
  if (graph_type == TamiBase::FORCE) {

    find_force_LR_vertices(gc, in_vv, out_vv);

    vertex_vector_t extern_vect_list;
    edge_vector_t extern_edge_list;
    find_external_vertices(gc, extern_vect_list, extern_edge_list);
    delete_legs(gc, extern_vect_list, extern_edge_list);
    has_external_legs = false;
  }

  // if( graph_type==TamiBase::density){
  // ord=graph_order(gc)+1;
  // }
  // if(graph_type==TamiBase::doubleocc){
  // ord=graph_order(gc)+2;
  // }

  // std::cout<<"Ord is "<<ord<<std::endl;

  std::vector<int> v;

  reset_g(gc);

  // print_all_edge_info(gc);

  // if(ord==1){

  // }

  // label_indep_edge(internal_F_edges[v[j]],gc);

  // }

  // std::cout<<"Graph was reset and is now"<<std::endl;
  // print_all_edge_info(gc);

  // label the external legs
  vertex_vector_t extern_vect_list;
  edge_vector_t extern_edge_list;

  if (has_external_legs) {
    // std::cout<<"Looking for external "<<std::endl;
    find_external_vertices(gc, extern_vect_list, extern_edge_list);
    // std::cout<<"exited"<<std::endl;
    // std::cout<<"Extern edge list of size
    // "<<extern_vect_list.size()<<std::endl;
    // if(ami_parameters.TYPE_!=TamiBase::density && ami_parameters.TYPE_!=
    // TamiBase::doubleocc){
    if (extern_vect_list.size() < 1) {
      throw std::runtime_error("Can't find external leg.");
    }
  }
  // }
  // if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find
  // external leg.");}

  // TODO: I think this is redundant
  // label_extern_legs(extern_edge_list,gc);

  // std::cout<<std::endl<<std::endl<<"-----external
  // labelled-------"<<std::endl; print_all_edge_info(gc);

  ///

  edge_vector_t internal_F_edges;
  find_internal_fermionic_edges(gc, internal_F_edges);

  // std::cout<<"-----Prelabelled------"<<std::endl;
  // print_all_edge_info(gc);
  // std::cout<<"-----------"<<std::endl;

  // std::cout<<"Number of internal edges is "<<
  // internal_F_edges.size()<<std::endl;

  vertex_vector_t internal_v;
  find_internal_vertices(gc, internal_v);

  // std::cout<<"Number of internal vertices is "<<
  // internal_v.size()<<std::endl;

  // fermionic loop stuff vs label along line stuff
  vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
  edge_vector_t labelled_edges, unlabelled_edges;

  int count = 0;

  // take entire list of F edges.
  for (int i = 0; i < internal_F_edges.size(); i++) {
    v.push_back(i);
  }

  ///
  do {
    count++;
    reset_g(gc);

    // std::cout<<"count is "<<count<<std::endl;

    // std::cout<<"-----Prelabelled------"<<std::endl;
    // print_all_edge_info(gc);
    // std::cout<<"-----------"<<std::endl;

    if (has_external_legs) {
      label_extern_legs(extern_edge_list, gc);
    }
    // std::cout<<"-----------"<<std::endl;
    // print_all_edge_info(gc);
    // std::cout<<"-----------"<<std::endl;
    // initialize some counters on the interior of the graph
    gc[boost::graph_bundle].n_indep = 0;
    gc[boost::graph_bundle].n_labelled = 0;

    // 2022-03-25 : Added tadpole labelling
    // Need function to find and label tadpoles
    // vertex_vector_t tp_vert, tp_conn_vert;
    // edge_vector_t tp_bose_edges;
    // label_and_find_tadpoles_ami(gc,tp_vert, tp_conn_vert, tp_bose_edges);

    // vertex_vector_t tp_vec;
    // edge_vector_t tp_bose_edges;
    // vertex_vector_t tp_conn_vec;

    // edge_vector_t edges_a, edges_b;

    // find_tadpoles(gc, tp_vec, tp_conn_vec, tp_bose_edges, edges_a, edges_b);
    // std::cout<<"Found number of tadpoles="<<tp_bose_edges.size()<<std::endl;
    // for(int m=0; m<tp_bose_edges.size(); m++){
    // mark_labelled(tp_bose_edges[m],gc);
    // }
    // vertex_vector_t tp_vert, tp_conn_vert;
    // edge_vector_t tp_bose_edges;

    // print_all_edge_info(gc);

    // find_and_label_tadpole_strings(gc,tp_vert, tp_conn_vert, tp_bose_edges);

    // std::cout<<"After tadpole labelling "<<std::endl;
    // print_all_edge_info(gc);
    // exit(0);

    // how is this even systematic?
    for (int j = 0; j < ord; j++) {

      label_indep_edge(internal_F_edges[v[j]], gc);
      // std::cout<<"Labelling independent edge "<< v[j]<<std::endl;
      // print_all_edge_info(gc);
    }

    // std::cout<<"Indep edges labelled"<<std::endl;
    // print_all_edge_info(gc);
    // now we have the independent edges labelled. So now we iterate through
    // vertices, looking for a vertex with two edges labelled

    bool cont = true;
    edge_vector_t in, out;
    int first, second;

    // print_all_edge_info(gc);

    do {
      find_unlabelled_edges(gc, in);
      first = in.size();

      for (int i = 0; i < internal_v.size(); i++) {

        // get labelled and unlabelled edges attached to the i'th vertex
        sort_labelled_unlabelled_adjacent(
            gc, internal_v[i], labelled_adj_vertices, unlabelled_adj_vertices,
            labelled_edges, unlabelled_edges);

        if (unlabelled_edges.size() != 1) {
          continue;
        } // immediately break function since we can't label.

        // Can only now do something if 2 of three edges are labelled.
        // if (unlabelled_edges.size()==1){

        if (labelled_edges.size() != 2) {

          continue; // immediatly break instead of error?
          // std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are
          // "<< unlabelled_edges.size()<<std::endl; print_all_edge_info(gc);
          // std::cout<<"In and out vertices are"<<std::endl;
          // std::cout<<gc[in_vv[0]].index_<<" "<<gc[in_vv[1]].index_<<" "<<
          // gc[out_vv[0]].index_<<" "<<gc[out_vv[1]].index_<<std::endl;
          throw std::runtime_error(
              "Syslabel: Number of labelled plus unlabelled not equal to total "
              "number? Can't be true :) ");
        }

        // std::cout<<"Assigning conserving label to edge
        // "<<print_edge(unlabelled_edges[0],g)<<  std::endl;
        // std::cout<<"Labelled edges are "<<print_edge(labelled_edges[0],g)<<"
        // "<<print_edge(labelled_edges[1],g)<<std::endl;
        assign_cons_label(gc, internal_v[i], labelled_edges[0],
                          labelled_edges[1], unlabelled_edges[0]);
        // std::cout<<"Done"<<std::endl;
        // }
      }

      find_unlabelled_edges(gc, out);
      second = out.size();

      // std::cout<<first<<" "<<second<<std::endl;

      if (first == second) {
        cont = false;
      }

    } while (cont);

    // at this point the graph cannot be labelled further. check two things

    // 1) Has every line been labelled? aka, out.size()==0
    // 2) is momentum conserved.

    // std::cout<<"out.size()="<<out.size()<<std::endl;
    bool append = false;
    find_unlabelled_edges(gc, out);
    if (out.size() == 0) {

      check_momentum_conservation(gc, append);
      // std::cout<<"Momentum conserved? "<<append<<std::endl;

      // print_all_edge_info(gc);

      if (append) {

        fix_epsilons(gc);
        // std::cout<<"Pushing into L"<<std::endl;

        if (graph_type == TamiBase::FORCE) {
          std::cout << "Fixing force labels" << std::endl;
          fix_force_labels(gc, in_vv, out_vv);
          put_back_legs(gc, in_vv, out_vv);
        }

        // std::cout<<"Final print before copy"<<std::endl;
        // print_all_edge_info(g);
        // std::cout<<std::endl;
        // print_all_edge_info(gc);

        g = gc;
        result = true;

        // std::cout<<"-----Postlabelled------"<<std::endl;
        // print_all_edge_info(gc);
        // std::cout<<"-----------"<<std::endl;

        return;
      }
    }
    // std::cout<<"one..";
    std::reverse(v.begin() + ord, v.end());
    // std::cout<<"two"<<std::endl;
  } while (std::next_permutation(v.begin(), v.end()));
}

/// This finds the vertices that are external - use with caution.
void TamiGraph::find_external_vertices(graph_t &g, vertex_vector_t &v,
                                       edge_vector_t &edges) {

  v.clear();
  edges.clear();

  boost::graph_traits<graph_t>::in_edge_iterator iei, iei_end;
  boost::graph_traits<graph_t>::out_edge_iterator oei, oei_end;

  boost::graph_traits<graph_t>::vertex_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = vertices(g); ei != ei_end; ++ei) {

    if (degree(*ei, g) == 1) {
      v.push_back(*ei);
      for (boost::tie(iei, iei_end) = in_edges(*ei, g); iei != iei_end; ++iei) {
        edges.push_back(*iei);
      }

      for (boost::tie(oei, oei_end) = out_edges(*ei, g); oei != oei_end;
           ++oei) {
        edges.push_back(*oei);
      }
    }
  }
}

void TamiGraph::find_tadpoles(graph_t &g, vertex_vector_t &tp_vec,
                              vertex_vector_t &tp_conn_vec,
                              edge_vector_t &tp_bose_edges,
                              edge_vector_t &edge_a, edge_vector_t &edge_b) {

  tp_vec.clear();
  tp_conn_vec.clear();
  tp_bose_edges.clear();
  edge_a.clear();
  edge_b.clear();

  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  boost::graph_traits<graph_t>::edge_iterator ei2, ei_end2;
  edge_t bose;

  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

    if (source(*ei, g) == target(*ei, g)) {
      tp_vec.push_back(source(*ei, g));
      find_bosonic_edge_on_three_pointvert(g, source(*ei, g), bose);
      tp_bose_edges.push_back(bose);

      // std::cout<<"found bose edge from "<<
      // g[source(bose,g)].index_<<"-"<<g[target(bose,g)].index_<<std::endl;

      if (source(bose, g) == source(*ei, g)) {
        tp_conn_vec.push_back(target(bose, g));
      }
      if (target(bose, g) == source(*ei, g)) {
        tp_conn_vec.push_back(source(bose, g));
      }

      for (boost::tie(ei2, ei_end2) = edges(g); ei2 != ei_end2; ++ei2) {

        if (g[*ei2].g_struct_.stat_ == TamiBase::Fermi) {

          if (target(*ei2, g) == tp_conn_vec.back()) {
            edge_a.push_back(*ei2);
          } // std::cout<<"pushback for a"<< std::endl;}
          if (source(*ei2, g) == tp_conn_vec.back()) {
            edge_b.push_back(*ei2);
          } // std::cout<<"pushback for b"<< std::endl;}
        }
      }
    }
  }
}

void TamiGraph::delete_legs(graph_t &g, vertex_vector_t &v,
                            edge_vector_t &edges) {

  for (int i = 0; i < v.size(); i++) {
    clear_vertex(v[i], g);
    remove_vertex(v[i], g);
  }
}

bool TamiGraph::label_consv_momentum(graph_t &g, vertex_t &vin,
                                     vertex_vector_t &labelled_adj_vertices,
                                     vertex_vector_t &unlabelled_adj_vertices,
                                     edge_vector_t &labelled_edges,
                                     edge_vector_t &unlabelled_edges) {

  edge_t bose;
  bool output = true;
  if (unlabelled_edges.size() == 3) {
    // std::cout<<"found 3 unlabelled edges"<<std::endl;
    output = false;
    return output;
  } // immediately break function since we can't label.
  if (unlabelled_edges.size() == 0) {
    // std::cout<<"no unlabelled edges"<<std::endl;
    output = false;
    return output;
  } // immediately break function since nothing to label.

  if (unlabelled_edges.size() == 2) {
    // sanity check
    if (labelled_edges.size() != 1) {
      std::cout << "Labelled are " << labelled_edges.size() << "unlabelled are "
                << unlabelled_edges.size() << std::endl;
      print_all_edge_info(g);
      throw std::runtime_error(
          "In Labelling: Number of labelled plus unlabelled not equal to total "
          "number? Can't be true :) ");
    }

    // both fermionic lines
    if (g[unlabelled_edges[0]].g_struct_.stat_ == TamiBase::Fermi &&
        g[unlabelled_edges[1]].g_struct_.stat_ == TamiBase::Fermi) {

      // label one of the two fermi edges independently and set the other
      // according to consv of momentum

      label_indep_edge(unlabelled_edges[0], g);
      assign_cons_label(g, vin, labelled_edges[0], unlabelled_edges[0],
                        unlabelled_edges[1]); // 0 is now labelled
    }

    // first fermi second bose
    if (g[unlabelled_edges[0]].g_struct_.stat_ == TamiBase::Fermi &&
        g[unlabelled_edges[1]].g_struct_.stat_ == TamiBase::Bose) {

      // bose=unlabelled_edges[1];
      // if(is_labelled_bose_zero(bose,g)){ std::cout<<"Found false
      // logic"<<std::endl; add_labels(unlabelled_edges[0],labelled_edges[0],g);
      // label_indep_epsilon(unlabelled_edges[0],g);
      // g[unlabelled_edges[0]].label=labelled;
      // }else{

      label_indep_edge(unlabelled_edges[0], g);
      assign_cons_label(g, vin, labelled_edges[0], unlabelled_edges[0],
                        unlabelled_edges[1]);
      // }
    }

    // first bose second fermi
    if (g[unlabelled_edges[0]].g_struct_.stat_ == TamiBase::Bose &&
        g[unlabelled_edges[1]].g_struct_.stat_ == TamiBase::Fermi) {

      // bose=unlabelled_edges[0];
      // if(is_labelled_bose_zero(bose,g)){std::cout<<"Found false
      // logic"<<std::endl; add_labels(unlabelled_edges[1],labelled_edges[0],g);
      // label_indep_epsilon(unlabelled_edges[1],g);
      // g[unlabelled_edges[0]].label=labelled;
      // }else{
      label_indep_edge(unlabelled_edges[1], g);
      assign_cons_label(g, vin, labelled_edges[0], unlabelled_edges[1],
                        unlabelled_edges[0]);
      // }
    }
  }

  if (unlabelled_edges.size() == 1) {

    if (labelled_edges.size() != 2) {
      // std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<<
      // unlabelled_edges.size()<<std::endl;
      throw std::runtime_error(
          "In Label_consv: Number of labelled plus unlabelled not equal to "
          "total number? Can't be true :) ");
    }

    assign_cons_label(g, vin, labelled_edges[0], labelled_edges[1],
                      unlabelled_edges[0]);
  }

  return output;
}

// this function is redundant, since it needs the same operations as the next
// part of the step... the order matters
void TamiGraph::get_adjacent_vertex_stat(graph_t &g, vertex_t &vin,
                                         vertex_t &next, edge_t &eout,
                                         TamiBase::stat_type stat) {

  edge_vector_t edges;

  boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

  for (boost::tie(ai, ai_end) = adjacent_vertices(vin, g); ai != ai_end; ++ai) {

    if (edge(vin, *ai, g).second) {
      if (g[edge(vin, *ai, g).first].g_struct_.stat_ == stat) {

        next = *ai;
        eout = edge(vin, *ai, g).first;
      }
    }
  }
}

void TamiGraph::sort_labelled_unlabelled_adjacent(
    graph_t &g, vertex_t &vin, vertex_vector_t &labelled_adj_vertices,
    vertex_vector_t &unlabelled_adj_vertices, edge_vector_t &labelled_edges,
    edge_vector_t &unlabelled_edges) {

  ////std::cout<<"sorting labelled/unlabelled"<<std::endl;

  // reset vectors of labelled/unlabelled sorting
  labelled_adj_vertices.clear();
  unlabelled_adj_vertices.clear();
  labelled_edges.clear();
  unlabelled_edges.clear();

  // in principle this should always be true
  // if(g[edge(last,vin,g).first].label==labelled){
  // labelled_adj_vertices.push_back(last);
  // labelled_edges.push_back(edge(last,vin,g).first);
  // }else{
  // unlabelled_adj_vertices.push_back(last);
  // unlabelled_edges.push_back(edge(last,vin,g).first);
  // }

  boost::graph_traits<graph_t>::in_edge_iterator iei, iei_end;
  boost::graph_traits<graph_t>::out_edge_iterator oei, oei_end;

  for (boost::tie(iei, iei_end) = in_edges(vin, g); iei != iei_end; ++iei) {

    if (g[*iei].label == labelled) {
      labelled_adj_vertices.push_back(source(*iei, g));
      labelled_edges.push_back(*iei);
    } else {
      unlabelled_adj_vertices.push_back(source(*iei, g));
      unlabelled_edges.push_back(*iei);
    }
  }

  for (boost::tie(oei, oei_end) = out_edges(vin, g); oei != oei_end; ++oei) {

    if (g[*oei].label == labelled) {
      labelled_adj_vertices.push_back(target(*oei, g));
      labelled_edges.push_back(*oei);
    } else {
      unlabelled_adj_vertices.push_back(target(*oei, g));
      unlabelled_edges.push_back(*oei);
    }
  }

  // boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

  // for (boost::tie(ai,ai_end) = adjacent_vertices(vin, g); ai != ai_end;
  // ++ai){

  // if( edge(vin, *ai, g).second){

  // if(g[edge(vin,*ai,g).first].label==labelled){
  // labelled_adj_vertices.push_back(target(edge(vin,*ai,g).first,g));
  // labelled_edges.push_back(edge(vin,*ai,g).first);

  // }else{
  // unlabelled_adj_vertices.push_back(target(edge(vin,*ai,g).first,g));
  // unlabelled_edges.push_back(edge(vin,*ai,g).first);

  // }

  // }else
  // if (edge(*ai, vin, g).second){
  // if(g[edge(*ai, vin,g).first].label==labelled){
  // labelled_adj_vertices.push_back(source(edge(*ai, vin,g).first,g));
  // labelled_edges.push_back(edge(*ai, vin,g).first);
  // }else {
  // unlabelled_adj_vertices.push_back(source(edge(*ai, vin,g).first,g));
  // unlabelled_edges.push_back(edge(*ai, vin,g).first);

  // }

  // }
  // }

  // std::cout<<"done sorting labelled/unlabelled"<<std::endl;
  // std::cout<<"Labelled is "<<  labelled_adj_vertices.size()<<"  unLabelled is
  // "<<  unlabelled_adj_vertices.size()<<std::endl;
}

void TamiGraph::reset_epsilons(TamiBase::g_prod_t &R0) {

  int this_size = R0.size(); // g[internal[i]].g_struct_.eps_.size();
  for (int i = 0; i < R0.size(); i++) {

    // this resets the epsilons
    // int this_size=g[internal[i]].g_struct_.eps_.size();
    R0[i].eps_.clear();
    R0[i].eps_.resize(this_size, 0);
    R0[i].eps_[i] = 1;
  }

  return;
}

void TamiGraph::add_labels(edge_t &one, edge_t &two, graph_t &g) {

  for (int i = 0; i < g[one].g_struct_.alpha_.size(); i++) {
    g[one].g_struct_.alpha_[i] += g[two].g_struct_.alpha_[i];
  }
}

void TamiGraph::subtract_labels(edge_t &one, edge_t &two, graph_t &g) {

  for (int i = 0; i < g[one].g_struct_.alpha_.size(); i++) {
    g[one].g_struct_.alpha_[i] -= g[two].g_struct_.alpha_[i];
  }
}

void TamiGraph::mark_labelled(edge_t &e, graph_t &g) {

  // update tracking integers
  g[e].label = labelled;
}

void TamiGraph::label_indep_epsilon(edge_t &e, graph_t &g) {

  // std::cout<<"Attempting to add to eps, in element "<<
  // g[boost::graph_bundle].n_labelled<<std::endl; std::cout<<"eps size is
  // "<<g[e].g_struct_.eps_.size()<<std::endl;

  //	g[e].g_struct_.alpha_[g[boost::graph_bundle].n_indep]=1;
  g[e].g_struct_.eps_[g[boost::graph_bundle].n_labelled] = 1;

  // update tracking integers
  g[e].label = labelled;
  //	g[boost::graph_bundle].n_indep++;
  g[boost::graph_bundle].n_labelled++;
}

bool TamiGraph::edge_alphas_are_equal(edge_t &one, edge_t &two, graph_t &g) {

  for (int i = 0; i < g[one].g_struct_.alpha_.size(); i++) {
    if (g[one].g_struct_.alpha_[i] != g[two].g_struct_.alpha_[i]) {
      return false;
    }
  }

  return true;
}

bool TamiGraph::edge_alphas_are_negative(edge_t &one, edge_t &two, graph_t &g) {

  for (int i = 0; i < g[one].g_struct_.alpha_.size(); i++) {
    if (g[one].g_struct_.alpha_[i] != -g[two].g_struct_.alpha_[i]) {
      return false;
    }
  }

  return true;
}

void TamiGraph::find_internal_fermionic_edges(graph_t &g,
                                              edge_vector_t &vector) {

  vector.clear();

  // std::cout<< "Finding Fermionic edges" <<std::endl;
  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  // Looking through all edges to find source and targets : this could be
  // improved?
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

    if (g[*ei].g_struct_.stat_ == TamiBase::Fermi) {
      // if(source(*ei,g)==target(*ei,g)){vector.push_back(*ei);}
      if (degree(source(*ei, g), g) != 1 && degree(target(*ei, g), g) != 1) {
        vector.push_back(*ei);
      }
    }
    // std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_
    // << std::endl;

    // std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<<
    // total<<std::endl;
  }
}

void TamiGraph::ggm_remove_pairs(gg_matrix_t &ggm, int min) {

  for (int ord = min; ord < ggm.size(); ord++) {
    for (int group = 0; group < ggm[ord].size(); group++) {
      for (int pair = 0; pair < ggm[ord][group].gp_vec.size(); pair++) {

        ggm[ord][group].graph_vec.push_back(ggm[ord][group].gp_vec[pair].g1_);
        ggm[ord][group].graph_vec.push_back(ggm[ord][group].gp_vec[pair].g2_);
      }
      // after all pairs pushing into graph_vec. we delete the gp_vecs
      ggm[ord][group].gp_vec.clear();
    }
  }
}
