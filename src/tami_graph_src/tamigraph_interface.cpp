#include "tami_graph.hpp"
#include <ifstream>
#include <string>
#include <stringstream>

using namespace boost;

void TamiGraph::graph_to_R0(graph_t &g, TamiBase::g_prod_t &R0) {

  R0.clear();

  // std::cout<<"Converting graph to R0_"<<std::endl;

  // print_all_edge_info(g);

  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  // Step through all edges and collect the alpha and epsilon labels of the 'g'
  // graph
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

    // g[*ei].g_struct_.stat_==TamiBase::Bose){

    if (g[*ei].g_struct_.stat_ == TamiBase::Fermi) {

      // special case for bare lines - since boost graph 'degree' function
      // returns 1 if you close a line on itself.
      /* 	if(num_edges(g)==1 && graph_type==TamiBase::density){
                      std::cout<<"Triggered on "<<std::endl;
              g[*ei].g_struct_.species_=(TamiBase::species_t)g[*ei].spin;
              R0.push_back(g[*ei].g_struct_);
              } */

      // if(graph_type==TamiBase::ENERGY){

      // if ( !(degree( source(*ei,g),g)==1)){

      // g[*ei].g_struct_.species_=(TamiBase::species_t)g[*ei].spin;
      // std::cout<<"G had spin "<<g[*ei].spin<<std::endl;
      // R0.push_back(g[*ei].g_struct_);
      // }

      // } else

      if (graph_type == TamiBase::ENERGY) {

        if (num_edges(g) == 1) {
          g[*ei].g_struct_.species_ = (TamiBase::species_t)g[*ei].spin;
          R0.push_back(g[*ei].g_struct_);

        } else if (!(degree(source(*ei, g), g) == 1)) {
          g[*ei].g_struct_.species_ = (TamiBase::species_t)g[*ei].spin;
          R0.push_back(g[*ei].g_struct_);
        }
      }

      if (graph_type == TamiBase::density || graph_type == TamiBase::Greens ||
          graph_type == TamiBase::DOS) {

        g[*ei].g_struct_.species_ = (TamiBase::species_t)g[*ei].spin;
        // std::cout<<"G had spin "<<g[*ei].spin<<std::endl;
        R0.push_back(g[*ei].g_struct_);

      } else {

        if (!(degree(source(*ei, g), g) == 1 ||
              degree(target(*ei, g), g) == 1)) {
          // std::cout<<"Pushing back"<<std::endl;
          g[*ei].g_struct_.species_ = (TamiBase::species_t)g[*ei].spin;
          // std::cout<<"G had spin "<<g[*ei].spin<<std::endl;
          R0.push_back(g[*ei].g_struct_);
        }
      }
    }
  }

  if (bose_alphas_in_R0) {

    // put the bosonic lines in here
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

      if (g[*ei].g_struct_.stat_ == TamiBase::Bose) {

        if (!(degree(source(*ei, g), g) == 1 ||
              degree(target(*ei, g), g) == 1)) {

          g[*ei].g_struct_.species_ = (TamiBase::species_t)g[*ei].spin;
          g[*ei].g_struct_.eff_stat_ = TamiBase::Bose;
          // std::transform(g[*ei].g_struct_.alpha_.begin(),
          // g[*ei].g_struct_.alpha_.end(), g[*ei].g_struct_.alpha_.begin(),
          // std::bind(std::multiplies<int>(), std::placeholders::_1, -1));
          // g[*ei].g_struct_.alpha_=-g[*ei].g_struct_.alpha_;
          R0.push_back(g[*ei].g_struct_);
        }
      }
    }

    // now reset the epsilons to include the extra green's functions.

    reset_epsilons(R0);
  }

  // std::cout<<"Finished with size "<<R0.size()<<std::endl;

  // return R0;
}

// TODO: this won't work for bare phonon(bosonic) propagators
void TamiGraph::extract_bose_alphas(graph_t g,
                                    std::vector<TamiBase::alpha_t> &bose) {

  graph_t gc = g;

  // std::cout<<"Here with num_edges "<< num_edges(gc)<<std::endl;
  // std::cout<<std::flush;
  if (num_edges(gc) > 1) {

    bose.clear();

    vertex_vector_t extern_vert_list;
    edge_vector_t extern_edge_list;

    // std::cout<<"Finding externals"<<std::endl;
    find_external_vertices(gc, extern_vert_list, extern_edge_list);

    // std::cout<<"Delete legs"<<std::endl;
    delete_legs(g, extern_vert_list, extern_edge_list);

    boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = edges(gc); ei != ei_end; ++ei) {

      if (g[*ei].g_struct_.stat_ == TamiBase::Bose) {
        bose.push_back(g[*ei].g_struct_.alpha_);
      }
    }
  }
}

void TamiGraph::read_ggmp(std::string folder, gg_matrix_t &ggm, int max_ord) {

  ggm.clear();
  // oX_gX_nX.graph
  std::stringstream ss;
  ss << std::filesystem::current_path().string() << "/" << folder;
  std::string path = ss.str();

  if (!std::filesystem::is_directory(path)) {
    throw std::runtime_error("Could not open ggm directory");
  }

  for (const auto &entry : std::filesystem::directory_iterator(path)) {

    // std::cout << entry.path() << std::endl;
    // std::cout << entry.path().extension() <<"
    // "<<entry.path().extension().string()<<  std::endl;

    if (entry.path().extension().string().compare(".graph") == 0) {
      // std::cout<<"Was true group?
      // "<<entry.path().extension().string().compare(".group")<<std::endl;

      std::string stem = entry.path().stem().string();

      std::size_t pos = stem.find("_");
      std::string second = stem.substr(1, pos - std::size_t(1));
      int ord = std::stoi(second);
      stem = stem.substr(pos + std::size_t(1));

      pos = stem.find("_");
      second = stem.substr(1, pos - std::size_t(1));
      int group = std::stoi(second);

      second = stem.substr(pos + std::size_t(2));
      int num = std::stoi(second);

      if (ord > max_ord) {
        continue;
      }

      // std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

      // std::cout << entry.path().filename() << std::endl;
      // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num
      // <<std::endl;

      if (ggm.size() < ord + 1) {
        ggm.resize(ord + 1);
      }
      if (ggm[ord].size() < group + 1) {
        ggm[ord].resize(group + 1);
      }
      if (ggm[ord][group].graph_vec.size() < num + 1) {
        ggm[ord][group].graph_vec.resize(num + 1);
      }

      graph_t temp;
      graph_read(entry.path().string(), temp);

      ggm[ord][group].graph_vec[num] = temp;
    }

    // std::strcmp(entry.path().extension().string(),".group")
    // if(entry.path().extension().string().compare("pair")==0){
    // if(entry.path().extension().string()==".pair"){
    if (entry.path().extension().string().compare(".pair") == 0) {
      // std::cout<<"Was true pair?
      // "<<entry.path().extension().string().compare(".group")<<std::endl;

      std::string stem = entry.path().stem().string();

      std::size_t pos = stem.find("_");
      std::string second = stem.substr(1, pos - std::size_t(1));
      int ord = std::stoi(second);
      stem = stem.substr(pos + std::size_t(1));

      pos = stem.find("_");
      second = stem.substr(1, pos - std::size_t(1));
      int group = std::stoi(second);

      second = stem.substr(pos + std::size_t(2));
      int num = std::stoi(second);

      if (ord > max_ord) {
        continue;
      }

      // std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

      // std::cout << entry.path().filename() << std::endl;
      // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num
      // <<std::endl;

      if (ggm.size() < ord + 1) {
        ggm.resize(ord + 1);
      }
      if (ggm[ord].size() < group + 1) {
        ggm[ord].resize(group + 1);
      }
      if (ggm[ord][group].gp_vec.size() < num + 1) {
        ggm[ord][group].gp_vec.resize(num + 1);
      }

      git_pair temppair;
      pair_read(entry.path().string(), temppair);

      ggm[ord][group].gp_vec[num] = temppair;
    }
  }

  // for(int ord=2; ord< ggm.size(); ord++){
  // for(int group=0; group< ggm[ord].size(); group++){
  // for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++)
}

void TamiGraph::read_ggmp(std::string folder, gg_matrix_t &ggm, int min_ord,
                          int max_ord) {

  ggm.clear();
  // oX_gX_nX.graph
  std::stringstream ss;
  ss << std::filesystem::current_path().string() << "/" << folder;
  std::string path = ss.str();

  if (!std::filesystem::is_directory(path)) {
    throw std::runtime_error("Could not open ggm directory");
  }

  for (const auto &entry : std::filesystem::directory_iterator(path)) {

    // std::cout << entry.path() << std::endl;
    // std::cout << entry.path().extension() <<"
    // "<<entry.path().extension().string()<<  std::endl;

    if (entry.path().extension().string().compare(".graph") == 0) {
      // std::cout<<"Was true group?
      // "<<entry.path().extension().string().compare(".group")<<std::endl;

      std::string stem = entry.path().stem().string();

      std::size_t pos = stem.find("_");
      std::string second = stem.substr(1, pos - std::size_t(1));
      int ord = std::stoi(second);
      stem = stem.substr(pos + std::size_t(1));

      pos = stem.find("_");
      second = stem.substr(1, pos - std::size_t(1));
      int group = std::stoi(second);

      second = stem.substr(pos + std::size_t(2));
      int num = std::stoi(second);

      if (ord > max_ord) {
        continue;
      }

      // std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

      // std::cout << entry.path().filename() << std::endl;
      // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num
      // <<std::endl;

      if (ggm.size() < ord + 1) {
        ggm.resize(ord + 1);
      }
      if (ggm[ord].size() < group + 1) {
        ggm[ord].resize(group + 1);
      }
      if (ggm[ord][group].graph_vec.size() < num + 1) {
        ggm[ord][group].graph_vec.resize(num + 1);
      }

      graph_t temp;
      graph_read(entry.path().string(), temp);

      ggm[ord][group].graph_vec[num] = temp;
    }

    // std::strcmp(entry.path().extension().string(),".group")
    // if(entry.path().extension().string().compare("pair")==0){
    // if(entry.path().extension().string()==".pair"){
    if (entry.path().extension().string().compare(".pair") == 0) {
      // std::cout<<"Was true pair?
      // "<<entry.path().extension().string().compare(".group")<<std::endl;

      std::string stem = entry.path().stem().string();

      std::size_t pos = stem.find("_");
      std::string second = stem.substr(1, pos - std::size_t(1));
      int ord = std::stoi(second);
      stem = stem.substr(pos + std::size_t(1));

      pos = stem.find("_");
      second = stem.substr(1, pos - std::size_t(1));
      int group = std::stoi(second);

      second = stem.substr(pos + std::size_t(2));
      int num = std::stoi(second);

      if (ord > max_ord) {
        continue;
      }
      if (ord < min_ord) {
        continue;
      }

      // std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

      // std::cout << entry.path().filename() << std::endl;
      // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num
      // <<std::endl;

      if (ggm.size() < ord + 1) {
        ggm.resize(ord + 1);
      }
      if (ggm[ord].size() < group + 1) {
        ggm[ord].resize(group + 1);
      }
      if (ggm[ord][group].gp_vec.size() < num + 1) {
        ggm[ord][group].gp_vec.resize(num + 1);
      }

      git_pair temppair;
      pair_read(entry.path().string(), temppair);

      ggm[ord][group].gp_vec[num] = temppair;
    }
  }

  // for(int ord=2; ord< ggm.size(); ord++){
  // for(int group=0; group< ggm[ord].size(); group++){
  // for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++)
}

void TamiGraph::graph_read(std::string filename, graph_t &g) {

  std::ifstream infile_stream;

  std::vector<int> svec, tvec, statvec, spinvec;

  infile_stream.open(filename);

  if (infile_stream.fail()) // checks to see if file opened
  {
    std::cout << filename << std::endl;
    throw std::runtime_error("Could not open input file");
  }

  std::string line;
  // std::getline(infile_stream,line);

  while (std::getline(infile_stream, line)) {

    std::stringstream ss(line);
    int source, target, stat, spin;

    std::string laststring;
    bool read = bool(ss >> laststring);
    if (read) {
      source = std::stoi(laststring);
      ss >> target >> stat >> spin; // >> kx >> ky >> realW>> imagW;

      svec.push_back(source);
      tvec.push_back(target);
      statvec.push_back(stat);
      spinvec.push_back(spin);

      // std::cout<<"s t stat spin were "<< source <<" "<<target <<" "<<stat<<"
      // "<<spin<<std::endl;
    }
  }

  infile_stream.close();

  // std::cout<<"File had num lines="<<svec.size()<<std::endl;

  int val = 0;

  for (int i = 0; i < svec.size(); i++) {
    if (svec[i] > val) {
      val = svec[i];
    }
    if (tvec[i] > val) {
      val = tvec[i];
    }
  }

  int n = val + 1;
  graph_t loaded(n);

  boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
  int id = 0;
  for (boost::tie(vi, vi_end) = vertices(loaded); vi != vi_end; ++vi, ++id) {
    loaded[*vi].index_ = id;
  }

  for (int edge = 0; edge < svec.size(); edge++) {

    boost::graph_traits<graph_t>::vertex_descriptor source, target;

    boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = vertices(loaded); vi != vi_end; ++vi, ++id) {
      if (loaded[*vi].index_ == svec[edge]) {
        source = *vi;
      }
      if (loaded[*vi].index_ == tvec[edge]) {
        target = *vi;
      }
    }

    add_edge(source, target,
             edge_info(TamiBase::stat_type(statvec[edge]), -1,
                       spin_type(spinvec[edge])),
             loaded);
  }

  g = loaded;
}

void TamiGraph::pair_read(std::string filename, git_pair &p) {

  std::ifstream infile_stream;

  std::vector<int> svec, tvec, statvec, spinvec;
  std::vector<int> pairID;
  std::vector<TamiBase::epsilon_t> epsvec;
  std::vector<TamiBase::alpha_t> alphavec;
  TamiBase::epsilon_t eps;
  TamiBase::alpha_t alpha;

  // graph_t g1, g2;
  git_perm_set_t pst;
  git_perm_t temp_pst;
  // git_pair temp_pair;

  infile_stream.open(filename);

  if (infile_stream.fail()) // checks to see if file opened
  {
    std::cout << filename << std::endl;
    throw std::runtime_error("Could not open input file");
  }

  std::string line;
  // std::getline(infile_stream,line);

  while (std::getline(infile_stream, line)) {

    std::stringstream ss(line);
    int source, target, stat, spin;

    std::string number;
    std::string left("{");
    std::string right("}");

    std::string laststring;
    bool read = bool(ss >> laststring);
    if (read) {
      pairID.push_back(stoi(laststring));
      // source=std::stoi(laststring);
      // std::cout<<"Last string was "<<laststring<<std::endl;
      if (std::stoi(laststring) < 2) {
        ss >> source >> target >> stat >> spin; // >> kx >> ky >> realW>> imagW;

        svec.push_back(source);
        tvec.push_back(target);
        statvec.push_back(stat);
        spinvec.push_back(spin);

        int par = 0;

        // ss>> laststring; // this is parenthesis

        while (ss >> number) {

          if (number == left || number == right) {
            par++;
          } else {
            if (par < 2) {
              eps.push_back(std::stoi(number));
            }
            if (par > 1) {
              alpha.push_back(std::stoi(number));
            }
          }

          // std::cout<<"par is "<<par<<" entry is "<< number<<std::endl;
          // if(is_number(number)){

          // if(par<2){
          // std::cout<<"par was "<<par<<" pushing back "<<number<<std::endl;
          // eps.push_back(std::stoi(number));
          // }
          // if(par>1 && par<4){
          // std::cout<<"par was "<<par<<" pushing back "<<number<<std::endl;
          // alpha.push_back(std::stoi(number));
          // }

          // }else{par++;}
        }
        // std::cout<<"s t stat spin were "<< source <<" "<<target <<"
        // "<<stat<<" "<<spin<<std::endl;
        epsvec.push_back(eps);
        alphavec.push_back(alpha);
        eps.clear();
        alpha.clear();
      } else {

        // temp_pst.push_back(std::stoi(laststring));
        while (ss >> number) {
          // std::cout<<number<<" "<<std::endl;
          temp_pst.push_back(std::stoi(number));
        }
        pst.push_back(temp_pst);
        temp_pst.clear();
      }
    }
  }

  // at this stage only the pst is loaded
  p.pst_ = pst;
  //

  infile_stream.close();

  // std::cout<<"File had num lines="<<svec.size()<<std::endl;

  int val = 0;

  for (int i = 0; i < svec.size(); i++) {
    if (svec[i] > val) {
      val = svec[i];
    }
    if (tvec[i] > val) {
      val = tvec[i];
    }
  }

  int n = val + 1;
  graph_t loaded(n);
  graph_t g2(n);

  boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
  int id = 0;
  for (boost::tie(vi, vi_end) = vertices(loaded); vi != vi_end; ++vi, ++id) {
    loaded[*vi].index_ = id;
  }

  id = 0;
  for (boost::tie(vi, vi_end) = vertices(g2); vi != vi_end; ++vi, ++id) {
    g2[*vi].index_ = id;
  }

  for (int edge = 0; edge < svec.size(); edge++) {

    if (pairID[edge] == 0) {
      boost::graph_traits<graph_t>::vertex_descriptor source, target;

      boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = vertices(loaded); vi != vi_end;
           ++vi, ++id) {
        if (loaded[*vi].index_ == svec[edge]) {
          source = *vi;
        }
        if (loaded[*vi].index_ == tvec[edge]) {
          target = *vi;
        }
      }

      add_edge(source, target,
               edge_info(epsvec[edge], alphavec[edge],
                         TamiBase::stat_type(statvec[edge]), -1,
                         spin_type(spinvec[edge])),
               loaded);
    }
    if (pairID[edge] == 1) {
      boost::graph_traits<graph_t>::vertex_descriptor source, target;

      boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = vertices(g2); vi != vi_end; ++vi, ++id) {
        if (g2[*vi].index_ == svec[edge]) {
          source = *vi;
        }
        if (g2[*vi].index_ == tvec[edge]) {
          target = *vi;
        }
      }

      add_edge(source, target,
               edge_info(epsvec[edge], alphavec[edge],
                         TamiBase::stat_type(statvec[edge]), -1,
                         spin_type(spinvec[edge])),
               g2);
    }
  }

  p.g1_ = loaded;
  p.g2_ = g2;
}

void TamiGraph::number_vertices(graph_t &g) {

  boost::graph_traits<graph_t>::vertex_iterator vi, v_end;
  // set_all_vertex_freq(g);
  int i = 0;

  for (boost::tie(vi, v_end) = vertices(g); vi != v_end; ++vi) {
    g[*vi].index_ = i;
    i++;
  }
}

void TamiGraph::print_all_edge_info(graph_t &g) {

  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  // Looking through all edges to find source and targets : this could be
  // improved?
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

    std::cout << "Edge (" << g[source(*ei, g)].index_ << ","
              << g[target(*ei, g)].index_ << ") with bm "
              << g[*ei].band_element.first << "," << g[*ei].band_element.second
              << " of stat_type " << g[*ei].g_struct_.stat_ << " with loop id "
              << g[*ei].fermi_loop_id << " with spin " << g[*ei].spin
              << " with label " << g[*ei].label << " label=[ "; // << std::endl;

    for (int i = 0; i < g[*ei].g_struct_.alpha_.size(); i++) {

      std::cout << g[*ei].g_struct_.alpha_[i] << " ";
    }
    std::cout << "] and eps=[";

    for (int i = 0; i < g[*ei].g_struct_.eps_.size(); i++) {

      std::cout << g[*ei].g_struct_.eps_[i] << " ";
    }

    std::cout << "]" << std::endl;

    // std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<<
    // total<<std::endl;
  }
}

void TamiGraph::generate_sigma_ct(graph_t &g_in, std::vector<graph_t> &ct_vec,
                                  int maxdots) {
  // max is the maximum number of dots to add to one graph

  graph_t gc = g_in;

  edge_vector_t fermi_edges;
  // find_fermionic_edges(g_in, fermi_edges);

  if (graph_type == TamiBase::density || graph_type == TamiBase::DOS ||
      graph_type == TamiBase::Greens) {
    find_fermionic_edges(g_in, fermi_edges);
  } else {
    find_internal_fermionic_edges(g_in, fermi_edges);
  }
  // we now create all combinations of
  int n = fermi_edges.size();
  for (int i = 1; i <= maxdots; i++) {

    int r = i;
    std::vector<std::vector<int>> list;

    combinations_repetition(n, r, list);

    // std::cout<<"Comb list has size "<<list.size()<<std::endl;
    // for(int i=0; i< list.size(); i++){

    // for(int j=0; j<list[i].size(); j++){
    // std::cout<<list[i][j];
    // }
    // std::cout<<std::endl;

    // }
    // now for each list item we will make a new graph that is then converted to
    // a ct graph

    for (int comb = 0; comb < list.size(); comb++) {

      // convert comb list to an index and length pair
      std::vector<std::pair<int, int>>
          chainlist; // vector of size n that contains how MANY dots to add to
                     // each current line for this particular combination
                     // std::cout<<"On Comb "<< comb<<std::endl;
                     // for(int j=0; j<list[comb].size(); j++){
                     // std::cout<<list[comb][j];
                     // }
                     // std::cout<<std::endl;

      // making chain lists
      for (int i = 0; i < n; i++) {

        std::pair<int, int> temp(i, 0);

        for (int j = 0; j < list[comb].size(); j++) {

          if (list[comb][j] == i) {
            temp.second++;
          }
        }

        if (temp.second > 0) {
          chainlist.push_back(temp);
        }
      }

      // for(int i=0; i<chainlist.size(); i++){

      // std::cout<<"i="<<i<<" "<<chainlist[i].first<<" "<<
      // chainlist[i].second<<std::endl;

      // }
      //////

      // next for each chain
      graph_t gtemp = gc;
      for (int i = 0; i < chainlist.size(); i++) {

        // std::cout<<"i="<<i<<" "<<chainlist[i].first<<" "<<
        // chainlist[i].second<<std::endl;
        insert_chain(gc, gtemp, fermi_edges[chainlist[i].first],
                     chainlist[i].second);

        // print_all_edge_info(gc);
      }

      ct_vec.push_back(gtemp);
    }
  }
}

// on first call presumes g_in==ctg and then modifies only the ctg
void TamiGraph::insert_chain(graph_t &g_in, graph_t &ctg, edge_t &e,
                             int length) {

  // std::cout<<"adding chain to edge from "<<g_in[source(e,g_in)].index_<<" -
  // to - "<<g_in[target(e,g_in)].index_<<std::endl;

  // first find the edge and vertices in gc with the correct index values

  vertex_t vA, vB;
  edge_t cedge;

  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

  for (boost::tie(ei, ei_end) = edges(ctg); ei != ei_end; ++ei) {

    if ((ctg[source(*ei, ctg)].index_ == g_in[source(e, g_in)].index_) &&
        (ctg[target(*ei, ctg)].index_ == g_in[target(e, g_in)].index_) &&
        g_in[*ei].g_struct_.stat_ == TamiBase::Fermi) {
      vA = source(*ei, ctg);
      vB = target(*ei, ctg);
      cedge = *ei;
    }
  }

  // add vertices and give them new index_ values
  int num_vert = num_vertices(ctg);
  vertex_vector_t vv(length);
  for (int i = 0; i < length; i++) {
    vv[i] = add_vertex(ctg);
    ctg[vv[i]].index_ = num_vert + i;
  }

  // connect inner vertices
  // std::cout<<"Connecting inner"<<std::endl;
  for (int i = 1; i < length; i++) {

    add_edge(vv[i - 1], vv[i],
             edge_info(g_in[e].g_struct_.eps_, g_in[e].g_struct_.alpha_,
                       g_in[e].g_struct_.stat_, g_in[e].fermi_loop_id,
                       g_in[e].spin),
             ctg);

    ctg[boost::graph_bundle].ct_alphas.push_back(g_in[e].g_struct_.alpha_);
    int ctindex=std::distance(  std::begin(g_in[e].g_struct_.eps_), 
                                std::find_if( std::begin(g_in[e].g_struct_.eps_), 
                                std::end(g_in[e].g_struct_.eps_), [](auto x) { return x != 0; }));
    std::cout<<"Found ctindex of "<< ctindex<<std::endl;
    ctg[boost::graph_bundle].ct_alpha_index.push_back(ctindex);
  }

  // connect chain to source and target
  // std::cout<<"Connecting outer"<<std::endl;

  add_edge(vA, vv[0],
           edge_info(g_in[e].g_struct_.eps_, g_in[e].g_struct_.alpha_,
                     g_in[e].g_struct_.stat_, g_in[e].fermi_loop_id,
                     g_in[e].spin),
           ctg);

    ctg[boost::graph_bundle].ct_alphas.push_back(g_in[e].g_struct_.alpha_);
    int ctindex=std::distance(  std::begin(g_in[e].g_struct_.eps_), 
                                std::find_if( std::begin(g_in[e].g_struct_.eps_), 
                                std::end(g_in[e].g_struct_.eps_), [](auto x) { return x != 0; }));
    std::cout<<"Found ctindex of "<< ctindex<<std::endl;
    ctg[boost::graph_bundle].ct_alpha_index.push_back(ctindex);

  add_edge(vv[vv.size() - 1], vB,
           edge_info(g_in[e].g_struct_.eps_, g_in[e].g_struct_.alpha_,
                     g_in[e].g_struct_.stat_, g_in[e].fermi_loop_id,
                     g_in[e].spin),
           ctg);

  // remove original edge
  // std::cout<<"Removing original"<<std::endl;

  remove_edge(cedge, ctg);
  // std::cout<<"Edge removed?"<<std::endl;

  ctg[boost::graph_bundle].sigma_ct_count += length;

  // std::cout<<"Exiting insert chain function"<<std::endl;
}

void TamiGraph::trojan_graph_to_R0(trojan_graph &tg, TamiBase::g_prod_t &R0) {
  TamiGraph::graph_t g = tg.graph;
  this->graph_to_R0(g, R0);
}

double TamiGraph::trojan_get_prefactor(trojan_graph &tg, int order) {
  TamiGraph::graph_t g = tg.graph;
  return this->get_prefactor(g, order);
}

void TamiGraph::trojan_generate_sigma_ct(trojan_graph &tg_in,
                                         std::vector<graph_t> &ct_vec,
                                         int maxdots) {
  TamiGraph::graph_t g_in = tg_in.graph;
  this->generate_sigma_ct(g_in, ct_vec, maxdots);
}

void TamiGraph::trojan_extract_bose_alphas(trojan_graph &tg,
                                           std::vector<TamiBase::alpha_t> &bose) {
  TamiGraph::graph_t g = tg.graph;
  this->extract_bose_alphas(g, bose);
}