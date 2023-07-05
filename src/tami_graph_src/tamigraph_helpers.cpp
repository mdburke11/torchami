#include "tami_graph.hpp"

using namespace boost;

void TamiGraph::find_fermionic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding Fermionic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==TamiBase::Fermi){
vector.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}

}

void TamiGraph::combinations_repetition(int n, int r, std::vector< std::vector<int>> &list ){
	
std::vector<int> temp_chosen(r);

combinations_repetition_util( temp_chosen, list, 0, r, 0, n-1);  	
	
	
}

void TamiGraph::combinations_repetition_util(std::vector<int> &chosen, std::vector< std::vector<int>> &list, int index, int r, int start, int end){
	
if(index==r){

list.push_back(chosen);
return;
}	

for (int i=start; i<=end; i++){
	
	chosen[index]=i;
	combinations_repetition_util(chosen,list, index+1, r, i, end); 
	
}
	
}


double TamiGraph::get_prefactor(graph_t &g, int order){

// std::cout<<"In get prefactor with order "<<order<<std::endl;
// print_all_edge_info(g);

double output;
int fermi_loops=fermi_connected_components(g);//count_fermi_loops(g);
if(ami_parameters.TYPE_==TamiBase::Pi_ppud || ami_parameters.TYPE_==TamiBase::Pi_ppuu){
  fermi_loops--;
}
if(ami_parameters.TYPE_==TamiBase::Sigma || ami_parameters.TYPE_==TamiBase::Greens || ami_parameters.TYPE_==TamiBase::density || ami_parameters.TYPE_==TamiBase::DOS || ami_parameters.TYPE_==TamiBase::Bare){
  if(order!=0){
  fermi_loops--;
  }
}

// if(ami_parameters.TYPE_==TamiBase::Greens ||  && order==0){
  // fermi_loops++;
// }

int sigma_ct=g[boost::graph_bundle].sigma_ct_count;

// std::cout<<"Found fermi loops"<< fermi_loops<<std::endl;

output=pow(-1, fermi_loops+order+sigma_ct);	
	
// std::cout<<"returning "<<output<<std::endl;
	
return output;	

}

void TamiGraph::print_ggm( gg_matrix_t &ggm){

std::cout<<"Graph group size is "<<std::endl;

for(int ord=0; ord< ggm.size(); ord++){

int check=0;
int gpcheck=0;

for(int m=0; m < ggm[ord].size();m++){
for(int i=0; i< ggm[ord][m].graph_vec.size(); i++){

// std::cout<<ord<<" "<<m<<" "<<i<<std::endl;
check++;
if(ord>2){
// std::cout<<"Graph ord "<<ord<<" group "<< m <<" num "<<i<<" has n_labels="<<ggm[ord][m].labels[i].size()<<std::endl;
}
}

for(int i=0; i< ggm[ord][m].gp_vec.size(); i++){

// std::cout<<ord<<" "<<m<<" "<<i<<std::endl;
gpcheck++;
if(ord>2){
// std::cout<<"Graph ord "<<ord<<" group "<< m <<" num "<<i<<" has n_labels="<<ggm[ord][m].labels[i].size()<<std::endl;
}
}


}	
	
std::cout<<"There were "<<check<<" graphs and "<< gpcheck <<" pairs at order "<< ord<< " in "<< ggm[ord].size()<<" groups." <<std::endl;	
// std::cout<<"There were "<<check<<" graphs of order "<< ord<< " in "<< ggm[ord].size()<<" groups."<<std::endl;	
}	
	
	
	
}

void TamiGraph::find_bose_fermi_edges(graph_t &g, edge_vector_t &bose, edge_vector_t &fermi){

bose.clear();
fermi.clear();

//std::cout<< "Finding bosonic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==TamiBase::Bose){
bose.push_back(*ei);
}

if( g[*ei].g_struct_.stat_==TamiBase::Fermi){
fermi.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}


	
	
}

void TamiGraph::find_internal_vertices(graph_t &g, vertex_vector_t &vector){
	
vector.clear();

boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

for (boost::tie(vi,vi_end)=vertices(g); vi!= vi_end; ++vi){

if (degree(*vi,g) !=1){
	vector.push_back(*vi);
}

}	
	
	
}


int TamiGraph::random_int( int min, int max){

std::uniform_int_distribution<int> int_dist(min,max);

return int_dist(rand_gen);

}

int TamiGraph::graph_order(graph_t &g){
edge_vector_t edges;

find_internal_edges_stat(g, edges, TamiBase::Bose);

// edge_vector_t fedges;
// find_internal_edges_stat(g,fedges,TamiBase::Fermi);

// int eff_bose=0;
// for(int i=0; i<fedges.size(); i++){

// if(g[fedges[i]].g_struct_.species_==1){eff_bose++;}
  
// }
// std::cout<<"In graph order function we have "<<edges.size()<<" "<<eff_bose<<std::endl;

return edges.size(); //+eff_bose;
	
}


void TamiGraph::find_force_LR_vertices(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
	
in_vv.clear();
out_vv.clear();

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!= ei_end; ++ei){

if(g[*ei].g_struct_.stat_==TamiBase::Bose){


if(degree(source(*ei,g),g)==1){
	
	in_vv.push_back(target(*ei,g));
	
}

if(degree(target(*ei,g),g)==1){
	
	out_vv.push_back(source(*ei,g));
	
}





}	
	
}

	
	return;
}


void TamiGraph::fix_force_labels(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
 
// std::cout<<"Trying to find a path on graph "<<std::endl;
// print_all_edge_info(g);

edge_vector_t ev;
vertex_vector_t vv;
bool success=false;

for(int i=0; i< in_vv.size(); i++){
  for(int j=0; j<out_vv.size(); j++){
    // std::cout<<i<<" "<<j<<std::endl;
find_path_between_vertices(g, in_vv[i], out_vv[j], success, ev, vv);
if(success){break;}
  }
  if(success){break;}
}

if(success){

for(int i=0; i< ev.size(); i++){
 
g[ev[i]].g_struct_.alpha_.back()=1; 
  
}
}else{ throw std::runtime_error("Failed to find force line");}
 
  
}

/* void TamiGraph::fix_force_labels(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
  
  std::cout<<"Fixing for labels on graph "<<std::endl;
  print_all_edge_info(g);
  
  std::cout<<"In vertices are "<< g[in_vv[0]].index_<<" "<<g[in_vv[1]].index_<<std::endl;
  std::cout<<"out vertices are "<< g[out_vv[0]].index_<<" "<<g[out_vv[1]].index_<<std::endl;
	
bool keep_going=true;
bool success=false;

edge_vector_t ev;

for (int m=0; m< in_vv.size(); m++){
  keep_going=true;
  success=false;
ev.clear();
vertex_t current=in_vv[m];
vertex_t next;
edge_t e;

std::cout<<"Starting on "<<g[in_vv[m]].index_<<std::endl;

do{
	
get_adjacent_vertex_stat(g, current, next,e, TamiBase::Fermi);
ev.push_back(e);
std::cout<<"Moved to "<<g[next].index_<<std::endl;
// g[e].g_struct_.alpha_.back()=1;

for(int i=0; i< out_vv.size(); i++){
	
	if(next==out_vv[i]){
    success=true;
		keep_going=false;
		break;
	}
  
  if(next==in_vv[m]){
    success=false;
		keep_going=false;
		break;
	}
	
}

current=next;	
	
	
}while(keep_going);

if(success){break;}

} // end m loop 	

// At this point it was succesful so we found a path and can set the external lines=1

if(success){

for(int i=0; i< ev.size(); i++){
 
g[ev[i]].g_struct_.alpha_.back()=1; 
  
}
}else{ throw std::runtime_error("Failed to find force line");}	
	
	
}
 */
// void TamiGraph::fix_force_labels(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
	
// bool keep_going=true;

// vertex_t current=in_vv[0];
// vertex_t next;
// edge_t e;

// do{
	
// get_adjacent_vertex_stat(g, current, next,e, TamiBase::Fermi);
// g[e].g_struct_.alpha_.back()=1;

// for(int i=0; i< out_vv.size(); i++){
	
	// if(next==out_vv[i]){
		// keep_going=false;
		// break;
	// }
	
// }

// current=next;	
	
	
// }while(keep_going);
	
	
	
	
// }

void TamiGraph::find_path_between_vertices(graph_t &g, vertex_t &v1, vertex_t &v2, bool &success, edge_vector_t &final_ev, vertex_vector_t &final_vv){
 success=false;
 bool keep_going=true;
 
 
edge_vector_list_t next_evl; 
edge_vector_list_t evl; 
vertex_vector_list_t vvl; 
vertex_vector_list_t next_vvl; 


vertex_vector_t v1start;

v1start.push_back(v1);
vvl.push_back(v1start);
evl.resize(vvl.size());

int steps=0;
do{
  
  steps++;
next_vvl.clear();
next_evl.clear();
	
for(int vv=0; vv< vvl.size(); vv++){

vertex_t current=vvl[vv].back();
vertex_t next;
edge_t e;

vertex_vector_t newvv=vvl[vv];
edge_vector_t newev=evl[vv];

get_adjacent_vertex_stat(g, current, next,e, TamiBase::Fermi);

// if(next==v1){keep_going=false; success=false;}
// std::cout<<"Next is "<<g[next].index_<<std::endl;
bool add=true;
for(int k=0; k<newvv.size(); k++){
  if(newvv[k]==next){ add=false;}
  
}

if(add){
newvv.push_back(next);
newev.push_back(e);

next_vvl.push_back(newvv);
next_evl.push_back(newev);
}

newvv=vvl[vv];
newev=evl[vv];

get_adjacent_vertex_stat(g, current, next,e, TamiBase::Bose);

// if(next==v1){keep_going=false; success=false;}
// if(next!=newvv[0]){
newvv.push_back(next);
newev.push_back(e);

next_vvl.push_back(newvv);
next_evl.push_back(newev);
// }


}

// std::cout<<next_vvl.size()<<std::endl;

vvl=next_vvl;
evl=next_evl;  

for(int i=0; i<vvl.size(); i++){
 
if (vvl[i].back()==v2){
keep_going=false;
final_vv=vvl[i];
final_ev=evl[i];
success=true;
}  
  
}

// std::cout<<steps<<" "<<vvl.size()<<std::endl;

if (steps>100){ keep_going=false; success=false;}

  
}while(keep_going); 


// std::cout<<"Keep going is "<< keep_going<<" success is "<< success<<std::endl; 

// std::cout<<"Vertex list is "<<std::endl;

// for (int i=0; i< final_vv.size(); i++){
// std::cout<<g[final_vv[i]].index_<<std::endl;
// }  

// std::cout<<"Edge list is "<<std::endl;

// for (int i=0; i< final_ev.size(); i++){
// std::cout<<g[source(final_ev[i],g)].index_ <<"-"<<g[target(final_ev[i],g)].index_ <<std::endl;
// }

 
  
}

void TamiGraph::put_back_legs(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
int alpha_size=0;
int eps_size=0;
// std::cout<<"Legs"<<std::endl;
do{
	// std::cout<<"In do"<<std::endl;

edge_t this_edge=random_edge(g,rand_gen);
	
alpha_size=g[this_edge].g_struct_.alpha_.size();
eps_size=g[this_edge].g_struct_.eps_.size();

// std::cout<<alpha_size<<" "<<eps_size<<std::endl;
	
}while(alpha_size==0 || eps_size==0);	

// std::cout<<"Determined sizes "<< alpha_size<< " and "<<eps_size<<std::endl;	
for(int i=0; i< in_vv.size(); i++){
	
vertex_t new_vert = add_vertex(g);
edge_t this_edge=add_edge(new_vert,in_vv[i], edge_info(TamiBase::Bose),g).first;

g[this_edge].g_struct_.alpha_.resize(alpha_size,0);
g[this_edge].g_struct_.eps_.resize(eps_size,0);
g[this_edge].label=labelled;

}	


for(int i=0; i< out_vv.size(); i++){
	
vertex_t new_vert = add_vertex(g);
edge_t this_edge=add_edge(out_vv[i],new_vert, edge_info(TamiBase::Bose),g).first;

g[this_edge].g_struct_.alpha_.resize(alpha_size,0);
g[this_edge].g_struct_.eps_.resize(eps_size,0);
}	
	
number_vertices(g);	
}


int TamiGraph::fermi_connected_components(graph_t &g){

int out	;
graph_t gc=g;

vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(gc, extern_vect_list, extern_edge_list);
// std::cout<<"Deleting legs"<<std::endl;
// print_all_edge_info(gc);
delete_legs(gc, extern_vect_list, extern_edge_list);
// std::cout<<"After Deleting legs"<<std::endl;
// print_all_edge_info(gc);

number_vertices(gc);
	
cc_graph_t ccg(num_vertices(gc));

std::vector<boost::graph_traits<cc_graph_t>::vertex_descriptor> v1(num_vertices(gc));

boost::graph_traits<cc_graph_t>::vertex_iterator i, end;
  int id = 0;
  for (boost::tie(i, end) = vertices(ccg); i != end; ++i, ++id) {
        v1[id] = *i;
  }
  
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  
  for(boost::tie(ei,ei_end)=edges(gc); ei!=ei_end; ++ei){
	
  if(gc[*ei].g_struct_.stat_==TamiBase::Fermi){  
  
	int sour=gc[source(*ei,gc)].index_;
	int targ=gc[target(*ei,gc)].index_;
	
	add_edge(v1[sour], v1[targ] , ccg);	  
  }
	  
  }
  
// conceptually now have a copy of the graph. in cc_graph_t 

std::vector<int> component(boost::num_vertices (ccg));
size_t num_components = boost::connected_components (ccg, &component[0]);  
  
out=(int)num_components;
	
	
return out;	
}


// TODO: This doesn't work for TamiBase::density graphs 
void TamiGraph::ggm_label(gg_matrix_t &ggm, int min){
	

for(int ord=min; ord< ggm.size(); ord++){
	std::cout<<"Producing ggm labels at order "<<ord<< " with size "<< ggm.size()<<std::endl;
	
	for(int group=0; group<ggm[ord].size(); group++){
		// if(group%5==0){std::cout<<"On group "<< group<<std::endl;}
	for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++){
	
	// systematic_vertex_label(ggm[ord][group].graph_vec[graph]);
// print_all_edge_info(	ggm[ord][group].graph_vec[graph]);
	 std::cout<<"Labelling graph ord"<<ord<<" group "<<group<<" num "<<graph<<std::endl;
bool success=true;

if(graph_type==TamiBase::density || graph_type==TamiBase::Greens|| graph_type==TamiBase::DOS|| graph_type==TamiBase::ENERGY || graph_type==TamiBase::FORCE || graph_type==TamiBase::Pi_ppud){

// NOT SURE why we wouldn't just try the systematic case anyways. 
	
sys_label(ggm[ord][group].graph_vec[graph], success);

if(!success){ repeated_labelling(ggm[ord][group].graph_vec[graph], success);}

}else{
repeated_labelling(ggm[ord][group].graph_vec[graph], success);
}



// sys_label(ggm[ord][group].graph_vec[graph], success);
if(!success){throw std::runtime_error("Failed to label a graph");}	
		
	}
	}
}	
	
	
}

void TamiGraph::find_internal_edges_stat(graph_t &g, edge_vector_t &vector, TamiBase::stat_type requested){
	
	
	
vector.clear();

//std::cout<< "Finding Fermionic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==requested){
//if(source(*ei,g)==target(*ei,g)){vector.push_back(*ei);}
if(degree(source(*ei,g),g) != 1 && degree(target(*ei,g),g)!=1){	
vector.push_back(*ei);
}

}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}
	
}

// NOTE: Due to force diagrams, need to not associate an external bosonic line with the actual tadpole 
void TamiGraph::find_bosonic_edge_on_three_pointvert(graph_t &g, vertex_t v, edge_t &bose_edge){


boost::graph_traits<graph_t>::in_edge_iterator ei, ei_end;
boost::graph_traits<graph_t>::out_edge_iterator eo, eo_end;

for (boost::tie(ei,ei_end) = in_edges(v,g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==TamiBase::Bose){
	
	
	

bose_edge=*ei;



}
}

for (boost::tie(eo,eo_end) = out_edges(v,g); eo != eo_end; ++eo){

if( g[*eo].g_struct_.stat_==TamiBase::Bose){
	
bose_edge=*eo;
	
}
}


}

