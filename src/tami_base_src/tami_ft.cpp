#include "tami_base.hpp"

// Empty constructor
TamiBase::FermiTree::FermiTree(){
  
}


void TamiBase::FermiTree::number_vertices(fermi_tree_t &g){
	
boost::graph_traits<fermi_tree_t>::vertex_iterator vi, v_end;
int i=0;

for (boost::tie(vi,v_end) = vertices(g); vi != v_end; ++vi){
g[*vi].index_=i;
i++;
   }
	
}

void TamiBase::FermiTree::initialize_ft( fermi_tree_t &ft, operation op){
  
  if(num_vertices(ft) >0){ throw std::runtime_error("Attempted to initialize a non-empty fermi tree");}

vertex_t seed=add_vertex(vertex_info(op,0),ft); 
number_vertices(ft);
  
  
}

bool TamiBase::FermiTree::is_empty_ft(fermi_tree_t &ft){
  
  if(num_vertices(ft)==1){
    vertex_t r=get_root(ft);
    if(ft[r].operation_==end){return false;}
    
  }
  
  
  if(num_edges(ft)==0){return true;}else{return false;}
  
}


void TamiBase::FermiTree::initialize_ft( fermi_tree_t &ft){
  
if(num_vertices(ft) >0){ throw std::runtime_error("Attempted to initialize a non-empty fermi tree");}

vertex_t seed=add_vertex(vertex_info(end,0),ft); 
number_vertices(ft);
}


void TamiBase::FermiTree::initialize_ft( fermi_tree_t &ft, TamiBase::pole_struct &pole){
  
if(num_vertices(ft) >0){ throw std::runtime_error("Attempted to initialize a non-empty fermi tree");}

vertex_t seed=add_vertex(vertex_info(pole,end),ft); 
number_vertices(ft);
}

void TamiBase::FermiTree::initialize_ft( fermi_tree_t &ft, TamiBase::pole_struct &pole, double &prefactor){
  
if(num_vertices(ft) >0){ throw std::runtime_error("Attempted to initialize a non-empty fermi tree");}

vertex_t seed=add_vertex(vertex_info(pole,prefactor,end),ft); 
number_vertices(ft);
}



TamiBase::FermiTree::vertex_t TamiBase::FermiTree::get_root(fermi_tree_t  &ft){


boost::graph_traits<fermi_tree_t>::vertex_iterator vi,vie;

vertex_t vout;

for( boost::tie(vi,vie)=vertices(ft); vi !=vie; ++vi){

// if( ft[*vi].depth_==0){std::cout<<"Found root "<<std::endl;
// vout=*vi;
// }  
  if( in_degree(*vi,ft)==0){
    // std::cout<<"Found root "<<std::endl;
vout=*vi;
}  
  
  
} 
  
  return vout;
  
}


void TamiBase::FermiTree::get_roots(fermi_tree_t  &ft, std::vector<vertex_t> &vv){

vv.clear();
boost::graph_traits<fermi_tree_t>::vertex_iterator vi,vie;

vertex_t vout;

for( boost::tie(vi,vie)=vertices(ft); vi !=vie; ++vi){

// if( ft[*vi].depth_==0){std::cout<<"Found root "<<std::endl;
// vout=*vi;
// }  
  if( in_degree(*vi,ft)==0){std::cout<<"Found root "<<std::endl;
vout=*vi;
vv.push_back(vout);
}  

  
  
} 
  
  
}

void TamiBase::FermiTree::plist_to_ft(TamiBase::pole_array_t &plist, TamiBase::FermiTree::fermi_tree_t &ft){
  TamiBase::FermiTree::fermi_tree_t newft;
  
  if(plist.size()!=0){
  initialize_ft(newft,plist[0]);
  }
  for(int i=1; i<plist.size(); i++){
    
    TamiBase::FermiTree::fermi_tree_t thisft;
    initialize_ft(thisft,plist[i]);
    newft=mult_ft(newft,thisft);   
  }
  
  
  ft=newft;
  
  
}

void TamiBase::FermiTree::plist_to_ft(TamiBase::pole_array_t &plist,double sign, TamiBase::FermiTree::fermi_tree_t &ft){
  TamiBase::FermiTree::fermi_tree_t newft;
  
  if(plist.size()!=0){
  initialize_ft(newft,plist[0],sign);
  }
  for(int i=1; i<plist.size(); i++){
    
    TamiBase::FermiTree::fermi_tree_t thisft;
    initialize_ft(thisft,plist[i],sign);
    newft=mult_ft(newft,thisft);   
  }
  
  
  ft=newft;
  
  
}


void TamiBase::FermiTree::get_next_level( fermi_tree_t &ft1, vertex_t &root1, std::vector<vertex_t > &next_level){
  next_level.clear();
  
  boost::graph_traits<fermi_tree_t>::edge_iterator ei,eie;
  
  
  for(boost::tie(ei,eie)=edges(ft1); ei!=eie; ++ei){
    
    if(source(*ei,ft1)==root1){
      vertex_t ft1_target=target(*ei,ft1);
      
      
      next_level.push_back(ft1_target);
            
    }
    
    
  }
  
}


void TamiBase::FermiTree::print_graph(fermi_tree_t &ft){
  
    if(num_edges(ft)==0){
      if(num_vertices(ft)!=1){ throw std::runtime_error("multiple disconnected vertices");}
      boost::graph_traits<fermi_tree_t>::vertex_iterator vi,vie;
      for(boost::tie(vi,vie)=vertices(ft); vi!=vie; ++vi){
      
      std::cout<<"|"<<ft[*vi].index_<<" "<<ft[*vi].operation_<<" pre:"<<ft[*vi].prefactor_<<"| "<<std::endl;
      
    }
      
      
    }
  
  
    boost::graph_traits<fermi_tree_t>::edge_iterator ei,eie;
    
    for(boost::tie(ei,eie)=edges(ft); ei!=eie; ++ei){
      
      if(ft[target(*ei,ft)].operation_!=2){
      std::cout<<"|"<<ft[source(*ei,ft)].index_<<" "<<ft[source(*ei,ft)].operation_<<"| "<<" pre:"<<ft[source(*ei,ft)].prefactor_<<" ---> |"<<ft[target(*ei,ft)].index_<<" "<<ft[target(*ei,ft)].operation_<<std::endl;
      }else{
        std::cout<<"|"<<ft[source(*ei,ft)].index_<<" "<<ft[source(*ei,ft)].operation_<<"| "<<" pre:"<<ft[source(*ei,ft)].prefactor_<<" ---> |"<<ft[target(*ei,ft)].index_<<" "<<ft[target(*ei,ft)].operation_<<" pole "<< pretty_print_pole(ft[target(*ei,ft)].pole_)<<" m="<<ft[target(*ei,ft)].pole_.der_<<" pre:"<<ft[target(*ei,ft)].prefactor_<<std::endl;
      }
      
    }
    
  
}


void TamiBase::FermiTree::copy_vertex(vertex_t &root1, vertex_t &root2, fermi_tree_t &ft1, fermi_tree_t &ft2, std::pair<vertex_t, vertex_t> &map){
  
  vertex_t newv=add_vertex(vertex_info(ft1[root1].pole_,ft1[root1].value_,ft1[root1].operation_,ft1[root1].depth_,ft1[root1].prefactor_),ft2);
  add_edge(root2, newv, edge_info(),ft2);
  std::pair<vertex_t, vertex_t> this_map(root1,newv);
  map=this_map;
  
}


void TamiBase::FermiTree::copy_level(vertex_t &root1, vertex_t &root2 ,fermi_tree_t &ft1, fermi_tree_t &ft2,std::vector<std::pair<vertex_t, vertex_t>> &map_vec){
  map_vec.clear();
  std::vector<vertex_t> level;
  
  get_next_level( ft1, root1, level);
  
  // for each vertex in the level, add it to ft2 below root 2 and store the map
  
  for(int i=0; i< level.size(); i++){
    std::pair<vertex_t, vertex_t> this_map;
    copy_vertex(level[i],root2,ft1,ft2, this_map);
    map_vec.push_back(this_map);
    
  }
  
  
}

void TamiBase::FermiTree::copy_tree(fermi_tree_t &ft1, fermi_tree_t &ft2,vertex_t &root1, vertex_t &root2){
  
 std::pair<vertex_t, vertex_t> initial_map;
 
// copy root1 below root2  
 copy_vertex(root1,root2,ft1,ft2,initial_map);
 
 std::vector<std::pair<vertex_t, vertex_t>> map_vec;
 map_vec.push_back(initial_map);
 std::vector<std::pair<vertex_t, vertex_t>> this_map_vec, next_map_vec;
 
bool cont=true; 
do{
  next_map_vec.clear();

for(int i=0; i<map_vec.size(); i++){


copy_level(map_vec[i].first, map_vec[i].second, ft1,ft2, this_map_vec);

for(int m=0; m< this_map_vec.size(); m++){ next_map_vec.push_back(this_map_vec[m]);} // for some reason can't use std::copy with std::pair 


}

if(next_map_vec.size()==0){cont=false;}
map_vec=next_map_vec;

}while(cont);

  
 
  
}



//// Operations

TamiBase::FermiTree::fermi_tree_t TamiBase::FermiTree::add_ft(fermi_tree_t ft1, fermi_tree_t ft2){
  
// std::cout<<"Attempting to add ft's"<<std::endl;

number_vertices(ft1);
number_vertices(ft2);
// std::cout<<"FT1 is "<<std::endl;
// print_graph(ft1);
// std::cout<<"FT2 is "<<std::endl;
// print_graph(ft2);  


fermi_tree_t output;
initialize_ft(output);  
vertex_t nr=get_root(output);  
output[nr].operation_=0; // I guess 0 is add and 1 is mult 
 
  
// check if root of ft1 or ft2 has operation==0 if so put copy the other tree there

vertex_t root1=get_root(ft1);
vertex_t root2=get_root(ft2); 

/* if(ft1[root1].operation_==0){
  
  copy_tree(ft1,ft2,root1,root2);
  
  output=ft1;
  
  number_vertices(output);
std::cout<<"1Result of the sum is "<<std::endl;
print_graph(output);
  return output;
  
}else if(ft2[root2].operation_==0){
  copy_tree(ft2,ft1,root2,root1);
  output=ft2;
  
  number_vertices(output);
std::cout<<"2Result of the sum is "<<std::endl;
print_graph(output);
  return output;
} */
  
// std::cout<<"Prefactors are "<<ft1[root1].prefactor_<<" and "  <<ft2[root2].prefactor_<<std::endl;

copy_tree(ft1,output,root1,nr);  
copy_tree(ft2,output,root2,nr);  

// number_vertices(output);
// std::cout<<"Result of the sum is "<<std::endl;
// print_graph(output);


return output;  
}

TamiBase::FermiTree::fermi_tree_t TamiBase::FermiTree::mult_ft(fermi_tree_t ft1, fermi_tree_t ft2){
  
  // std::cout<<"Entering mult"<<std::endl;
  
fermi_tree_t output;
initialize_ft(output);  
vertex_t nr=get_root(output);  
output[nr].operation_=1; // I guess 0 is add and 1 is mult 

vertex_t root1=get_root(ft1);
vertex_t root2=get_root(ft2);  
  
copy_tree(ft1,output,root1,nr);  
copy_tree(ft2,output,root2,nr);  

// std::cout<<"Exiting mult"<<std::endl;
  
return output;  
}

void TamiBase::FermiTree::update_prefactors(fermi_tree_t &ft, double sign){
  
  
boost::graph_traits<fermi_tree_t>::vertex_iterator vi,vie;

vertex_t vout;

for( boost::tie(vi,vie)=vertices(ft); vi !=vie; ++vi){
  
  if(ft[*vi].operation_==end){ 
  ft[*vi].prefactor_=ft[*vi].prefactor_*sign;
  }
  
  
}
  
  
  
}

void TamiBase::FermiTree::mult_prefactor(fermi_tree_t &ft, double sign){
  
vertex_t root=get_root(ft);
// std::cout<<"In mult we have root with prefactor "<<ft[root].prefactor_<<std::endl;
if(ft[root].operation_==1){
ft[root].prefactor_=ft[root].prefactor_*sign;  
}
if(ft[root].operation_==0){
  std::vector<vertex_t> level;
  
  get_next_level( ft, root, level);
  
  for(int i=0; i< level.size(); i++){
    // std::cout<<"Setting prefactor to "<<ft[level[i]].prefactor_*sign<<std::endl;
    ft[level[i]].prefactor_=ft[level[i]].prefactor_*sign;
  }
  
}


  
}

std::string TamiBase::FermiTree::pretty_print(fermi_tree_t &ft){
  
  vertex_t root=get_root(ft);
  
  std::stringstream ss;
  
  ss<<"[";
  ss<<pretty_print(ft,root);
  ss<<"]";
  
  
  return ss.str();
  
}

std::string TamiBase::FermiTree::pretty_print(fermi_tree_t &ft, vertex_t &v){
  
  std::stringstream ss;
  std::vector<vertex_t> level;
  
  get_next_level( ft, v, level);
  // bool is_root= get_root(ft)==v;
  
    if(ft[v].operation_==add){
      // if(!is_root){ss<<"+";}
      
    ss<<"("<<pretty_print_level(ft,level,ft[v].operation_)<<")";
    }else{
      bool need_brackets=false;
      
      if(ft[v].prefactor_<0){ss<<"("<<ft[v].prefactor_<<")";}
      
      for(int i=0; i<level.size(); i++){
        
        if(ft[level[i]].operation_==add){need_brackets=true; break;}
        if(ft[level[0]].prefactor_<0){need_brackets=true; break;}
        // if(v==get_root(ft) && ft[get_root(ft)].operation_==mult){need_brackets=true; break;}
      }
      // need_brackets=true;
      if(need_brackets){
        ss<<"("<<pretty_print_level(ft,level,ft[v].operation_)<<")";
      }else{
      
      ss<<pretty_print_level(ft,level,ft[v].operation_);
      }
    }
  
  return ss.str();
  
}
std::string TamiBase::FermiTree::pretty_print_level(fermi_tree_t &ft, std::vector<vertex_t> &vv, int source_op){
  
  std::stringstream ss;
  
  for(int i=0; i<vv.size(); i++){
    // bool is_root= get_root(ft)==vv[i];
    
    switch(ft[vv[i]].operation_){
      
      case 0:
        // if(!is_root){ss<<"+";}
        
        ss<<pretty_print(ft,vv[i]);
        break;
      case 1:
        
        if(source_op==0 && i!=0){ ss<<"+";}
        
        ss<<pretty_print(ft,vv[i]);
        break;
      case 2:
        if(source_op==1){ss<<"(";}
        if(ft[vv[i]].prefactor_==-1){ss<<"-";}
        if(ft[vv[i]].prefactor_==1 && i!=0){ss<<"+";}
        if(ft[vv[i]].pole_.der_>2){
        ss<<"\frac{1}{"<< int(ft[vv[i]].pole_.der_-1)<<"!}";          
        }
        ss<<"f";
        for(int m=0;m<ft[vv[i]].pole_.der_; m++){
          ss<<"'";
        }
        if(!is_empty_ft(ft)){
        ss<<"("<< pretty_print_pole(ft[vv[i]].pole_)<<")";
        }
        if(source_op==1){ss<<")";}
        // output=ft1[v].value_*ft1[v].prefactor_;
        break;
      
      
    }
    
    
  }
  
  return ss.str();
}

std::string TamiBase::FermiTree::pretty_print(fermi_tree_t &ft, std::vector<vertex_t> &vv, int op){
  
  std::stringstream ss;
  
  for(int i=0; i< vv.size(); i++){
      
    switch(op){
    
    case 0:
      ss<<"+";
      ss<< pretty_print(ft,vv[i]);
      
      
      break;
    case 1:
      ss<<pretty_print(ft,vv[i]);
    case 2:
      if(ft[vv[i]].prefactor_==-1){ss<<"-";}
      if(ft[vv[i]].prefactor_==1){ss<<"+";}
      ss<<"f";
      for(int m=0;m<ft[vv[i]].pole_.der_; m++){
        ss<<"'";
      }
      if(!is_empty_ft(ft)){
      ss<<"("<< pretty_print_pole(ft[vv[i]].pole_)<<")";
      }
      // output=ft1[v].value_*ft1[v].prefactor_;
      break;
    
    
  }
    
    
    
  }
  
  return ss.str(); 
  
  
}


std::string TamiBase::FermiTree::pretty_print_ft(fermi_tree_t &ft1, vertex_t &v){
  
  
 
  std::stringstream ss;
 
  std::vector<vertex_t> level;
  
  get_next_level( ft1, v, level);
  
  bool is_root=get_root(ft1)==v;

  std::complex<double> output;
  std::complex<double> aoutput(0,0);
  std::complex<double> moutput(1,0);
  
  // if(is_root){
  
      // if(ft1[v].prefactor_==-1){ss<<"-";}
      // if(ft1[v].prefactor_==1){ss<<"+";}
  
    // if(ft1[v].prefactor_==-1){ss<<"(-1)";}
    // if(ft1[v].prefactor_==1){ss<<"(1)";}
      
  // }
  
  
  switch(ft1[v].operation_){
    
    case 0:
      for(int i=0;i< level.size(); i++){
      
        
        ss<< pretty_print_ft(ft1,level[i]);
        
        // aoutput+=eval_ft(ft1,level[i]);
        
      }
      output=aoutput;
      break;
    case 1:
      for(int i=0;i< level.size(); i++){
       if(is_root){ss<<"("<<pretty_print_ft(ft1,level[i])<<")";}
       else{
         
       ss<<")("<<pretty_print_ft(ft1,level[i]);
        // ss<<"("<<pretty_print_ft(ft1,level[i]);
       }

       // moutput=moutput*eval_ft(ft1,level[i]);
        
      }
      // output=moutput;
      break;
    case 2: 
      
      // std::cout<<"DER IS "<<ft1[v].pole_.der_<<std::endl;
      if(ft1[v].prefactor_==-1){ss<<"-";}
      if(ft1[v].prefactor_==1){ss<<"+";}
      ss<<"f";
      for(int m=0;m<ft1[v].pole_.der_; m++){
        ss<<"'";
      }
      if(!is_empty_ft(ft1)){
      ss<<"("<< pretty_print_pole(ft1[v].pole_)<<")";
      }
      // output=ft1[v].value_*ft1[v].prefactor_;
      break;
    
  }
  
  // if(level.size()==0){ ss<<")";}
  

// std::cout<<std::endl<<"This part is "<<ss.str()<<std::endl;

return ss.str();  
}

std::string TamiBase::FermiTree::pretty_print_pole(TamiBase::pole_struct &pole){
  
  std::stringstream ss;
  
  bool first=true;
  
  for(int i=0; i< pole.alpha_.size(); i++){
    if(pole.alpha_[i]!=0){
    if(first){
    if(pole.alpha_[i]==1){ ss<<"\\nu_"<<i;}
    if(pole.alpha_[i]==-1){ ss<<"-\\nu_"<<i;}
    first=false;
    }else{
      if(pole.alpha_[i]==1){ ss<<"+\\nu_"<<i;}
    if(pole.alpha_[i]==-1){ ss<<"-\\nu_"<<i;}
    }
    }
    
  }
  
  for(int i=0; i<pole.eps_.size(); i++){
    if(pole.eps_[i]!=0){
    
    if(first){
    if(pole.eps_[i]==1){ ss<<"x_"<<i;}
    if(pole.eps_[i]==-1){ ss<<"-x_"<<i;}
    first=false;
    }else{
    if(pole.eps_[i]==1){ ss<<"+x_"<<i;}
    if(pole.eps_[i]==-1){ ss<<"-x_"<<i;}
      
      
    }
    
    
    }
  }
  
  // std::cout<<std::endl<<"This pole is "<<ss.str()<<std::endl;
  
  return ss.str();  
}


void TamiBase::FermiTree::print_vertex(vertex_t &v, fermi_tree_t &ft){
  
std::cout<<" | "<<  ft[v].index_;
if( ft[v].operation_==0){ std::cout<<" + ";}
if(ft[v].operation_==1){ std::cout<<" x ";}
if(ft[v].operation_==2){std::cout<<" end ";}  
 
std::cout<<" | ";
 
}

