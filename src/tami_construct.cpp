#include "tami_ft.hpp"
#include "tami_base.hpp"
/**
 * This is the primary AMI symbolic integration function.  It takes a starting
 * integrand defined by `g_prod_t` R0 and `ami_parms` object, and returns the
 * S, P and R arrays necessary for symbolic evaluation.
 * @param[in] parms : `ami_parms` object, basic parameters for AMI.
 * @param[in] R0 : `g_prod_t` integrand to be processed.
 * @param[out] R_array: Resultant `R_t`.
 * @param[out] P_array : Resultant `P_t`.
 * @param[out] S_array : Resultant `S_t`.
 */


void TamiBase::construct(int N_INT, TamiBase::g_prod_t R0, ft_terms &terms_out) {
  terms_out.clear();

  FermiTree::fermi_tree_t ft;
  FT.initialize_ft(ft,FermiTree::mult);
  ft_term start_term(1.0, ft, R0);
  ft_terms these_terms;
  these_terms.push_back(start_term);
  
  // FT.print_graph(these_terms[0].ft_);
  // these_terms[0].ft_
  // exit(0);

  ft_terms next_terms;

  for (int index = 0; index < N_INT; index++) {
    std::cout<<"On integration step "<<index<<std::endl;
    integrate_step(index, these_terms, next_terms);
    // std::cout<<"Pretty!"<<std::endl;
    std::cout<<"Next size is "<<next_terms.size()<<std::endl;
    // std::cout<<pretty_print_ft_terms(next_terms);
    these_terms = next_terms;
  }

  // put new terms in the list
  terms_out.insert(terms_out.end(), these_terms.begin(), these_terms.end());
}

// 1. 


void TamiBase::integrate_step(int index, ft_terms &in_terms, ft_terms &out_terms) {
  
  out_terms.clear();
  ft_terms collect, this_terms;
  for (int t_index = 0; t_index < in_terms.size(); t_index++) {
    this_terms.clear();
    term_integrate_step(index, in_terms[t_index], this_terms);
    
    
    collect.insert(collect.end(), this_terms.begin(), this_terms.end());
  }
  
  // additional factorization is possible at this step? 
  ft_terms factorized_terms;
  
  factorize(collect, factorized_terms);
  out_terms=factorized_terms;//collect;// factorized_terms;
  
}

void TamiBase::term_integrate_step(int index, ft_term &in_term, ft_terms &out_terms){
  
  // create a blank term that has just the g_prod from the ft_term 
  
  TamiBase::term blank_term;
  blank_term.sign=1.0;
  // Empty pole list is ok 
  // blank_term.p_list
  //
  blank_term.g_list= in_term.g_prod_;
  
  TamiBase::terms terms, newterms;
  terms.push_back(blank_term);
  
  // amibase.print_terms(terms);
  
  integrate_Mat_ind_step(index, terms, newterms);
  
  // std::cout<<"newterms has size "<< newterms.size()<<std::endl;
  // amibase.print_terms(newterms);
  terms_to_ftterms(newterms,out_terms);
  // std::cout<<"After conversion to ft_format"<<std::endl;
  // std::cout<<"Outterms size is "<< out_terms.size()<<std::endl;
  
  // std::cout<<"Terms are "<< out_terms.size()<<"------------------"<<std::endl;
  // std::cout<<pretty_print_ft_terms(out_terms)<<std::endl;
  
  // Need factorization 
  // factorize BEFORE adding the ft_ back 
  ft_terms factorized_terms;
  
  factorize(out_terms, factorized_terms);
  
  out_terms=factorized_terms;
  
  // if(FT.is_empty_ft(in_term.ft_)){std::cout<<"Input term is empty"<<std::endl;
  
  // FT.print_graph(in_term.ft_);
  // std::cout<<pretty_print_ft_term(in_term);
  // }
  if(!FT.is_empty_ft(in_term.ft_)){
  
   for(int i=0; i< out_terms.size(); i++){
    // std::cout<<"On term i "<<i<<std::endl;
    
    out_terms[i].ft_=FT.mult_ft(out_terms[i].ft_,in_term.ft_);
   
    
  }
  }
  
  
  
  
  // std::cout<<"Factorized Terms are "<< factorized_terms.size()<<"------------------"<<std::endl;
  // std::cout<<pretty_print_ft_terms(out_terms)<<std::endl;
  
  
  // std::cout<<"Exiting term integrate step "<<std::endl;
}

void TamiBase::factorize(ft_terms &in_terms, ft_terms &out_terms){
  
  out_terms.clear();
  
  std::vector<int> used(in_terms.size(),0);
  
  // each term 
  for (int i=0; i< in_terms.size(); i++){
    // compare to other terms
    // skip terms already used 
    if(used[i]==1){continue;}
    
    ft_term this_ft;
    this_ft=in_terms[i];
    
    for(int j=i+1; j< in_terms.size(); j++){
    
      // std::cout<<"Comparing ft terms "<<i<<" "<<j<<std::endl;
      int junk=0;
      
      if( g_prod_equiv(in_terms[i].g_prod_,in_terms[j].g_prod_,junk)){
        // std::cout<<"Terms equiv "<< i<<" and "<<j<<" with sign="<<junk<<std::endl;
        
  
        // in_terms[j].sign_=in_terms[j].sign_*junk;
        // std::cout<<"Before mult in factorize"<<std::endl;
        // FT.number_vertices(in_terms[j].ft_);
        // FT.number_vertices(this_ft.ft_);
        // FT.print_graph(in_terms[j].ft_);
        FT.mult_prefactor(in_terms[j].ft_, junk);
        // FT.print_graph(in_terms[j].ft_);
        
        // std::cout<<FT.pretty_print(in_terms[j].ft_)<<std::endl;
        
        // exit(0);
        // std::cout<<FT.pretty_print(this_ft.ft_)<<std::endl;
        
        this_ft.ft_=FT.add_ft(in_terms[j].ft_,this_ft.ft_);
        // std::cout<<FT.pretty_print(this_ft.ft_)<<std::endl;
        
        used[j]=1;
        
        // exit(0);
        
      }
    
    }
    
    
    out_terms.push_back(this_ft);
    
  }
  
  
  
  
  
}


// check if two g_prods contain the same terms within overall minus sign 
bool TamiBase::g_prod_equiv(TamiBase::g_prod_t &gp1, TamiBase::g_prod_t &gp2, int &sign){
  if(gp1.size()!= gp2.size()){ return false;}
  sign=1;
  // for each G in gp1 find one in gp2 
  std::vector<int> used1(gp1.size(),0);
  std::vector<int> used2(gp2.size(),0);
  
  // don't need to track 
  for(int i=0; i< gp1.size(); i++){
    bool found=false;
    for(int j=0; j<gp2.size(); j++){
  
      // std::cout<<"Comparing G "<<i<<" "<<j<<std::endl; 
      // skip any in gp2 that were already used 
      if(used2[j]==1){continue;}
        int junk=0;
      if(g_struct_equiv(gp1[i],gp2[j],junk)){
        // std::cout<<"Equal"<<std::endl;
        sign=sign*junk;
        used2[j]=1;
        found=true;
        break;
      }
    
    }
    
    // if found is not true it means that the i'th term of gp1 is not in gp2 
    if(!found){return false;}
  }
  
  // if it clears all the terms then it is true these are equivalent within an overall prefactor of sign 
 
  return true;
  
}

std::string TamiBase::pretty_print_ft_terms(ft_terms &fts){
  
  
  std::stringstream ss;
  for(int i=0; i<fts.size(); i++){
    ss<<"Term["<<i<<"]: ";
  ss<<pretty_print_ft_term(fts[i])<<std::endl;;
  }
  
  
  
return ss.str();  
  
}

std::string TamiBase::pretty_print_ft_term(ft_term &ft){
  
  std::stringstream ss;
  FermiTree::vertex_t r=FT.get_root(ft.ft_);
  ss<<"["<<ft.sign_<<"]"<<FT.pretty_print(ft.ft_)<<pretty_print_gprod(ft.g_prod_);
  
  
return ss.str();  
  
}

std::string TamiBase::pretty_print_gprod(TamiBase::g_prod_t &gp){
  
  std::stringstream ss;
  
  for(int i=0; i<gp.size(); i++){
  ss<<"\\frac{1}{"<<pretty_print_g(gp[i])<<"}";
  
  }
  
  return ss.str(); 
  
}

std::string TamiBase::pretty_print_g(TamiBase::g_struct &g){
  
  std::stringstream ss;
  bool first=true;
  
  for(int i=0; i< g.alpha_.size(); i++){
    if(g.alpha_[i]!=0){
      if(first){
    
    if(g.alpha_[i]==1){ ss<<"\\nu_"<<i;}
    if(g.alpha_[i]==-1){ ss<<"-\\nu_"<<i;}
      first=false;
      }else{
    if(g.alpha_[i]==1){ ss<<"+\\nu_"<<i;}
    if(g.alpha_[i]==-1){ ss<<"-\\nu_"<<i;}
        
        
      }
    }
  
  
  }
  
  for(int i=0; i<g.eps_.size(); i++){
    
    if(g.eps_[i]!=0){
    
    if(first){
    if(g.eps_[i]==1){ ss<<"x_"<<i;}
    if(g.eps_[i]==-1){ ss<<"-x_"<<i;}
    first=false;
    }else{
    if(g.eps_[i]==1){ ss<<"+x_"<<i;}
    if(g.eps_[i]==-1){ ss<<"-x_"<<i;}
      
    }   
    
    
    }
    
    
    
  }
  
  
  
  return ss.str();  
}


void TamiBase::terms_to_ftterms(TamiBase::terms &in_terms, ft_terms &out_terms){
  
  // look at each term and find unique G's and references and their prefactors.  This is already done by the terms factorization functions.
  // Going to do it again to get my head around what is happening. 
  out_terms.clear();

// if only one term then nothing to factor   
  if(in_terms.size()==1){
    // std::cout<<"Size was 1 with sign "<< in_terms[0].sign<<std::endl;
    
    out_terms.resize(1);
    out_terms[0].g_prod_=in_terms[0].g_list;
    
    FT.plist_to_ft(in_terms[0].p_list,in_terms[0].sign,out_terms[0].ft_);
    out_terms[0].sign_=in_terms[0].sign;
    // std::cout<<"Exit1"<<std::endl;
    return ;
  }
  
  // if size is more than one then we need to identify common factors - this was already done with the factorize_Rn function 
  
  TamiBase::Ri_t Ri;

  convert_terms_to_ri(in_terms, Ri); // this should take the terms and pull out just the producs of G's. 

  // std::cout<<"Ri size is"<< Ri.size()<<std::endl;

  TamiBase::g_prod_t unique_g;
  TamiBase::R_ref_t Rref;
  TamiBase::ref_eval_t Eval_list;
  factorize_Rn(Ri, unique_g, Rref, Eval_list);
  
  for(int i=0; i< Eval_list.size(); i++){
    ft_term newft;
    // FT.initialize_ft(newft.ft_,FermiTree::add);
    for(int j=0; j< Eval_list[i].size(); j++){
    
    // std::cout<<"i,j "<< i<<","<<j<<std::endl;
    // std::cout<<"Prefactor mult is "<<in_terms[Eval_list[i][j].first].sign*double(Eval_list[i][j].second)<<std::endl;
    
    FermiTree::fermi_tree_t ft;
    
    FT.plist_to_ft(in_terms[Eval_list[i][j].first].p_list,in_terms[Eval_list[i][j].first].sign*double(Eval_list[i][j].second), ft);
    // FT.print_graph(ft);
    if(j==0){
      newft.ft_=ft;
    }else{
    newft.ft_=FT.add_ft(newft.ft_,ft);
    }
    // newft.g_prod_.push_back(unique_g[Eval_list[i][j].first]);
        
    }
    
    TamiBase::ref_v_t pair_vec = Rref[i];
    for (int j = 0; j < pair_vec.size(); j++) {
      
      newft.g_prod_.push_back(unique_g[pair_vec[j].first]);
      
    }
    
       
    out_terms.push_back(newft);
  }


/**

* Simple function to scan through `g_prod_t` and collect an std::vector of
poles with respect to index.

*/
TamiBase::pole_array_t TamiBase::find_poles(int index, TamiBase::g_prod_t &R) {
  pole_array_t pole_array;
  pole_struct pole;

  pole.alpha_.reserve(R[0].alpha_.size());
  pole.eps_.reserve(R[0].eps_.size());

  pole.index_ = index;

  if (pole.index_ >= R[0].alpha_.size()) {
    std::cerr << "WARNING: Pole index exceeds length of g_prod_t: This may "
                 "result in errors.  Exiting find_poles for safety."
              << std::endl;
    return pole_array;
  }

  for (int i = 0; i < R.size(); i++) {
    if (R[i].alpha_[index] != 0) {
      pole.which_g_.push_back(i);

      for (int j = 0; j < R[i].alpha_.size(); j++) {
        if (j != index) {
          pole.alpha_.push_back(-R[i].alpha_[index] * R[i].alpha_[j]);
        } else {
          pole.alpha_.push_back(0);
        }
      }

      for (int j = 0; j < R[i].eps_.size(); j++) {
        pole.eps_.push_back(-R[i].alpha_[index] * R[i].eps_[j]);
      }

      // Context:  So now you have a pole... need to check if it is already
      // in the array, and then if it isn't add it.  if it is, modify the
      // multiplicity of the pole in the array it is equal to. and add it's
      // which_g_ index to the pole_array with multiplicity

      bool duplicate = false;
      for (int ploop = 0; ploop < pole_array.size(); ploop++) {
        if (pole_equiv(
                pole_array[ploop],
                pole)) { // std::cerr<<"Duplicate pole detected!"<<std::endl;
          duplicate = true;

          pole_array[ploop].multiplicity_ += 1;
          pole_array[ploop].which_g_.push_back(pole.which_g_[0]);

          break;
        }
      }

      
      // Context: extra check if the pole is only a fermionic frequency

      bool non_zero = false;
      if (drop_matsubara_poles) {
        // std::cout<<"In drop mat poles function"<<std::endl;
        int zcount = std::count(pole.eps_.begin(), pole.eps_.end(), 0);
        
        if (zcount < pole.eps_.size()) {
          non_zero = true;
        }
        
        // At this point if the energy is all zeros then if the nucount is even, then it should be ok
        if(!non_zero){
          int nucount=0;
          for(int i=0; i<pole.alpha_.size(); i++){
            if(pole.alpha_[i]!=0){ nucount++;}
          }
          // int nucount = std::count(pole.alpha_.begin(), pole.alpha_.end(), 0);
          if(nucount%2==0){ non_zero=true;}
        }
        
      } else {
        non_zero = true;
      }

      if (non_zero) {
        // if (!duplicate && non_zero ){
        if (!duplicate) {
          pole_array.push_back(pole);
        }
      }

      // }

      pole.alpha_.clear();
      pole.eps_.clear();
      pole.which_g_.clear();
    }
  }

  return pole_array;
}


/**
 * Used in pole multiplicity == 1 - aka, simple pole.  Performs residue theorem
 *on G_in and produces a new `g_prod_t` as output.
 *
 *@param[in] G_in: `g_prod_t` to evaluate residues with respect to pole.
 *@param[in] pole: `pole_struct` defining the pole.
 *
 *@return is of type `g_prod_t`.
 */
TamiBase::g_prod_t TamiBase::simple_residue(TamiBase::g_prod_t G_in,
                                          TamiBase::pole_struct pole) {
  g_prod_t residue;

  if (pole.multiplicity_ != 1) {
    std::cerr << "Simple residue called for pole with M!=1" << std::endl;
  }

  for (int i = 0; i < G_in.size(); i++) {
    if (i != pole.which_g_[0]) {
      residue.push_back(update_G_pole(G_in[i], pole));
    }
  }

  return residue;
}

double TamiBase::get_simple_sign(int index, g_prod_t &R, pole_struct pole) {
  double sign;

  sign = R[pole.which_g_[0]].alpha_[index];

  return sign;
}

/**
 * Manipulates a Green's function to replace residue variable 'z' with complex
 * pole.
 *
 */
TamiBase::g_struct TamiBase::update_G_pole(TamiBase::g_struct g_in,
                                         TamiBase::pole_struct pole) {
  TamiBase::g_struct g_new;

  for (int i = 0; i < g_in.alpha_.size(); i++) {
    if (i != pole.index_) {
      g_new.alpha_.push_back(g_in.alpha_[i] +
                             g_in.alpha_[pole.index_] * pole.alpha_[i]);
    } else {
      g_new.alpha_.push_back(0);
    }
  }

  for (int i = 0; i < g_in.eps_.size(); i++) {
    g_new.eps_.push_back(g_in.eps_[i] +
                         g_in.alpha_[pole.index_] * pole.eps_[i]);
  }

  if (g_new.alpha_.size() != g_in.alpha_.size()) {
    std::cerr << "Error: Something may be wrong. Alphas of G and pole not the "
                 "same size"
              << std::endl;
  }

  return g_new;
}

/**
 * Obtains the elements of `S_t` arrays, modified for poles with
 * multiplicity!=1.
 */
double TamiBase::get_starting_sign(TamiBase::g_prod_t G_in,
                                  TamiBase::pole_struct pole) {
  double result = 1.0;

  for (int i = 0; i < pole.which_g_.size(); i++) {
    result = result * (double)G_in[pole.which_g_[i]].alpha_[pole.index_];
  }

  result = result / (double)factorial(pole.multiplicity_ - 1);

  return result;
}

int TamiBase::factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

/**
 * Removes green's functions from `g_prod_t` to account for numerator of
 * residue (z-z0)^m
 */
TamiBase::g_prod_t TamiBase::reduce_gprod(TamiBase::g_prod_t G_in,
                                        TamiBase::pole_struct pole) {
  g_prod_t reduced;

  for (int i = 0; i < G_in.size(); i++) {
    bool add = true;

    for (int j = 0; j < pole.which_g_.size(); j++) {
      if (pole.which_g_[j] == i) {
        add = false;
        break;
      }
    }

    if (add) {
      reduced.push_back(G_in[i]);
    }
  }

  return reduced;
}
  
  // std::cout<<"Out terms is size "<< out_terms.size()<<std::endl;
  
  
  // std::cout<<"Unique_g size is "<< unique_g.size()<<std::endl;
  // std::cout<<"Rref size is "<< Rref.size()<<std::endl;
  // std::cout<<"Eval_list size is "<< Eval_list.size()<<std::endl;
  
  // for(int i=0; i< Eval_list.size(); i++){
    // std::cout<<"On i "<<i<<" eval_list size is "<< Eval_list[i].size()<<std::endl;
    // for(int m=0; m< Eval_list[i].size(); m++){
    // std::cout<<Eval_list[i][m].first<<" "<<Eval_list[i][m].second<<std::endl;
    // }
  // }
  
  // std::cout<<"Rref"<<std::endl;
  // for(int i=0; i< Rref.size(); i++){
    // std::cout<<"On i "<<i<<" rref size is "<< Rref[i].size()<<std::endl;
    // for(int m=0; m< Rref[i].size(); m++){
    // std::cout<<Rref[i][m].first<<" "<<Rref[i][m].second<<std::endl;
    // }
  // }
  
  // std::cout<<"Exiting "<<std::endl;
  // exit(0);
}

TamiBase::g_prod_t construct_example2(){

TamiBase::g_prod_t g;

// Setting up G array
// defining alpha's /// std::vector<int>


TamiBase::alpha_t alpha_1={1,0,0};
TamiBase::alpha_t alpha_2={0,1,0};
TamiBase::alpha_t alpha_3={-1,1,1};

//defining epsilon's
TamiBase::epsilon_t epsilon_1={1,0,0};
TamiBase::epsilon_t epsilon_2={0,1,0};
TamiBase::epsilon_t epsilon_3={0,0,1};

TamiBase::g_struct g1(epsilon_1,alpha_1);
TamiBase::g_struct g2(epsilon_2,alpha_2);
TamiBase::g_struct g3(epsilon_3,alpha_3);

// TamiBase::g_prod_t R0={g1,g2,g3};

// OR
TamiBase::g_prod_t R0;
R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);




return R0;

}

TamiBase::ami_vars construct_ext_example2(){


TamiBase::energy_t energy={-4,0.1,-1};

TamiBase::frequency_t frequency;

for(int i=0;i<2;i++){ frequency.push_back(std::complex<double>(0,0));}

frequency.push_back(std::complex<double>(0,M_PI/5));

double BETA=5.0;
TamiBase::ami_vars external(energy, frequency,BETA);

return external;

}





TamiBase::ami_vars construct_ext_example_J(){

TamiBase::energy_t energy={-4.64,1.02,1.04,1.05,1.06,1.07,1.08,1.09,1.11,1.23,-4.43,-4.52,1.5,1.6,1.7,1.8,1.9};
		//{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//energy_t energy={-4,1,-1,1,1,-4,1,1,1,1,1,1,1,1,1,1,1};

TamiBase::frequency_t frequency= {std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

//for(int i=0;i<9;i++){ frequency.push_back(std::complex<double>(0,0));}//frequency.push_back(std::complex<double>(0,0));}

//frequency.push_back(std::complex<double>(0, M_PI));//(0,M_PI));

double BETA=1.0;
TamiBase::ami_vars external(energy, frequency,BETA);

return external;

}


TamiBase::g_prod_t construct_example_J(){

TamiBase::g_prod_t g;

// Setting up G array
// defining alpha's


TamiBase::alpha_t alpha_1={1,0,0,0,0,0,0,0,0,0};
TamiBase::alpha_t alpha_2={0,1,0,0,0,0,0,0,0,0};
TamiBase::alpha_t alpha_3={1,1,0,0,0,0,0,0,0,-1};
TamiBase::alpha_t alpha_4={0,0,1,0,0,0,0,0,0,0};
TamiBase::alpha_t alpha_5={0,0,0,1,0,0,0,0,0,0};
TamiBase::alpha_t alpha_6={1,1,1,-1,0,0,0,0,0,-1};
TamiBase::alpha_t alpha_7={0,0,0,0,1,0,0,0,0,0};
TamiBase::alpha_t alpha_8={0,0,0,0,0,0,0,1,0,0};
TamiBase::alpha_t alpha_9={0,-1,0,0,1,0,0,1,0,0};
TamiBase::alpha_t alpha_10={0,0,1,-1,1,0,0,0,0,0};
TamiBase::alpha_t alpha_11={0,0,0,0,0,1,0,0,0,0};
TamiBase::alpha_t alpha_12={0,0,0,0,0,0,0,0,1,0};
TamiBase::alpha_t alpha_13={0,0,0,0,0,0,1,0,0,0};
TamiBase::alpha_t alpha_14={0,0,-1,1,-1,1,1,0,0,0};
TamiBase::alpha_t alpha_15={0,-1,-1,1,0,0,0,1,0,1};
TamiBase::alpha_t alpha_16={0,0,1,-1,1,-1,0,0,1,0};
TamiBase::alpha_t alpha_17={0,-1,0,0,1,-1,0,1,1,0};


//defining epsilon's
TamiBase::epsilon_t epsilon_1= {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_2= {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_3= {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_4= {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_5= {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_6= {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_7= {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_8= {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_9= {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_10={0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_11={0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_12={0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0};
TamiBase::epsilon_t epsilon_13={0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
TamiBase::epsilon_t epsilon_14={0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0};
TamiBase::epsilon_t epsilon_15={0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0};
TamiBase::epsilon_t epsilon_16={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0};
TamiBase::epsilon_t epsilon_17={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};


TamiBase::g_struct g1(epsilon_1,alpha_1);
TamiBase::g_struct g2(epsilon_2,alpha_2);
TamiBase::g_struct g3(epsilon_3,alpha_3);
TamiBase::g_struct g4(epsilon_4,alpha_4);
TamiBase::g_struct g5(epsilon_5,alpha_5);
TamiBase::g_struct g6(epsilon_6,alpha_6);
TamiBase::g_struct g7(epsilon_7,alpha_7);
TamiBase::g_struct g8(epsilon_8,alpha_8);
TamiBase::g_struct g9(epsilon_9,alpha_9);
TamiBase::g_struct g10(epsilon_10,alpha_10);
TamiBase::g_struct g11(epsilon_11,alpha_11);
TamiBase::g_struct g12(epsilon_12,alpha_12);
TamiBase::g_struct g13(epsilon_13,alpha_13);
TamiBase::g_struct g14(epsilon_14,alpha_14);
TamiBase::g_struct g15(epsilon_15,alpha_15);
TamiBase::g_struct g16(epsilon_16,alpha_16);
TamiBase::g_struct g17(epsilon_17,alpha_17);




TamiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);
R0.push_back(g4);
R0.push_back(g5);
R0.push_back(g6);
R0.push_back(g7);
R0.push_back(g8);
R0.push_back(g9);
R0.push_back(g10);
R0.push_back(g11);
R0.push_back(g12);
R0.push_back(g13);
R0.push_back(g14);
R0.push_back(g15);
R0.push_back(g16);
R0.push_back(g17);


return R0;

}




// So this has some issues once you hit the multipole problem. 
// Lets fully use the terms construction functions but then strip away the ft part on each integration step. Then translate the terms into the ft terms 
// question: Should the poles be collected in pole array or in ft?
/* void TamiBase::integrate_step(int index, ft_terms &in_terms, ft_terms &out_terms) {
  out_terms.clear();

  for (int t_index = 0; t_index < in_terms.size(); t_index++) {
    TamiBase::pole_array_t poles;
    poles = amibase.find_poles(index, in_terms[t_index].g_prod_);

    for (int i = 0; i < poles.size(); i++) {
      if (poles[i].multiplicity_ == 1) {
        ft_term new_term;
        new_term.g_prod_ = amibase.simple_residue(in_terms[t_index].g_prod_, poles[i]);
        // take sign from original term and multiply by new one
        new_term.sign_ =
            in_terms[t_index].sign_ *
            amibase.get_simple_sign(index, in_terms[t_index].g_prod_, poles[i]);
        // take poles from originating term
        
        
        FermiTree::fermi_tree_t existing= in_terms[t_index].ft_;
        FermiTree::fermi_tree_t newtree;
        FT.initialize_ft(newtree,poles[i]);
        
        new_term.ft_=FT.mult_ft(existing, newtree);
        
        
        out_terms.push_back(new_term);
      } else {
        TamiBase::terms new_terms;
        terms_general_residue(in_terms[t_index], poles[i], new_terms);

      for(int j=0; j< new_terms.size(); j++){
        FermiTree::fermi_tree_t blank=in_terms[t_index].ft_;
        FermiTree::fermi_tree_t nt;
        FT.initialize_ft(nt, new_terms
        
        
      }



        // put new terms in the list
        out_terms.insert(out_terms.end(), new_terms.begin(), new_terms.end());
      }
    }
  }
}
 */




TamiBase::ami_vars construct_4ord_ext_multipole_example(){



TamiBase::energy_t energy={1,1.1,1.2,1.31,1.4,0.01, 0.1}; //{1,1.1,1.2,1.3,1.4,0.01, 0.1};

TamiBase::frequency_t frequency= {std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

double BETA=1.0;
TamiBase::ami_vars external(energy, frequency,BETA);

return external;

}


TamiBase::g_prod_t construct_multipole_example(){

TamiBase::g_prod_t g;


// Setting up G array
// defining alpha's
// TamiBase::alpha_t alpha_1={1,0,0,1,-1};
// TamiBase::alpha_t alpha_2={0,0,0,1,0};
// TamiBase::alpha_t alpha_3={1,0,0,0,0};
// TamiBase::alpha_t alpha_4={1,0,0,0,0};
// TamiBase::alpha_t alpha_5={0,1,0,0,0};
// TamiBase::alpha_t alpha_6={-1,1,1,0,0};
// TamiBase::alpha_t alpha_7={0,0,1,0,0};

TamiBase::alpha_t alpha_1={0,0,1,1,-1};
TamiBase::alpha_t alpha_2={0,0,0,1,0};
TamiBase::alpha_t alpha_3={0,0,1,0,0};
TamiBase::alpha_t alpha_4={0,0,1,0,0};
TamiBase::alpha_t alpha_5={0,1,0,0,0};
TamiBase::alpha_t alpha_6={1,1,-1,0,0};
TamiBase::alpha_t alpha_7={1,0,0,0,0};


//defining epsilon's
TamiBase::epsilon_t epsilon_1={0,0,0,0,0,1,0};
TamiBase::epsilon_t epsilon_2={0,0,0,0,1,0,0};
TamiBase::epsilon_t epsilon_3={1,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_4={1,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_5={0,1,0,0,0,0,0};
TamiBase::epsilon_t epsilon_6={0,0,0,1,0,0,0};
TamiBase::epsilon_t epsilon_7={0,0,1,0,0,0,0};


TamiBase::g_struct g1(epsilon_1,alpha_1);
TamiBase::g_struct g2(epsilon_2,alpha_2);
TamiBase::g_struct g3(epsilon_3,alpha_3);
TamiBase::g_struct g4(epsilon_4,alpha_4);
TamiBase::g_struct g5(epsilon_5,alpha_5);
TamiBase::g_struct g6(epsilon_6,alpha_6);
TamiBase::g_struct g7(epsilon_7,alpha_7);


TamiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);
R0.push_back(g4);
R0.push_back(g5);
R0.push_back(g6);
R0.push_back(g7);



return R0;
}



TamiBase::g_prod_t construct_example1_bose(){

TamiBase::g_prod_t g;

// Setting up G array
// defining alpha's /// std::vector<int>


TamiBase::alpha_t alpha_1={1,0,0};
TamiBase::alpha_t alpha_2={1,1,0};

//defining epsilon's
TamiBase::epsilon_t epsilon_1={1,0,0};
TamiBase::epsilon_t epsilon_2={0,1,0};

TamiBase::g_struct g1(epsilon_1,alpha_1);
TamiBase::g_struct g2(epsilon_2,alpha_2);


TamiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);


return R0;

}

TamiBase::ami_vars construct_ext_example1_bose(){


TamiBase::energy_t energy={-4,0.1};

TamiBase::frequency_t frequency;

for(int i=0;i<1;i++){ frequency.push_back(std::complex<double>(0,0));}

frequency.push_back(std::complex<double>(0,M_PI));// This frequency is expected to be a bosonic matsubara or real frequency. There is no catch if this is untrue. 

double BETA=1.0;
TamiBase::ami_vars external(energy, frequency,BETA);



return external;

}


TamiBase::ami_vars construct_ext_example_Y(){



TamiBase::energy_t energy={4,-1,-1,-1,-2,2,1};

TamiBase::frequency_t frequency= {std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

double BETA=1.0;
TamiBase::ami_vars external(energy, frequency,BETA);

return external;

}


TamiBase::g_prod_t construct_example_Y(){

TamiBase::g_prod_t g;


// Setting up G array
// defining alpha's

TamiBase::alpha_t alpha_1={1,0,0,0,0};
TamiBase::alpha_t alpha_2={0,1,0,0,0};
TamiBase::alpha_t alpha_3={0,0,1,0,0};
TamiBase::alpha_t alpha_4={0,0,0,1,0};
TamiBase::alpha_t alpha_5={-1,1,0,0,1};
TamiBase::alpha_t alpha_6={0,1,1,-1,0};
TamiBase::alpha_t alpha_7={1,0,1,0,-1};

//defining epsilon's
TamiBase::epsilon_t epsilon_1={1,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_2={0,1,0,0,0,0,0};
TamiBase::epsilon_t epsilon_3={0,0,1,0,0,0,0};
TamiBase::epsilon_t epsilon_4={0,0,0,1,0,0,0};
TamiBase::epsilon_t epsilon_5={0,0,0,0,1,0,0};
TamiBase::epsilon_t epsilon_6={0,0,0,0,0,1,0};
TamiBase::epsilon_t epsilon_7={0,0,0,0,0,0,1};


TamiBase::g_struct g1(epsilon_1,alpha_1);
TamiBase::g_struct g2(epsilon_2,alpha_2);
TamiBase::g_struct g3(epsilon_3,alpha_3);
TamiBase::g_struct g4(epsilon_4,alpha_4);
TamiBase::g_struct g5(epsilon_5,alpha_5);
TamiBase::g_struct g6(epsilon_6,alpha_6);
TamiBase::g_struct g7(epsilon_7,alpha_7);


TamiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);
R0.push_back(g4);
R0.push_back(g5);
R0.push_back(g6);
R0.push_back(g7);





return R0;

}
