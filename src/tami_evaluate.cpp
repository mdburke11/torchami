#include "tami_ft.hpp"
#include "tami_base.hpp"



std::complex<double> TamiBase::evaluate(TamiBase::ami_parms &parms, ft_terms &ft_terms,
                                            TamiBase::ami_vars &external){

  std::complex<double> term;
  std::complex<double> output(0,0);


for (int i = 0; i < ft_terms.size(); i++) {
  // std::cout<<"Term "<<i<<std::endl;
    term=evaluate_term(parms,ft_terms[i],external);
    output+=term;
  
}


return output;

}
 
std::complex<double> TamiBase::evaluate_term(TamiBase::ami_parms &parms, ft_term &ft_term,
                                            TamiBase::ami_vars &external){

std::complex<double> gprod,fprod;

  gprod = eval_gprod(parms, ft_term.g_prod_, external);
  FermiTree::vertex_t r=FT.get_root(ft_term.ft_);
  fprod = eval_ft(parms, ft_term.ft_,r, external);
  
  std::complex<double> output(0, 0);

  output = ft_term.sign_ * gprod * fprod;




return output;

}

// std::complex<double> TamiBase::eval_fprod(ami_parms &parms, pole_array_t &p_list,
                                         // ami_vars &external)

std::complex<double> TamiBase::eval_ft(TamiBase::ami_parms &parms, FermiTree::fermi_tree_t &ft1,  FermiTree::vertex_t &v, TamiBase::ami_vars &external){
  
 
  std::vector<FermiTree::vertex_t> level;
  
  FT.get_next_level( ft1, v, level);

  std::complex<double> output;
  std::complex<double> aoutput(0,0);
  std::complex<double> moutput(1,0);
  
  // std::cout<<"In Evaluate the Prefactor is "<<ft1[v].prefactor_<<" on vertex "<< ft1[v].index_<< " with operation "<< ft1[v].operation_<<std::endl;
  
  switch(ft1[v].operation_){
    
    case 0:
      for(int i=0;i< level.size(); i++){
      
        aoutput+=eval_ft(parms,ft1,level[i],external)*ft1[v].prefactor_;
        
      }
      output=aoutput;
      break;
    case 1:
      for(int i=0;i< level.size(); i++){
        // std::cout<<"On evaluate prefactor is "<<ft1[v].prefactor_<<" on vertex "<< ft1[v].index_<< std::endl;
        moutput=moutput*eval_ft(parms,ft1,level[i],external); // added in prefactor to the mult term. 
        
      }
      moutput=moutput*ft1[v].prefactor_;// only one instance of the prefactor if there is a multiplication
      output=moutput;
      break;
    case 2: 
      
      ft1[v].value_=fermi_pole(parms, ft1[v].pole_, external);
      
      output=ft1[v].value_*ft1[v].prefactor_;
      break;
    
  }
  
  

return output;  
}

/**
 *
 * Numerical evaluation of a product of Green's functions. Used both in `terms`
 * and `R_t` constructions.
 * @param[in] parms : `ami_parms` object, basic parameters for AMI.
 * @param[in] g_prod : `g_prod_t` a list of `g_struct` that is interpretted as
 * \f$ \prod{G_i} \f$.
 * @param[in] external : Input external variables in a `ami_vars` struct.
 * @return Single value for product of Green's functions.
 */
std::complex<double> TamiBase::eval_gprod(ami_parms &parms, g_prod_t g_prod,
                                         ami_vars external) {
  std::complex<double> output(0, 0);

  std::complex<double> denom_prod(1, 0);
  double prefactor = external.prefactor;


  double E_REG = parms.E_REG_;

  bool verbose = false;


  for (int i = 0; i < g_prod.size(); i++) {
    std::complex<double> alphadenom(0, 0);
    std::complex<double> epsdenom(0, 0);

    for (int a = 0; a < g_prod[i].alpha_.size(); a++) {
      alphadenom += double(g_prod[i].alpha_[a]) * external.frequency_[a];
    }

    std::complex<double> zero(0, 0);
    std::complex<double> im(0, 1);

    for (int a = 0; a < g_prod[i].eps_.size(); a++) {
      epsdenom += double(g_prod[i].eps_[a]) * external.energy_[a];

    }

  
    denom_prod = denom_prod * (alphadenom + epsdenom);
   
  }

  output = 1.0 / denom_prod * prefactor;

  return output;
}

