#include "ft.hpp"
#include "pt_ami.hpp"



std::complex<double> PtAmiBase::evaluate(AmiBase::ami_parms &parms, ft_terms &ft_terms,
                                            AmiBase::ami_vars &external){

  std::complex<double> term;
  std::complex<double> output(0,0);


for (int i = 0; i < ft_terms.size(); i++) {
  // std::cout<<"Term "<<i<<std::endl;
    term=evaluate_term(parms,ft_terms[i],external);
    output+=term;
  
}


return output;

}
 
std::complex<double> PtAmiBase::evaluate_term(AmiBase::ami_parms &parms, ft_term &ft_term,
                                            AmiBase::ami_vars &external){

std::complex<double> gprod,fprod;

  gprod = amibase.eval_gprod(parms, ft_term.g_prod_, external);
  FermiTree::vertex_t r=FT.get_root(ft_term.ft_);
  fprod = eval_ft(parms, ft_term.ft_,r, external);
  
  std::complex<double> output(0, 0);

  output = ft_term.sign_ * gprod * fprod;




return output;

}

// std::complex<double> AmiBase::eval_fprod(ami_parms &parms, pole_array_t &p_list,
                                         // ami_vars &external)

std::complex<double> PtAmiBase::eval_ft(AmiBase::ami_parms &parms, FermiTree::fermi_tree_t &ft1,  FermiTree::vertex_t &v, AmiBase::ami_vars &external){
  
 
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
      
      ft1[v].value_=amibase.fermi_pole(parms, ft1[v].pole_, external);
      
      output=ft1[v].value_*ft1[v].prefactor_;
      break;
    
  }
  
  

return output;  
}

