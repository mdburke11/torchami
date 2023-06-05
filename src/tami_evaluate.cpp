#include "tami_base.hpp"



at::Tensor TamiBase::evaluate(TamiBase::ami_parms &parms, ft_terms &ft_terms,
                                            TamiBase::ami_vars &external){

  int batch_size = external.energy_.size(0);                           

  at::Tensor term;
  at::Tensor output = at::zeros(batch_size, at::kComplexDouble).to(device);


  for (int i = 0; i < ft_terms.size(); i++) {
    // std::cout<<"Term "<<i<<std::endl;
      term=evaluate_term(parms,ft_terms[i],external).to(device);
      output+=term;
  
}

return output;

}
 
at::Tensor TamiBase::evaluate_term(TamiBase::ami_parms &parms, ft_term &ft_term,
                                            TamiBase::ami_vars &external){

  at::Tensor gprod,fprod;

  gprod = eval_gprod(parms, ft_term.g_prod_, external).to(device);
  TamiBase::FermiTree::vertex_t r=FT.get_root(ft_term.ft_);
  fprod = eval_ft(parms, ft_term.ft_,r, external).to(device);
  
  at::Tensor output = ft_term.sign_ * torch::multiply(gprod, fprod);

return output;

}

// TamiBase::complex_double TamiBase::eval_fprod(ami_parms &parms, pole_array_t &p_list,
                                         // ami_vars &external)

at::Tensor TamiBase::eval_ft(TamiBase::ami_parms &parms, TamiBase::FermiTree::fermi_tree_t &ft1,  TamiBase::FermiTree::vertex_t &v, TamiBase::ami_vars &external){

// TODO: memory optimization here? do we copy or move (aoutput or moutput) into output? and is it still on the same device?

  std::vector<TamiBase::FermiTree::vertex_t> level;
  
  FT.get_next_level( ft1, v, level);

  int batch_size = external.energy_.size(0);
  at::Tensor output = at::zeros(batch_size, at::kComplexDouble).to(device); // SEE TODO above maybe a better way to write this?
  at::Tensor aoutput = at::zeros(batch_size, at::kComplexDouble).to(device);
  at::Tensor moutput = at::ones(batch_size, at::kComplexDouble).to(device);
  
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
 * @param[in] parms : `ami_parms` object, basic parameters for TAMI.
 * @param[in] g_prod : `g_prod_t` a list of `g_struct` that is interpretted as
 * \f$ \prod{G_i} \f$.
 * @param[in] external : Input external variables in a `ami_vars` struct.
 * @return Value for product of Green's functions for all energies in external
 */
at::Tensor TamiBase::eval_gprod(ami_parms &parms, g_prod_t g_prod,
                                         ami_vars external) {
  
  int batch_size = external.energy_.size(0);
  
  at::Tensor output = at::zeros(batch_size, at::kComplexDouble).to(device); // get len(energies) answers

  at::Tensor denom_prod = at::ones(batch_size, at::kComplexDouble).to(device);
  double prefactor = external.prefactor;

  double E_REG = parms.E_REG_;
  bool verbose = false;

  for (int i = 0; i < g_prod.size(); i++) {
    TamiBase::complex_double alphadenom(0, 0); 
    at::Tensor epsdenom = at::zeros(batch_size, at::kComplexDouble).to(device);

    for (int a = 0; a < g_prod[i].alpha_.size(); a++) {
      alphadenom += double(g_prod[i].alpha_[a]) * external.frequency_[a];
    }

    for (int a = 0; a < g_prod[i].eps_.size(); a++) {
      epsdenom += double(g_prod[i].eps_[a]) * external.energy_.index({torch::indexing::Slice(),a});
    }

    //TODO: What device is this scalar on ... or does it matter? -- Hopefully this works
    //c10::Scalar alpha_denom = c10::complex<double>(alphadenom.real(), alphadenom.imag()); // cast to a tensor friendly complex value
  
    denom_prod = torch::multiply(denom_prod, (alphadenom + epsdenom)); // now c10::Scalar adds pairwisely to epsdenom
   
  }

  output = 1.0 / denom_prod * prefactor; // scalar multiplication

  return output;
}

/// Evaluation of a single pole in Fermi/Bose functions for a given `ami_vars`.
at::Tensor TamiBase::fermi_pole(ami_parms &parms, pole_struct pole,
                                         ami_vars external) {
  if (verbose) {
    std::cout << "Working on pole" << std::endl;
    print_pole_struct_info(pole);
  }

  at::Tensor output = at::zeros(external.energy_.size(0), at::kComplexDouble).to(device);
  int eta = 0;

  double beta = external.BETA_;
  double E_REG = parms.E_REG_;

  // TODO: should we use at::scalars (c10::complex<double>) instead of std::complex<double> -- Hopefully this is fixed
  // Spectral evaluation only
  TamiBase::complex_double freq_shift(0, 0);
  if (pole.x_alpha_.size() != 0) {
    for (int i = 0; i < pole.x_alpha_.size(); i++) {
      freq_shift = external.frequency_[i] * (double)pole.x_alpha_[i];
    }
  }

  // Future Development: 
  // In order to generalize to have fermi and bose lines, from here to 'sigma' needs
  // to be considered.

  // example fixed
  //
  // 1) create a stat map for the frequencies std::vector<int> stat_map:
  // populate 1 for fermi and 0 for bose. length is equal to alpha.size() 2)
  // simply replace eta=eta + 1*stat_map[i]
  //

  // alternate fix.  parms.TYPE_ is 0 for sigma, 1 for Pi etc.  So if
  // parms.TYPE_==1 and pole.alpha_.back()==1 (or -1), don't add one. else add
  // one to eta.


  for (int i = 0; i < pole.alpha_.size() - 1; i++) {
   
    if (pole.alpha_[i] != 0) {
      eta++;
    }
  }

  // handle external based on graph type
  if (pole.alpha_.back() != 0 && parms.TYPE_ != TamiBase::Pi_phuu &&
      parms.TYPE_ != TamiBase::Pi_phud && parms.TYPE_ != TamiBase::doubleocc &&
      parms.TYPE_ != TamiBase::Pi_ppuu && parms.TYPE_ != TamiBase::Pi_ppud &&
      parms.TYPE_ != TamiBase::FORCE) {
    eta++;
   
  }

  // if this is a double occupancy graph then the external line is bosonic. so
  // it is a bosonic matsubara integral. so eta needs to be incremented IF the
  // pole is for the last integration index
  double docc_sign=1;
  if (parms.TYPE_ == TamiBase::doubleocc || bosonic_external) {
    if (pole.index_ == pole.alpha_.size() - 1) {
      eta++;
      docc_sign=-1;
    
    }
  }

  // END TODO

  // could put infor into ami_vars external as to what the state type of the
  // external variables is.
  at::Tensor E = get_energy_from_pole(pole, external).to(device);

  // In the case of spectral poles the freq_shift might not be zero
  // TODO: Need to convert freqshift into a c10 complex or scalar to add pairwise or just always use them
  c10::Scalar freqshift = c10::complex<double>(freq_shift.real(), freq_shift.imag()); // cast to a tensor friendly complex value
  E = E + freqshift;


  double sigma = pow(-1.0, double(eta));

  TamiBase::complex_double zero(0, 0);
  TamiBase::complex_double im(0, 1);

  // TODO: Figure out what to do about the regulation code here - for now we proceed commented (E_REG is usually 0 anyways)

  /*
  // If energy denominator would be zero attempts to regulate if bosonic
  // function
  if (E == zero && sigma == -1) { // && pole.der_==0    Not sure if pole
    // derivative plays any role

    if (drop_bosonic_diverge) {

      return zero; // This is a dangerous approximation. 
    }
    E += E_REG;
  } else {
    if (sgn(E.real()) != 0) {
      E += E_REG * sgn(E.real());
    } else {
      E += E_REG;
    }
  }

  if (drop_der && pole.der_ != 0) {
    return zero;
  }
  */

  int m = pole.der_;

  // compute m'th derivative
  output = fermi_bose(m, sigma, beta, E);

/* 
  // TODO: print max here as well ? commented until proper printint formating is implemented
  if (verbose) {
    std::cout << "Fermi Pole Term evaluated to " << output << " at energy " << E
              << " with sigma " << sigma << " betaE is " << beta * E
              << " in exponent " << at::exp(beta * (E)) << std::endl;
			  
  // TODO: No idea how to print this .. I guess its only a verbose output so its fine to print the whole batch
	std::cout<<"Energy list is :(";
		for(int ii=0; ii< external.energy_.size(1); ii++){
			std::cout<<std::setprecision(20)<<" "<<external.energy_.index({torch::indexing::Slice(),a})<<" ,";
		}
		std::cout<<")"<<std::endl;
  }
*/

// TODO: Unclear if this is correct 
  if (parms.TYPE_ == TamiBase::doubleocc || bosonic_external) {
    output = docc_sign * output;
  }

  return output;
}

/// This computes the mth order derivative of the Fermi function or the
/// negative of the Bose distribution functions given by \f$\frac{1}{\sigma
/// \exp^{\beta E}+1} \f$ at \f$ \beta\f$, for energy \f$ E\f$. \f$
/// \sigma=1.0\f$ for Fermi and -1.0 for Bose.
at::Tensor TamiBase::fermi_bose(int m, double sigma, double beta,
                                         at::Tensor E) {

  int batch_size = E.size(0);
  at::Tensor output = at::zeros(batch_size, at::kComplexDouble).to(device);
  at::Tensor term = at::zeros(batch_size, at::kComplexDouble).to(device);

  if (m == 0) {
              // Note: Disabled on 09/09/2022 for overflow testing 
              // double arg = std::real(beta * E);
              // double arg_amp = std::abs(arg);

              // if (arg > exp_max_arg) {
                // double arg_sign = (double)sgn(arg);
                // if (arg_sign > 0) {
                  // output = 0;
                // } else {
                  // output = 1;
                // }
              // } else {
      output = 1.0 / (sigma * at::exp(beta * (E)) + 1.0);
                // }
    // return output;
  } else { // compute m'th derivative

    for (int k = 0; k < m + 1; k++) {
      // depricated: Original format
      // term = frk(m, k) * std::exp(k * beta * (E)) * std::pow(sigma, k) *
             // std::pow(-1.0, k + 1) /
             // std::pow(sigma * std::exp(beta * (E)) + 1.0, k + 1);
	// Reformatted with fewer exponentials
      term= frk(m,k)*std::pow(sigma,k) *std::pow(-1.0,
      k+1)*(1.0/(sigma*at::exp(beta*(E))+1.0)/at::pow(sigma
      +at::exp(-beta*(E)), k)) ;
      output += term;

      // TODO: fix the printing function then print the vector correctly
      /*
      if (verbose) {
        std::cout << "On k'th derivative " << k << std::endl;
        // TODO: Right now, this will print the whole tensor of numbers - not sure if imaginary part will print (usually shows an error)
        std::cout << "Fermi-bose Term evaluated to " << term << " at energy "
                  << E << " with sigma " << sigma << " betaE is " << beta * E
                  << " in exponent " << at::exp(beta * (E)) << std::endl;
      }
      */
    }

    output = output * std::pow(beta, m) * (-1.0);

    // TODO: figure out what to do about this ... it kind of defeats the purpose of the tensor. I guess we just use max(tensor) ... come back to this
    /*
    if ((std::abs(std::real(output)) > precision_cutoff)|| (std::abs(std::imag(output)) > precision_cutoff)) {
    
      overflow_detected = true;
    }
    */
  }

  return output;
}

at::Tensor TamiBase::get_energy_from_pole(pole_struct pole,
                                                   ami_vars external) {
  at::Tensor output = at::zeros(external.energy_.size(0), at::kComplexDouble).to(device);

  // Evaluating energies for pole
  for (int i = 0; i < pole.eps_.size(); i++) {
    output += double(pole.eps_[i]) * external.energy_.index({torch::indexing::Slice(),i});
  }

  return output;
}

/// Using notation to match https://doi.org/10.1103/PhysRevB.99.035120.
/// They produced coefficients to the fermi functions and put them in a table.
/// We derive a general expression for those coefficients - we believe this to
/// be general but have only checked up to 6th order.
double TamiBase::frk(int r, int k) {
  double output = 0.0;

  for (int m = 0; m < k + 1; m++) {
    output += binomialCoeff(k, m) * std::pow(m, r) * (std::pow(-1, k - m));
  }

  
  return output;
}


/// Returns value of Binomial Coefficient C(n, k).
int TamiBase::binomialCoeff(int n, int k) {
  int res = 1;

  // Since C(n, k) = C(n, n-k)
  if (k > n - k)
    k = n - k;

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}