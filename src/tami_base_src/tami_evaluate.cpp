#include "tami_base.hpp"

at::Tensor TamiBase::evaluate(TamiBase::ami_parms &parms, ft_terms &ft_terms,
                              TamiBase::ami_vars &external) {

  int ebatch_size = external.energy_.size(0);
  int fbatch_size = external.frequency_.size(0);

  at::Tensor term =
      at::zeros({fbatch_size, ebatch_size},
                options); // dummy intialization on device in options
  at::Tensor output = at::zeros({fbatch_size, ebatch_size}, options);

  for (int i = 0; i < ft_terms.size(); i++) {
    term = evaluate_term(parms, ft_terms[i], external);
    output += term;
  }

  return output;
}

at::Tensor TamiBase::evaluate_term(TamiBase::ami_parms &parms, ft_term &ft_term,
                                   TamiBase::ami_vars &external) {

  at::Tensor gprod, fprod;

  gprod =
      eval_gprod(parms, ft_term.g_prod_, external); // n_freq by n_energy tensor
  TamiBase::FermiTree::vertex_t r = FT.get_root(ft_term.ft_);
  fprod = eval_ft(parms, ft_term.ft_, r,
                  external); // 1 by n_energy tensor - it needs to be blown up
                             // to match the dimension of g_prod

  int fbatch_size = external.frequency_.size(0);
  at::Tensor blownupfprod = fprod.repeat(
      {fbatch_size, 1}); // now n_freq by n_energy and can be multiplied
  at::Tensor output = ft_term.sign_ * torch::multiply(gprod, blownupfprod);

  return output;
}

// TamiBase::complex_double TamiBase::eval_fprod(ami_parms &parms, pole_array_t
// &p_list, ami_vars &external)

at::Tensor TamiBase::eval_ft(TamiBase::ami_parms &parms,
                             TamiBase::FermiTree::fermi_tree_t &ft1,
                             TamiBase::FermiTree::vertex_t &v,
                             TamiBase::ami_vars &external) {

  std::vector<TamiBase::FermiTree::vertex_t> level;

  FT.get_next_level(ft1, v, level);

  int ebatch_size = external.energy_.size(0);
  at::Tensor output = at::zeros({1, ebatch_size}, options);
  at::Tensor aoutput = at::zeros({1, ebatch_size}, options);
  at::Tensor moutput = at::ones({1, ebatch_size}, options);

  switch (ft1[v].operation_) {

  case 0:
    for (int i = 0; i < level.size(); i++) {

      aoutput += eval_ft(parms, ft1, level[i], external) * ft1[v].prefactor_;
    }
    output = aoutput;
    break;
  case 1:
    for (int i = 0; i < level.size(); i++) {
      // std::cout<<"On evaluate prefactor is "<<ft1[v].prefactor_<<" on vertex
      // "<< ft1[v].index_<< std::endl;
      moutput =
          moutput * eval_ft(parms, ft1, level[i],
                            external); // added in prefactor to the mult term.
    }
    moutput = moutput * ft1[v].prefactor_; // only one instance of the prefactor
                                           // if there is a multiplication
    output = moutput;
    break;
  case 2:

    ft1[v].value_ = fermi_pole(parms, ft1[v].pole_, external);

    output = ft1[v].value_ * ft1[v].prefactor_;
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

  int ebatch_size = external.energy_.size(0);
  int fbatch_size = external.frequency_.size(0);

  at::Tensor output =
      at::zeros({fbatch_size, ebatch_size},
                options); // get len(energies) * len(frequencies) answers

  at::Tensor denom_prod = at::ones({fbatch_size, ebatch_size}, options);
  double prefactor = external.prefactor;

  double E_REG = parms.E_REG_;
  bool verbose = false;

  for (int i = 0; i < g_prod.size(); i++) {
    at::Tensor alphadenom = at::zeros({1, fbatch_size}, options);
    at::Tensor epsdenom = at::zeros({1, ebatch_size}, options);

    for (int a = 0; a < g_prod[i].alpha_.size(); a++) {
      alphadenom +=
          double(g_prod[i].alpha_[a]) *
          external.frequency_.index({torch::indexing::Slice(),
                                     a}); // I think we should have the eps and
                                          // alpha on GPU as well - TODO!
    }

    for (int a = 0; a < g_prod[i].eps_.size(); a++) {
      epsdenom +=
          double(g_prod[i].eps_[a]) *
          external.energy_.index({torch::indexing::Slice(),
                                  a}); // I think we should have the eps and
                                       // alpha on GPU as well - TODO!
    }

    // Now do "Matrix multiplication" but with addition opperation for these
    // vectors: eg epsdenom = [E1, E2, E3], alphadenom = [A1, A2, A3],
    // {alphadenom^T} "*" {epsdenom} = [[A1+E1, A1+E2, A1+E3], [A2+E1, A2+E2,
    // A2+E3], [A3+E1, A3+E2, A3+E3]] This evaluates all frequencies at all the
    // energies in the batch simultaneously

    at::Tensor blownupEps = epsdenom.repeat({fbatch_size, 1});
    at::Tensor blownupAlpha =
        alphadenom.transpose(0, 1).repeat({1, ebatch_size});

    denom_prod = torch::multiply(
        denom_prod,
        (blownupEps + blownupAlpha)); // Still the same pairwise multiply but
                                      // now there is n_f * n_e results that
                                      // then need to be added down the columns
  }

  output = 1.0 / denom_prod *
           prefactor; // scalar multiplication: 1/tensor = [1/t_ij]_ij

  return output;
}

/// Evaluation of a single pole in Fermi/Bose functions for a given `ami_vars`.
at::Tensor TamiBase::fermi_pole(ami_parms &parms, pole_struct pole,
                                ami_vars external) {
  if (verbose) {
    std::cout << "Working on pole" << std::endl;
    print_pole_struct_info(pole);
  }

  int ebatch_size = external.energy_.size(0);

  at::Tensor output =
      at::zeros({1, ebatch_size}, options); // get 1 by len(energies) answers
  int eta = 0;

  double beta = external.BETA_;
  double E_REG = parms.E_REG_;

  // Future Development:
  // In order to generalize to have fermi and bose lines, from here to 'sigma'
  // needs to be considered.

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
  double docc_sign = 1;
  if (parms.TYPE_ == TamiBase::doubleocc || bosonic_external) {
    if (pole.index_ == pole.alpha_.size() - 1) {
      eta++;
      docc_sign = -1;
    }
  }

  // END TODO

  // could put infor into ami_vars external as to what the state type of the
  // external variables is.
  at::Tensor E = get_energy_from_pole(
      pole, external); // list of eval'd lin combs of energys

  /*

  // No longer needed because the freq_shift is a tensor of type
  TamiBase::complex_double

  // In the case of spectral poles the freq_shift might not be zero
  // TODO: Need to convert freqshift into a c10 complex or scalar to add
  pairwise or just always use them c10::Scalar freqshift =
  c10::complex<double>(freq_shift.real(), freq_shift.imag()); // cast to a
  tensor friendly complex value
  */

  double sigma = pow(-1.0, double(eta));

  TamiBase::complex_double zero(0, 0);
  TamiBase::complex_double im(0, 1);

  // TODO: Figure out what to do about the regulation code here - for now we
  // proceed commented (E_REG is usually 0 anyways)

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
    if (mathUtils::sgn(E.real()) != 0) {
      E += E_REG * mathUtils::sgn(E.real());
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
    // TODO: print max here as well ? commented until proper printint formating
    is implemented if (verbose) { std::cout << "Fermi Pole Term evaluated to "
    << output << " at energy " << E
                << " with sigma " << sigma << " betaE is " << beta * E
                << " in exponent " << at::exp(beta * (E)) << std::endl;

    // TODO: No idea how to print this .. I guess its only a verbose output so
    its fine to print the whole batch std::cout<<"Energy list is :("; for(int
    ii=0; ii< external.energy_.size(1); ii++){
                          std::cout<<std::setprecision(20)<<"
    "<<external.energy_.index({torch::indexing::Slice(),a})<<" ,";
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

/// This is the original fermi_bose function that is not necessarily numerically stable
// at::Tensor TamiBase::fermi_bose(int m, double sigma, double beta,
//                                 at::Tensor E) {

//   int ebatch_size = E.size(1);
//   at::Tensor output = at::zeros({1, ebatch_size}, options);
//   at::Tensor term = at::zeros({1, ebatch_size}, options);

//   if (m == 0) {
    
//     output = 1.0 / (sigma * at::exp(beta * (E)) + 1.0);
    
//   } else { // compute m'th derivative

    // for (int k = 0; k < m + 1; k++) {
      
    //   // Reformatted with fewer exponentials
    //   term = mathUtils::frk(m, k) * std::pow(sigma, k) * std::pow(-1.0, k + 1) *
    //          (1.0 / (sigma * at::exp(beta * (E)) + 1.0) /
    //           at::pow(sigma + at::exp(-beta * (E)), k));
    //   output += term;

      
    // }

    // output = output * std::pow(beta, m) * (-1.0);

   
//   }

//   return output;
// }

at::Tensor TamiBase::fermi_bose(int m, double sigma, double beta, at::Tensor E){

if(m<3){
return fermi_bose_mlte2(m,sigma,beta,E);
}
else{

  return fermi_bose_agrad(m,sigma,beta,E);
}
}

/// This is a version that should be numerically stable based on stable primitives. 
at::Tensor TamiBase::fermi_bose_mlte2(int m, double sigma, double beta, at::Tensor E) {
    TORCH_CHECK(E.is_complex(), "Input tensor E must be complex-valued");

    at::Tensor x = beta * E;
    at::Tensor x_real = at::real(x);  // Use only real part for distribution function

    at::Tensor out;

    if (m == 0) {
        // Fermi: f = 1 / (1 + exp(beta * E)) = sigmoid(-x)
        // Bose: f = -1 / (exp(beta * E) - 1)
        if (sigma == 1.0) {
            out = 1.0 / (1.0 + at::exp(x_real));  // Fermi-Dirac
        } else {
            out = -1.0 / (at::exp(x_real) - 1.0); // Negative Bose-Einstein
        }
    } else if (m == 1) {
        if (sigma == 1.0) {
            // f' = -β / (4 cosh²(βx/2))
            at::Tensor half_x = 0.5 * x_real;
            out = -beta / (4.0 * at::pow(at::cosh(half_x), 2));
        } else {
            // f' = β / (4 sinh²(βx/2))
            at::Tensor half_x = 0.5 * x_real;
            out = beta / (4.0 * at::pow(at::sinh(half_x), 2));
        }
    } else if (m == 2) {
        if (sigma == 1.0) {
            // f'' = β² sinh(βx/2) / (4 cosh³(βx/2))
            at::Tensor half_x = 0.5 * x_real;
            out = beta * beta * at::sinh(half_x) / (4.0 * at::pow(at::cosh(half_x), 3));
        } else {
            // f'' = -β² cosh(βx/2) / (4 sinh³(βx/2))
            at::Tensor half_x = 0.5 * x_real;
            out = -beta * beta * at::cosh(half_x) / (4.0 * at::pow(at::sinh(half_x), 3));
        }
    } else {
        TORCH_CHECK(false, "Only m = 0, 1, or 2 is supported in fermi_bose()");
    }

    return out.to(E.options());
}

// /// @brief  This is a more advanced auto differentiation version of the fermi_bose function. It should be tested more thoroughly for performance guarantees.
// /// @param m 
// /// @param sigma 
// /// @param beta 
// /// @param E 
// /// @return 
at::Tensor TamiBase::fermi_bose_agrad(int m, double sigma, double beta, at::Tensor E) {
    TORCH_CHECK(m >= 0, "Derivative order must be non-negative");

    // Ensure E has correct type and grad enabled
    at::Tensor E_real = at::real(E);
    E_real = E_real.detach().clone().set_requires_grad(true);
    at::Tensor x = beta * E_real;

    at::Tensor f0;
    if (sigma == 1.0) {
        f0 = 1.0 / (1.0 + at::exp(x));  // Fermi
    } else {
        f0 = -1.0 / (at::exp(x) - 1.0); // negative Bose
    }

    at::Tensor result = f0;
    for (int i = 0; i < m; ++i) {
        std::vector<at::Tensor> grads = torch::autograd::grad(
            /* outputs = */ {result},
            /* inputs  = */ {E_real},
            /* grad_outputs = */ {torch::ones_like(E_real)},
            /* retain_graph = */ true,
            /* create_graph = */ true
        );
        result = grads[0];
    }

    return result.detach().to(E.options());
}



// at::Tensor TamiBase::fermi_bose(int m, double sigma, double beta, at::Tensor E) {
//   at::Tensor x = beta * E;
//   at::Tensor sx = -sigma * x;

//   // Sigmoid(-σx) is the stable Fermi function
//   at::Tensor f = torch::sigmoid(sx);

//   if (m == 0) {
//     // fermi: f = sigmoid(-x), bose: f = -sigmoid(x)
//     return f;
//   }

//   if (m == 1) {
//     // Derivative: -β f (1 - f) for Fermi; β f (1 + f) for Bose
//     at::Tensor deriv;
//     if (sigma > 0) {
//       deriv = -beta * f * (1.0 - f);  // Fermi case
//     } else {
//       deriv = beta * f * (1.0 + f);   // Bose case
//     }
//     return deriv;
//   }

//   throw std::runtime_error("Only m=0 or m=1 supported in stable fermi_bose.");
// }



// at::Tensor TamiBase::fermi_bose(int m, double sigma, double beta,
//                                 at::Tensor E) {
//   constexpr double EXP_CUTOFF = 500.0;

//   at::Tensor x = beta * E;
//   at::Tensor x_real = at::real(x);  // use this for all comparisons

//   int ebatch_size = E.size(1);
//   at::Tensor output;

//   if (m == 0) {
//     if (sigma == 1.0) {
//       // Fermi case: numerically stable
//       output = torch::sigmoid(-x);
//     } else {
//       // Bose case: piecewise handling to prevent overflow
//       at::Tensor mask_large_pos = (x_real > EXP_CUTOFF);
//       at::Tensor mask_large_neg = (x_real < -EXP_CUTOFF);
//       at::Tensor mask_safe = ~(mask_large_pos | mask_large_neg);

//       output = at::zeros_like(E);
//       output.masked_fill_(mask_large_pos, 0.0);    // lim x→∞: f_B → 0
//       output.masked_fill_(mask_large_neg, -1.0);   // lim x→-∞: f_B → -1

//       if (mask_safe.any().item<bool>()) {
//         at::Tensor x_safe = x.masked_select(mask_safe);
//         at::Tensor val_safe = -1.0 / (1.0 + sigma * at::exp(x_safe));
//         output.masked_scatter_(mask_safe, val_safe);
//       }
//     }

//   } else {
//     // Higher-order derivatives of Fermi/Bose
//     output = at::zeros({1, ebatch_size}, E.options());

//     for (int k = 0; k < m + 1; k++) {
//       double prefactor = mathUtils::frk(m, k) * std::pow(sigma, k) * std::pow(-1.0, k + 1);

//       at::Tensor mask_large_pos = (x_real > EXP_CUTOFF);
//       at::Tensor mask_large_neg = (x_real < -EXP_CUTOFF);
//       at::Tensor mask_safe = ~(mask_large_pos | mask_large_neg);

//       at::Tensor term = at::zeros({1, ebatch_size}, E.options());

//       if (mask_safe.any().item<bool>()) {
//         at::Tensor x_safe = x.masked_select(mask_safe);
//         at::Tensor exp_x = at::exp(x_safe);
//         at::Tensor exp_neg_x = at::exp(-x_safe);

//         at::Tensor numerator = 1.0 / (sigma * exp_x + 1.0);
//         at::Tensor denom = at::pow(sigma + exp_neg_x, k);
//         at::Tensor val = prefactor * numerator / denom;

//         term.masked_scatter_(mask_safe, val);
//       }

//       // Optional: ensure defined behavior in extreme regions
//       term.masked_fill_(mask_large_pos | mask_large_neg, 0.0);

//       output += term;
//     }

//     output *= std::pow(beta, m) * (-1.0);
//   }

//   return output;
// }



at::Tensor TamiBase::get_energy_from_pole(pole_struct pole, ami_vars external) {
  at::Tensor output = at::zeros({1, external.energy_.size(0)}, options);

  // Evaluating energies for pole
  for (int i = 0; i < pole.eps_.size(); i++) {
    output += double(pole.eps_[i]) *
              external.energy_.index({torch::indexing::Slice(), i});
  }

  return output;
}

/// Using notation to match https://doi.org/10.1103/PhysRevB.99.035120.
/// They produced coefficients to the fermi functions and put them in a table.
/// We derive a general expression for those coefficients - we believe this to
/// be general but have only checked up to 6th order.
double mathUtils::frk(int r, int k) {
  double output = 0.0;

  for (int m = 0; m < k + 1; m++) {
    output += mathUtils::binomialCoeff(k, m) * std::pow(m, r) * (std::pow(-1, k - m));
  }

  return output;
}

/// Returns value of Binomial Coefficient C(n, k).
int mathUtils::binomialCoeff(int n, int k) {
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