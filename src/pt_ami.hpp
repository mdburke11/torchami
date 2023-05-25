/*
 * Copyright (C) 2018 JPF LeBlanc jleblanc@mun.ca  See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#pragma once

#include <algorithm>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

 AmiBase::g_prod_t construct_example2();
 AmiBase::ami_vars construct_ext_example2();
 
 AmiBase::ami_vars construct_4ord_ext_multipole_example();
 AmiBase::g_prod_t construct_multipole_example();
 
 AmiBase::ami_vars construct_ext_example1_bose();
 AmiBase::g_prod_t construct_example1_bose();
 
 AmiBase::g_prod_t construct_example_Y();
 AmiBase::ami_vars construct_ext_example_Y();
 AmiBase::g_prod_t construct_example_J();
 AmiBase::ami_vars construct_ext_example_J();

/**
 * @class PtAmiBase
 *
 *
 *
 * @brief  The primary class of libami with pole treee factorization 
 *
 * @note See https://github.com/jpfleblanc/libami
 *
 * @author James P.F. LeBlanc
 *
 * @version Revision: 0.8
 *
 * @date Date: 2023/02/20
 *
 *
 * Contact: jleblanc@mun.ca
 *
 *
 *
 *
 */
 
 
 class PtAmiBase {
public:
  /// Returns the sign of a value - or zero if it is uniquely zero.
  template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
  
  // template for checking equality of vectors
  template<typename T>
  bool isEqual(std::vector<T> const &v1, std::vector<T> const &v2)
  {
      return (v1.size() == v2.size() &&
              std::equal(v1.begin(), v1.end(), v2.begin()));
  }
  
  // this replaces all of the signs and poles used in previous AMI versions
  
  /** Term Structure for term-by-term evaluation.  Conceptually simpler than SPR
   * construction. Storage translates to \f$ \prod{f(p_i)}\prod{G_j}\times sign
   * \f$.
   *
   */
  struct ft_term {
    ft_term() {}

    ft_term(double s, FermiTree::fermi_tree_t ft, AmiBase::g_prod_t g) {
      sign_ = s;
      ft_ = ft;
      g_prod_ = g;
    }

    /// Sign prefactor
    double sign_ = 1;
    /// List of poles, \f$ \prod{f(p_i)}\f$.
    FermiTree::fermi_tree_t ft_;
    /// List of Green's functions, \f$ \prod{G_j}\f$.
    AmiBase::g_prod_t g_prod_;
  };  
  
  
  typedef std::vector<ft_term> ft_terms;
  typedef std::vector<FermiTree::fermi_tree_t> ft_list;
  
  void construct(int N_INT, AmiBase::g_prod_t R0, ft_terms &terms_out);
  void construct(AmiBase::ami_parms &parms, AmiBase::g_prod_t R0, ft_terms &terms_out);
  
  void factorize(ft_terms &in_terms, ft_terms &out_terms);
  bool g_prod_equiv(AmiBase::g_prod_t &gp1, AmiBase::g_prod_t &gp2, int &sign);
  // Need evaluate functions - worry about these later 
  
  /// Integrates a single Matsubara index.
  void integrate_step(int index, ft_terms &in_terms, ft_terms &out_terms);
  void term_integrate_step(int index, ft_term &in_term, ft_terms &out_terms);
  void terms_general_residue(ft_term &this_term, AmiBase::pole_struct this_pole,ft_terms &out_terms);
  
  void terms_to_ftterms(AmiBase::terms &in_terms, ft_terms &out_terms);
  
  void plist_to_ft(AmiBase::pole_array_t &plist, FermiTree::fermi_tree_t &ft);
                             
  void take_term_derivative(ft_term &in_term, AmiBase::pole_struct &pole, ft_terms &out_terms);
  
  /// Screen IO for debugging.
  void print_terms(ft_terms &t);
  /// Screen IO for debugging.
  void print_term(ft_term &t);
  
  std::string pretty_print_ft_term(ft_term &ft);
  std::string pretty_print_gprod(AmiBase::g_prod_t &gp);
  std::string pretty_print_g(AmiBase::g_struct &g);
  std::string pretty_print_ft_terms(ft_terms &fts);
  
  
  std::complex<double> eval_ft(AmiBase::ami_parms &parms, FermiTree::fermi_tree_t &ft1,  FermiTree::vertex_t &v, AmiBase::ami_vars &external);
  std::complex<double> evaluate_term(AmiBase::ami_parms &parms, ft_term &ft_term,
                                            AmiBase::ami_vars &external);
  std::complex<double> evaluate(AmiBase::ami_parms &parms, ft_terms &ft_terms,
                                            AmiBase::ami_vars &external);
  AmiBase amibase;
  FermiTree FT;
  
 
  
  
private:
};
  
 
  
  
  