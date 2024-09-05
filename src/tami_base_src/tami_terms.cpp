// Catchall source file for terms helper functions in the torchami development -
// contains functions from libami's ami_base_terms.cpp, ami_base_optimize.cpp,
// ami_base_terms_optimize.cpp

#include "tami_base.hpp"

void TamiBase::convert_terms_to_ri(terms &ami_terms, Ri_t &Ri) {
  Ri.clear();

  for (int i = 0; i < ami_terms.size(); i++) {
    Ri.push_back(ami_terms[i].g_list);
  }

  return;
}

void TamiBase::integrate_Mat_ind_step(int index, terms &in_terms,
                                      terms &out_terms) {
  out_terms.clear();

  for (int t_index = 0; t_index < in_terms.size(); t_index++) {
    pole_array_t poles;
    poles = find_poles(index, in_terms[t_index].g_list);

    for (int i = 0; i < poles.size(); i++) {
      if (poles[i].multiplicity_ == 1) {
        term new_term;
        new_term.g_list = simple_residue(in_terms[t_index].g_list, poles[i]);
        // take sign from original term and multiply by new one
        new_term.sign =
            in_terms[t_index].sign *
            get_simple_sign(index, in_terms[t_index].g_list, poles[i]);
        // take poles from originating term
        new_term.p_list = in_terms[t_index].p_list;
        new_term.p_list.push_back(poles[i]);

        out_terms.push_back(new_term);
      } else {
        terms new_terms;
        terms_general_residue(in_terms[t_index], poles[i], new_terms);

        // put new terms in the list
        out_terms.insert(out_terms.end(), new_terms.begin(), new_terms.end());
      }
    }
  }
}

void TamiBase::terms_general_residue(term &this_term, pole_struct this_pole,
                                     terms &out_terms) {
  out_terms.clear();

  double starting_sign;
  starting_sign =
      get_starting_sign(this_term.g_list, this_pole) * this_term.sign;

  term W;
  W.g_list = reduce_gprod(this_term.g_list, this_pole);

  W.p_list.push_back(this_pole);

  W.sign = starting_sign;
  // the W

  terms int_terms;
  int_terms.push_back(W);

  for (int m = 0; m < this_pole.multiplicity_ - 1; m++) {
    terms temp_terms;
    for (int i = 0; i < int_terms.size(); i++) {
      take_term_derivative(int_terms[i], this_pole, temp_terms);
    }

    int_terms = temp_terms;
  }

  // now we have all of the terms and their derivatives. so now it is safe to
  // sub in the poles

  for (int i = 0; i < int_terms.size(); i++) {
    for (int j = 0; j < int_terms[i].g_list.size(); j++) {
      int_terms[i].g_list[j] = update_G_pole(int_terms[i].g_list[j], this_pole);
    }

    // for every term have to put the pole list back

    int_terms[i].p_list.insert(int_terms[i].p_list.end(),
                               this_term.p_list.begin(),
                               this_term.p_list.end());
  }

  // if this gets slow should swap instead
  out_terms = int_terms;
}

void TamiBase::take_term_derivative(term &in_term, pole_struct &pole,
                                    terms &out_terms) {

  terms fd_terms;
  terms gd_terms;

  for (int one = 0; one < in_term.p_list.size(); one++) {
    term temp;
    temp.g_list = in_term.g_list;
    temp.p_list = in_term.p_list;
    temp.sign = in_term.sign;

    temp.p_list[one].der_++;

    fd_terms.push_back(temp);
  }

  // now do the gprod terms

  for (int i = 0; i < in_term.g_list.size(); i++) {
    int alpha = in_term.g_list[i].alpha_[pole.index_];

    term temp_gd_term;

    temp_gd_term.p_list = in_term.p_list;

    if (alpha != 0) {
      temp_gd_term.sign = -(double)alpha * in_term.sign;
      for (int m = 0; m < in_term.g_list.size(); m++) {
        if (i == m) {
          temp_gd_term.g_list.push_back(in_term.g_list[m]);
          temp_gd_term.g_list.push_back(in_term.g_list[m]);
        } else {
          temp_gd_term.g_list.push_back(in_term.g_list[m]);
        }
      }

      gd_terms.push_back(temp_gd_term);
    }
  }

  out_terms.insert(out_terms.end(), fd_terms.begin(), fd_terms.end());

  out_terms.insert(out_terms.end(), gd_terms.begin(), gd_terms.end());
}

void TamiBase::print_pole_struct_info(pole_struct g) {
  std::cout << "Eps=(";
  print_epsilon_info(g.eps_);
  std::cout << ")";
  std::cout << std::endl;
  std::cout << "Alpha=(";
  print_alpha_info(g.alpha_);
  std::cout << ")";
  std::cout << "X_Alpha=(";
  print_alpha_info(g.x_alpha_);
  std::cout << ")";
  std::cout << std::endl;
  std::cout << "Der: " << g.der_ << std::endl;
}

void TamiBase::print_epsilon_info(TamiBase::epsilon_t eps) {

  for (std::vector<int>::iterator it = eps.begin(); it != eps.end(); ++it) {
    std::cout << *it << ' ';
  }
}

void TamiBase::print_alpha_info(TamiBase::alpha_t alpha) {
  for (std::vector<int>::iterator it = alpha.begin(); it != alpha.end(); ++it) {
    std::cout << *it << ' ';
  }
}

void TamiBase::factorize_Rn(Ri_t &Rn, g_prod_t &unique_g, R_ref_t &Rref,
                            ref_eval_t &Eval_list) {
  unique_g.clear();
  Rref.clear();
  Eval_list.clear();

  // i is the entry of Rn
  // first step is to generate
  for (int i = 0; i < Rn.size(); i++) {
    ref_v_t ref_v;
    // Rn[i] is a g_prod_t, so Rn[i][j] is a G -> g_struct
    for (int j = 0; j < Rn[i].size(); j++) {
      ref_t ref;

      // if unique list is empty then add entry to it.
      if (unique_g.size() == 0) {
        unique_g.push_back(Rn[i][j]);
        ref = std::make_pair(0, 1);
        ref_v.push_back(ref);
        continue;
      }
      // for each G search through unique_g
      bool found = false;
      for (int m = 0; m < unique_g.size(); m++) {
        int this_sign = 0;
        bool eq = g_struct_equiv(Rn[i][j], unique_g[m], this_sign);

        if (eq) {
          ref = std::make_pair(m, this_sign);
          ref_v.push_back(ref);
          found = true;
          break; // break out of the unique_g loop since we found an
                 // equivalent Rn[i][j] g_struct
        }
      }

      // if made it here it means this is a unique G and so needs to be
      // stored and a ref created
      if (!found) {
        unique_g.push_back(Rn[i][j]);
        ref = std::make_pair(unique_g.size() - 1, 1);
        ref_v.push_back(ref);
      }
    }

    if (ref_v.size() != 0) {
      Rref.push_back(ref_v);
    }
  }

  reduce_rref(Rref, Eval_list);
}

bool TamiBase::g_struct_equiv(g_struct &g1, g_struct &g2, int &sign) {

  sign = 0; // by default zero sign means they are not equiv

  bool result = true;

  if (g1.eps_.size() != g2.eps_.size()) {
    return false;
  }
  if (g1.alpha_.size() != g2.alpha_.size()) {
    return false;
  }

  if (g1.species_ != g2.species_) {
    return false;
  }

  std::vector<int> signs;

  for (int i = 0; i < g1.eps_.size(); i++) {
    if (std::abs(g1.eps_[i]) != std::abs(g2.eps_[i])) {
      return false;
      break;
    } else {
      if (g1.eps_[i] != 0) {
        signs.push_back(mathUtils::sgn(g1.eps_[i] / g2.eps_[i]));
      }
    }
  }

  for (int i = 0; i < g1.alpha_.size(); i++) {
    if (std::abs(g1.alpha_[i]) != std::abs(g2.alpha_[i])) {
      return false;
      break;
    } else {
      if (g1.alpha_[i] != 0) {
        signs.push_back(mathUtils::sgn(g1.alpha_[i] / g2.alpha_[i]));
      }
    }
  }

  // in principle if the code gets here the size of signs is >0 so don't need
  // to check that

  if (std::adjacent_find(signs.begin(), signs.end(),
                         std::not_equal_to<int>()) == signs.end()) {
    sign = signs.back();
  } else {
    return false;
  }

  return result;
}

// this version kicks out the extra Rref entries
void TamiBase::reduce_rref(R_ref_t &Rref, ref_eval_t &Eval_list) {
  std::vector<int> used;

  for (int i = 0; i < Rref.size(); i++) {

    ref_v_t this_vec;
    if (std::find(used.begin(), used.end(), i) != used.end()) {
      continue;
    } else {
      used.push_back(i);
      int this_sign = 1;
      for (int pair = 0; pair < Rref[i].size(); pair++) {

        this_sign = this_sign * Rref[i][pair].second;
      }
      ref_t this_ref = std::make_pair(i, this_sign);
      this_vec.push_back(this_ref);
    }

    for (int j = i + 1; j < Rref.size(); j++) {

      if (std::find(used.begin(), used.end(), j) != used.end()) {
        continue;
      } else {
        int sign1, sign2;
        bool ditto = pair_v_equiv(Rref[i], Rref[j], sign1, sign2);

        if (ditto) {
          used.push_back(j);
          ref_t this_ref = std::make_pair(j, sign2);

          this_vec.push_back(this_ref);
        }
      }
    }
    // don't want empty though don't think it is possible
    if (this_vec.size() > 0) {
      Eval_list.push_back(this_vec);
    }
  }

  // now make a minimal Rref_v

  R_ref_t newref;

  for (int i = 0; i < Eval_list.size(); i++) {
    ref_v_t newentry;

    for (int j = 0; j < Rref[Eval_list[i][0].first].size(); j++) {
      ref_t this_pair;
      this_pair.first = Rref[Eval_list[i][0].first][j].first;
      this_pair.second = 1;
      newentry.push_back(this_pair);
    }

    newref.push_back(newentry);
  }

  Rref = newref;
}

bool TamiBase::pair_v_equiv(ref_v_t &r1, ref_v_t &r2, int &r1sign,
                            int &r2sign) {
  r1sign = 0;
  r2sign = 0;

  if (r1.size() != r2.size()) {
    return false;
  }

  int this_sign = 1;
  int this_sign2 = 1;
  for (int i = 0; i < r1.size(); i++) {
    if (r1[i].first != r2[i].first) {
      return false;
    }
    this_sign = this_sign * r1[i].second;
    this_sign2 = this_sign2 * r2[i].second;
  }

  r1sign = this_sign;
  r2sign = this_sign2;
  return true;
}
