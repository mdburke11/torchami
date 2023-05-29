#include "examples.hpp"


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
