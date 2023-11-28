#include "examples.hpp"


int main( int argc , char *argv[] )
{
  int mode=1;
  if (argc >= 2){
    std::istringstream ss(argv[1]);

    if (!(ss >> mode)){
      std::cerr << "Using default" << std::endl;
      mode = 99;
    }
  }

  switch(mode) {
    case 1:
      example2();
      break;
    case 2:
      //example1_bose();
      break;
    case 3:
      example4();
      break;
    case 4:
      example9();
      break;
    case 5:
      graph_library_example();
      break;
    case 6:
      renorm_PT_graph_example();
      break;
    default:
      example2();
      //example1_bose(); // not working
      example4();
      example6();
      example9();
      break;
  }
}

void default_example(){

  at::Device myDev = at::kCPU;//at::kCUDA;
  int freq_batchsize = 1;
  int energy_batchsize = 10;

  TamiBase PT(myDev);
  TamiBase::ft_terms ftout;

  TamiBase::g_prod_t R02=construct_multipole_example();//construct_example_J();//construct_multipole_example();//construct_example1_bose();
  TamiBase::ami_vars avars2=construct_4ord_ext_multipole_example(PT, energy_batchsize, freq_batchsize);//construct_ext_example_J();//construct_4ord_ext_multipole_example();//construct_ext_example1_bose();

  std::cout << "Starting energy tensor: " << std::endl;
  std::cout << format_r2_tensor(avars2.energy_) << std::endl;

  TamiBase::ft_terms ftout2;

  PT.construct(4, R02, ftout2);
  TamiBase::ami_parms parms2(0, 0);

  auto t2=std::chrono::high_resolution_clock::now();
  at::Tensor result2=PT.evaluate(parms2,ftout2,avars2);

  auto t_end=std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> diff2=t_end-t2;

  // std::chrono::nanoseconds d=std::chrono::duration_cast<std::chrono::nanoseconds>(diff1);

  
  std::cout<<"Result 4th order MP "<<format_r2_tensor(result2)<<std::endl;
  std::chrono::nanoseconds d2=std::chrono::duration_cast<std::chrono::nanoseconds>(diff2);
  std::cout<<"Evaluation took "<< d2.count()<<" nanoseconds"<<std::endl;	

}

void graph_library_example(){

  int fbatch_size = 1;
  int ebatch_size = 100;

  // this is an example of using the libtami_graph library to evaluate the graphs
  // in the ggm_examples folder.

  std::string foldername = "../ggm_examples"; // folder that holds all the graph files

  TamiBase::graph_type g_type = TamiBase::Sigma; // Specify graph type: self energy diagram

  int seed = 0;

  TamiGraph g(g_type, seed); // initiate the TamiGraph object with intial seed

  TamiGraph::gg_matrix_t ggm;
  int max_ord = 6; // max order diagram in the graph folder if order 6 by default

  g.read_ggmp(foldername, ggm, max_ord); 

  g.print_ggm(ggm);

  std::cout<<"Now lets label the graphs "<<std::endl;
	g.ggm_label(ggm,0);
	std::cout<<"All done! That was easy!"<<std::endl;

  TamiBase::g_prod_t R0;

  for (int i=2; i<7; i+=2){ // looping through the orders 2, 4, 6
    g.graph_to_R0(ggm[i][0].graph_vec[0], R0); // converting each graph into R0
    std::cout<< "\norder: " << i << std::endl;
    std::cout<< "alpha: " << std::endl;

    for (auto x : R0){ // Print the alphas
      for (int j=0; j < x.alpha_.size(); ++j){
        std::cout << x.alpha_[j] << " ";
      }
      std::cout << std::endl;
    }

    std::cout<< "epsilon: " << std::endl; // Print the epsilons
    for (auto x : R0){
      for (int j=0; j < x.eps_.size(); ++j){
        std::cout << x.eps_[j] << " ";
      }
      std::cout << std::endl;
    }
  }
}

void renorm_PT_graph_example(){

  int fbatch_size = 1;
  int ebatch_size = 100;

  // this is an example of using the libtami_graph library to evaluate the graphs
  // in the ggm_examples folder.

  std::string foldername = "../ggm_examples"; // folder that holds all the graph files

  TamiBase::graph_type g_type = TamiBase::Sigma; // Specify graph type: self energy diagram

  int seed = 0;

  TamiGraph g(g_type, seed); // initiate the TamiGraph object with intial seed

  TamiGraph::gg_matrix_t ggm;
  int max_ord = 2; // Only do the second order diagram's Ct diagrams
  int max_insertions = 1; // include all diagrams up to 2 insertions

  g.read_ggmp(foldername, ggm, max_ord); 

  g.print_ggm(ggm);

  std::cout<<"Now lets label the graphs "<<std::endl;
	g.ggm_label(ggm,0);
	std::cout<<"All done! That was easy!"<<std::endl;

  std::vector<TamiGraph::graph_t> temp_ct;
  std::vector<TamiBase::g_prod_t> R0_ct; // catch all of the counter term R0's
  TamiBase::g_prod_t R0;

  TamiGraph::graph_t second_ord_Sigma = ggm[2][0].graph_vec[0];
  g.generate_sigma_ct(second_ord_Sigma, temp_ct, max_insertions);
  
  for (auto x : temp_ct){
    g.graph_to_R0(x, R0);
    double pf = g.get_prefactor(x, 2); // fermionic loop counter
    std::cout<< "Prefactor: " << pf << std::endl;
    R0_ct.push_back(R0);
  }

  for (auto x : R0_ct){ // looping through all the CT R0's
    std::cout<< "\nalpha: " << std::endl;

    for (auto y : x){ // Print the alphas
      for (int j=0; j < y.alpha_.size(); ++j){
        std::cout << y.alpha_[j] << " ";
      }
      std::cout << std::endl;
    }

    std::cout<< "epsilon: " << std::endl; // Print the epsilons
    for (auto y : x){
      for (int j=0; j < y.eps_.size(); ++j){
        std::cout << y.eps_[j] << " ";
      }
      std::cout << std::endl;
    }
  }
}

void example2(){

// Number of parameters to evaluate - number of frequencies and num of energies to eval at a time
int fbatch_size = 1;
int ebatch_size = 10;

std::cout<<std::endl<<"-_-_-_ Example - Second Order _-_-_-"<<std::endl<<std::endl;	
	
//START Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
// class instance
at::Device myDev = at::kCPU;
TamiBase ami(myDev);

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_example2(); // Sets initial integrand 
TamiBase::ami_vars avars=construct_ext_example2(ami, ebatch_size, fbatch_size); // Sets 'external' parameter values 

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=2;  // Number of Matsubara sums to perform
TamiBase::ami_parms test_amiparms(N_INT, E_REG);

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

//simplified storage type 
TamiBase::ft_terms amiterms;

// Construct solution for problem defined in R0
ami.construct(N_INT, R0, amiterms);

	//timing info 
	auto t2=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff1=t2-t1;
	std::chrono::nanoseconds d1=std::chrono::duration_cast<std::chrono::nanoseconds>(diff1);

std::cout<<"Construction took "<<d1.count()<<" nanoseconds"<<std::endl;

	//timing info 
	auto t3=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
at::Tensor term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'. Note that the test_amiparms is the same as the first case 
	
	//timing info
	auto t4=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::chrono::nanoseconds d2=std::chrono::duration_cast<std::chrono::nanoseconds>(diff2);

std::cout<<"Term result was "<< format_r2_tensor(term_val)<< "Shape: " << term_val.sizes() <<std::endl;
std::cout<<"Evaluation took "<<d2.count()<<" nanoseconds"<<std::endl;
std::cout<<"Length of terms object was " << amiterms.size() << std::endl;
	
}

void example1_bose(){

// Number of parameters to evaluate - number of frequencies and num of energies to eval at a time
int fbatch_size = 1;
int ebatch_size = 100;

std::cout<<std::endl<<"-_-_-_ Example - First Order For Bosonic _-_-_-"<<std::endl<<std::endl;	
	
//START Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
// class instance
at::Device myDev = at::kCPU;
TamiBase ami(myDev);

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_example1_bose(); // Sets initial integrand 
TamiBase::ami_vars avars=construct_ext_example1_bose(ami, ebatch_size, fbatch_size); // Sets 'external' parameter values 

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=1;  // Number of Matsubara sums to perform
TamiBase::graph_type bose=TamiBase::Pi_phuu;
TamiBase::ami_parms test_amiparms(N_INT, E_REG, bose);

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

//simplified storage type 
TamiBase::ft_terms amiterms;

// Construct solution for problem defined in R0
ami.construct(N_INT, R0, amiterms);

	//timing info 
	auto t2=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff1=t2-t1;
	std::chrono::nanoseconds d1=std::chrono::duration_cast<std::chrono::nanoseconds>(diff1);

std::cout<<"Construction took "<<d1.count()<<" nanoseconds"<<std::endl;

	//timing info 
	auto t3=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
at::Tensor term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'. Note that the test_amiparms is the same as the first case 
	
	//timing info
	auto t4=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::chrono::nanoseconds d2=std::chrono::duration_cast<std::chrono::nanoseconds>(diff2);

std::cout<<"Term result was "<< format_r2_tensor(term_val)<<std::endl;
std::cout<<"Evaluation took "<<d2.count()<<" nanoseconds"<<std::endl;
std::cout<<"Length of terms object was " << amiterms.size() << std::endl;

 }


void example4(){

// Number of parameters to evaluate - number of frequencies and num of energies to eval at a time
int fbatch_size = 1;
int ebatch_size = 100;

std::cout<<std::endl<<"-_-_-_ Example - Fourth Order _-_-_-"<<std::endl<<std::endl;	
	
//START Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
// class instance
at::Device myDev = at::kCUDA;
TamiBase ami(myDev);

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_multipole_example(); // Sets initial integrand 
TamiBase::ami_vars avars=construct_4ord_ext_multipole_example(ami, ebatch_size, fbatch_size); // Sets 'external' parameter values 

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=4;  // Number of Matsubara sums to perform
TamiBase::ami_parms test_amiparms(N_INT, E_REG);

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

//simplified storage type 
TamiBase::ft_terms amiterms;

// Construct solution for problem defined in R0
ami.construct(N_INT, R0, amiterms);

	//timing info 
	auto t2=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff1=t2-t1;
	std::chrono::nanoseconds d1=std::chrono::duration_cast<std::chrono::nanoseconds>(diff1);

std::cout<<"Construction took "<<d1.count()<<" nanoseconds"<<std::endl;

	//timing info 
	auto t3=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
at::Tensor term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'.
	
	//timing info
	auto t4=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::chrono::nanoseconds d2=std::chrono::duration_cast<std::chrono::nanoseconds>(diff2);

std::cout<<"Term result was "<< format_r2_tensor(term_val)<<std::endl;
std::cout<<"Evaluation took "<<d2.count()<<" nanoseconds"<<std::endl;
std::cout<<"Length of terms object was " << amiterms.size() << std::endl;
 }


  void example6(){

  // Number of parameters to evaluate - number of frequencies and num of energies to eval at a time
  int fbatch_size = 1;
  int ebatch_size = 100;

  std::cout<<std::endl<<"-_-_-_ Example - Sixth Order _-_-_-"<<std::endl<<std::endl;	
	
  //START Example 
  // Same Problem, using ami_term storage type
  std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
  // class instance
  at::Device myDev = at::kCPU;
  TamiBase ami(myDev);

  // Problem setup (see ami_example.cpp)
  TamiBase::g_prod_t R0=construct_example6(); // Sets initial integrand 
  TamiBase::ami_vars avars=construct_ext_example6(ami, ebatch_size, fbatch_size); // Sets 'external' parameter values 

  // Integration/Evaluation parameters
  double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
  int N_INT=6;  // Number of Matsubara sums to perform
  TamiBase::ami_parms test_amiparms(N_INT, E_REG);

    //timing info
    auto t1=std::chrono::high_resolution_clock::now();

  //simplified storage type 
  TamiBase::ft_terms amiterms;

  // Construct solution for problem defined in R0
  ami.construct(N_INT, R0, amiterms);

    //timing info 
    auto t2=std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff1=t2-t1;
    std::chrono::nanoseconds d1=std::chrono::duration_cast<std::chrono::nanoseconds>(diff1);

  std::cout<<"Construction took "<<d1.count()<<" nanoseconds"<<std::endl;

    //timing info 
    auto t3=std::chrono::high_resolution_clock::now();

  //Evaluate term-by-term solution 
  at::Tensor term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'. 
    
    //timing info
    auto t4=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff2=t4-t3;
    std::chrono::nanoseconds d2=std::chrono::duration_cast<std::chrono::nanoseconds>(diff2);

  std::cout<<"Term result was "<< format_r2_tensor(term_val)<<std::endl;
  std::cout<<"Evaluation took "<<d2.count()<<" nanoseconds"<<std::endl;
  std::cout<<"Length of terms object was " << amiterms.size() << std::endl;

  }


 void example9(){

// Number of parameters to evaluate - number of frequencies and num of energies to eval at a time
int fbatch_size = 1;
int ebatch_size = 100;

std::cout<<std::endl<<"-_-_-_ Example - Ninth Order _-_-_-"<<std::endl<<std::endl;	
	
//START Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
// class instance
at::Device myDev = at::kCPU;
TamiBase ami(myDev);

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_example_J(); // Sets initial integrand 
TamiBase::ami_vars avars=construct_ext_example_J(ami, ebatch_size, fbatch_size); // Sets 'external' parameter values 

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=9;  // Number of Matsubara sums to perform
TamiBase::ami_parms test_amiparms(N_INT, E_REG);

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

//simplified storage type 
TamiBase::ft_terms amiterms;

// Construct solution for problem defined in R0
ami.construct(N_INT, R0, amiterms);

	//timing info 
	auto t2=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff1=t2-t1;
	std::chrono::nanoseconds d1=std::chrono::duration_cast<std::chrono::nanoseconds>(diff1);

std::cout<<"Construction took "<<d1.count()<<" nanoseconds"<<std::endl;

	//timing info 
	auto t3=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
at::Tensor term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'. 
	
	//timing info
	auto t4=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::chrono::nanoseconds d2=std::chrono::duration_cast<std::chrono::nanoseconds>(diff2);

std::cout<<"Term result was "<< format_r2_tensor(term_val)<<std::endl;
std::cout<<"Evaluation took "<<d2.count()<<" nanoseconds"<<std::endl;
std::cout<<"Length of terms object was " << amiterms.size() << std::endl;
}
