#include "tami_base.hpp"
#include "examples.hpp"


int main( int argc , char *argv[] )
{
  int mode=0;
  if (argc >= 2){
    std::istringstream ss(argv[1]);

    if (!(ss >> mode)){
      std::cerr << "Using default" << std::endl;
      mode = 99;
    }
  }

  switch(mode) {
    case 0:
      default_example();
      break;
    case 1:
      example_1();
      break;
    case 2:
      example2();
      break;
    case 3:
      example1_bose();
      break;
    case 4:
      example4();
      break;
    case 5:
      example9();
      break;
    default:
      example2();
      example1_bose();
      example4();
      example9();
      break;
  }

  return 0;
}



void default_example(){

  // Terms example

  TamiBase::g_prod_t R0=construct_example2();

  TamiBase PT;
  TamiBase::ft_terms ftout;

  /* PT.construct(2, R0, ftout);

  FermiTree::vertex_t r=PT.FT.get_root(ftout[0].ft_);
  std::cout<<PT.FT.pretty_print_ft(ftout[0].ft_,r)<<std::endl;

  PT.FT.number_vertices(ftout[0].ft_);
  PT.FT.print_graph(ftout[0].ft_);

  std::cout<<PT.pretty_print_ft_terms(ftout)<<std::endl;
  TamiBase::ami_vars ext=construct_ext_example2();

  double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
  int N_INT=2;  // Number of Matsubara sums to perform
  TamiBase::ami_parms parms(N_INT, E_REG);

  std::complex<double> result=PT.evaluate(parms,ftout,ext);

  std::cout<<"Result is "<<result<<std::endl; */

  TamiBase::g_prod_t R02=construct_multipole_example();//construct_example_J();//construct_multipole_example();//construct_example1_bose();
  TamiBase::ami_vars avars2=construct_4ord_ext_multipole_example();//construct_ext_example_J();//construct_4ord_ext_multipole_example();//construct_ext_example1_bose();

  // construct_4ord_ext_multipole_example();
  // TamiBase::g_prod_t construct_multipole_example();

  TamiBase::ft_terms ftout2;

  PT.construct(4, R02, ftout2);
  TamiBase::ami_parms parms2(0, 0);

  auto t2=std::chrono::high_resolution_clock::now();
  std::complex<double> result2=PT.evaluate(parms2,ftout2,avars2);

  auto t_end=std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> diff2=t_end-t2;

  // std::chrono::microseconds d=std::chrono::duration_cast<std::chrono::microseconds>(diff1);



  std::cout<<"--------------------"<<std::endl;
  // PT.FT.number_vertices(ftout2[0].ft_);
  // PT.FT.print_graph(ftout2[0].ft_);
  // std::cout<<PT.FT.pretty_print(ftout2[0].ft_)<<std::endl;

  std::cout<<PT.pretty_print_ft_terms(ftout2)<<std::endl;
  std::cout<<"Result 4th order MP "<<result2<<std::endl;
  std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);
  std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	

  // PtTamiBase::ft_terms ftout3;
  // ftout3.push_back(ftout2[1]);
  // ftout3.push_back(ftout2[6]);

  // PT.FT.number_vertices(ftout3[0].ft_);
  // PT.FT.number_vertices(ftout3[1].ft_);
  // ftout3.push_back(ftout2[3]);

  // std::complex<double> result3=PT.evaluate(parms2,ftout3,avars2);
  // std::cout<<"First term gives "<< result3<<std::endl;
  // std::cout<<PT.pretty_print_ft_terms(ftout3)<<std::endl;

  // PT.FT.print_graph(ftout3[1].ft_);
  // PT.FT.mult_prefactor(ftout3[0].ft_,-1);
  // PT.FT.print_graph(ftout3[1].ft_);
  // result3=PT.evaluate(parms2,ftout3,avars2);
  // std::cout<<"After mult First term gives "<< result3<<std::endl;
  // std::cout<<PT.pretty_print_ft_terms(ftout3)<<std::endl;

  // std::cout<<"Factorize!----------"<<std::endl;

  // PtTamiBase::ft_terms factorized;
  // PT.factorize(ftout3,factorized);
  // std::cout<<PT.pretty_print_ft_terms(factorized)<<std::endl;
  // std::complex<double> result4=PT.evaluate(parms2,factorized,avars2);
  // std::cout<<"Factorized two terms give "<< result4<<std::endl;

  // PT.FT.number_vertices(ftout3[1].ft_);
  // PT.FT.print_graph(ftout3[1].ft_);

  // result4=PT.evaluate(parms2,ftout3,avars2);

  // PT.FT.print_graph(ftout3[1].ft_);
  // result4=PT.evaluate(parms2,ftout3,avars2);
  // std::cout<<"NonFactorized two terms give "<< result4<<std::endl;
  // PT.FT.number_vertices(ftout2[0].ft_);
  // PT.FT.print_graph(ftout2[0].ft_);

  // TamiBase::epsilon_t eps={0,1};
  // TamiBase::alpha_t alpha={1,0};

  // TamiBase::pole_struct p(eps,alpha);

  // FermiTree FT;

  // FermiTree::fermi_tree_t ft;
  // FermiTree::vertex_t root_ft=FT.get_root(ft);
  // std::cout<<ft[root_ft].index_<<" op: "<< ft[root_ft].operation_<<std::endl;

  // FT.initialize_ft(ft);

  // FermiTree::fermi_tree_t ft2;
  // FT.initialize_ft(ft2);



  // std::cout<< num_vertices(ft)<<std::endl;
  // std::cout<< num_vertices(ft2)<<std::endl;

  // root_ft=FT.get_root(ft);
  // std::cout<<ft[root_ft].index_<<" op: "<< ft[root_ft].operation_<<std::endl;

  // FermiTree::vertex_t root_ft2=FT.get_root(ft2);
  // std::cout<<ft2[root_ft2].index_<<" op: "<< ft2[root_ft2].operation_<<std::endl;




  // std::pair<vertex_t,vertex_t> map;
  // copy_vertex(root_ft,root_ft2,ft,ft2,map);
  // std::cout<<"ft2 now has "<< num_vertices(ft2)<<std::endl;
  // print_ft(ft2);
  // copy_vertex(root_ft,root_ft2,ft,ft2,map);
  // std::cout<<"ft2 now has "<< num_vertices(ft2)<<std::endl;
  // print_ft(ft2);

  /* 
  FermiTree::fermi_tree_t ft3=FT.add_ft(ft,ft2);
  // fermi_tree_t ft5=ft3;



  FermiTree::fermi_tree_t ft4=FT.mult_ft(ft3,ft3);
  // number_vertices(ft4);
  // print_graph(ft4);

  FermiTree::fermi_tree_t ft5=FT.add_ft(FT.mult_ft(ft4,ft3),ft4);
  // fermi_tree_t ft5=ft4;
  FT.number_vertices(ft5);
  FT.print_graph(ft5);


  boost::graph_traits<FermiTree::fermi_tree_t>::vertex_iterator vi,vie;

  for(boost::tie(vi,vie)=vertices(ft5); vi!=vie; ++vi){

  if(ft5[*vi].operation_==2){ ft5[*vi].value_=2.0;}

  }

  FT.number_vertices(ft5);
  FermiTree::vertex_t r=FT.get_root(ft5);
  std::complex<double> result=FT.eval_ft(ft5,r);
  std::cout<<"Result="<<result<<std::endl;

  */
}

 void example_1(){

  std::cout << "This is example1" << std::endl;


 }

void example2(){
std::cout<<std::endl<<"-_-_-_ Example - Second Order _-_-_-"<<std::endl<<std::endl;	
	
//START Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
// class instance
TamiBase ami;

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_example2(); // Sets initial integrand 
TamiBase::ami_vars avars=construct_ext_example2(); // Sets 'external' parameter values 

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
	std::chrono::microseconds d1=std::chrono::duration_cast<std::chrono::microseconds>(diff1);

std::cout<<"Construction took "<<d1.count()<<" microseconds"<<std::endl;

	//timing info 
	auto t3=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'. Note that the test_amiparms is the same as the first case 
	
	//timing info
	auto t4=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"Evaluation took "<<d2.count()<<" microseconds"<<std::endl;
	
}

 void example1_bose(){

std::cout<<std::endl<<"-_-_-_ Example - First Order For Bosonic _-_-_-"<<std::endl<<std::endl;	
	
//START Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
// class instance
TamiBase ami;

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_example1_bose(); // Sets initial integrand 
TamiBase::ami_vars avars=construct_ext_example1_bose(); // Sets 'external' parameter values 

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=1;  // Number of Matsubara sums to perform
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
	std::chrono::microseconds d1=std::chrono::duration_cast<std::chrono::microseconds>(diff1);

std::cout<<"Construction took "<<d1.count()<<" microseconds"<<std::endl;

	//timing info 
	auto t3=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'. Note that the test_amiparms is the same as the first case 
	
	//timing info
	auto t4=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"Evaluation took "<<d2.count()<<" microseconds"<<std::endl;


 }


 void example4(){

std::cout<<std::endl<<"-_-_-_ Example - Fourth Order _-_-_-"<<std::endl<<std::endl;	
	
//START Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
// class instance
TamiBase ami;

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_multipole_example(); // Sets initial integrand 
TamiBase::ami_vars avars=construct_4ord_ext_multipole_example(); // Sets 'external' parameter values 

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
	std::chrono::microseconds d1=std::chrono::duration_cast<std::chrono::microseconds>(diff1);

std::cout<<"Construction took "<<d1.count()<<" microseconds"<<std::endl;

	//timing info 
	auto t3=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'.
	
	//timing info
	auto t4=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"Evaluation took "<<d2.count()<<" microseconds"<<std::endl;

 }


 void example9(){

std::cout<<std::endl<<"-_-_-_ Example - Ninth Order _-_-_-"<<std::endl<<std::endl;	
	
//START Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing TAMI term by term-----"<<std::endl;
// class instance
TamiBase ami;

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_example_J(); // Sets initial integrand 
TamiBase::ami_vars avars=construct_ext_example_J(); // Sets 'external' parameter values 

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
	std::chrono::microseconds d1=std::chrono::duration_cast<std::chrono::microseconds>(diff1);

std::cout<<"Construction took "<<d1.count()<<" microseconds"<<std::endl;

	//timing info 
	auto t3=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'. 
	
	//timing info
	auto t4=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"Evaluation took "<<d2.count()<<" microseconds"<<std::endl;

}