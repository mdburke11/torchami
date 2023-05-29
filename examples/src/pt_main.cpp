#include "tami_base.hpp"
#include "examples.hpp"


int main( int argc , char *argv[] )
{
  int mode=0;
  if (argc >= 2){
    std::istringstream ss(argv[1]);

    if (!(ss >> mode)){
      std::cerr << "Something weird happened" << std::endl;
    }
  }

  switch(mode) {
    case 0:
      default_example();
      break;
    case 1:
      example_1();
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