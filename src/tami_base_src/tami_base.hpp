#pragma once

#include <vector>
#include <complex>
#include <stdlib.h>
#include <math.h>
#include <chrono>

#include <boost/graph/adjacency_list.hpp>

#include <torch/extension.h>
#include <torch/python.h>
#include <torch/torch.h>

class TamiBase {

public:
  /// Typedef for the type of complex numbers to be used throughout - c10 has
  /// torch's complex scalar class c10::complex<T> is the prefered complex
  /// datatype to build tensors with and therefore the better option to use
  /// throughout
  typedef c10::complex<double> complex_double;

  // Each instance will hold a torch device so that all internal calculartions
  // are performed on this device
  // at::Device device = at::kCPU; // default behaviour is to run on cpu - see
  // new c'tor

  /// Hold the device options inside the base class - specifically the
  /// options.dtype and options.device
  at::TensorOptions options =
      at::TensorOptions().dtype(at::kComplexDouble).device(at::kCPU);

  at::Device getDevice() { return options.device(); }

  // Maximum argument allowed for fermi and bose functions - to prevent inf due
  // to double precision numbers
  // Note: Possible that this prevents catching overflow issues
  double exp_max_arg = 500.0;

  bool drop_bosonic_diverge =
      false; // if set to true, then E_REG not needed because bosonic
             // divergences will simply be set to zero.  Rigorously this may
             // not be correct.
  bool drop_der = false; // if set to true, then all fermi/bose derivative terms
                         // will be dropped
  bool drop_matsubara_poles =
      true; // if set to true, ignores Matsubara poles with zero energy
  // bool is_real_external=false;
  bool zero_external_w = false;
  bool overflow_detected = false;
  bool verbose = false; // flag for verbose output - primarily for debugging

  double precision_cutoff =
      1e15; // By default this is set to roughly the precision of double.  If
            // values exceed this then numerical overflow is virtually
            // guaranteed

  bool bosonic_external = 0;

  // External list of energies and frequencies
  /// The energy of each denominator will always appear as a linear combination
  /// of these initial (pre integration) energies, \f$\epsilon_1, \epsilon_2\f$
  /// ..etc By convention, the energy_t contains the NEGATIVE of the energy of a
  /// given Green's function line, \f$ 1/(X+E) \f$ where \f$ E=-\epsilon \f$.
  ///  The `energy_t` is a tensor of sets of appropriate energy lists, allowing
  ///  for batch evaluation.  This is particularly useful for lattice problems
  ///  where one wants to sum over sets of momenta which convert to energies.
  typedef at::Tensor
      energy_t; // TODO: is there a way to declare this so it contains strictly
                // c10::complex<double> = complex_double

  /// This is the list of internal and external frequencies values (each column)
  /// and batches of such values to be evalated simultaneously (each row).
  /// Typically only the last elements for external frequencies are non-zero -
  /// but one can evaluate intermediate steps where multiple external
  /// frequencies are non-zero. (Warning: the fermi/bose statistics are not
  /// generalized to multiple internal Matsubara frequencies.  If this
  /// functionality is needed please contact the authors.)
  ///
  typedef at::Tensor frequency_t;

  // Fundamental objects

  // the symbolic epsilon

  /// Vector of type `int` with elements \f$ a_i\f$.  This is the symbolic
  /// representation of the energy \f$E=-\sum\limits_{i}a_i\epsilon_i\f$
  /// described in AMI paper (https://doi.org/10.1103/PhysRevB.99.035120).  We
  /// use the convention \f$G=\frac{1}{X+E}\f$.  It is the coefficients for a
  /// linear combination of a set of possible values.
  typedef std::vector<int>
      epsilon_t; // ASK JAMES: why not but these on GPU to use matrix
                 // multiplication in eval_prod as well + avoid loops

  /// Vector of type `int` with elements \f$ \alpha_i\f$.  This is the symbolic
  /// representation of the frequency, as a linear combination of possible
  /// entries.  Typically contains only values of 0, -1 and +1. Other values at
  /// intermediate steps typically represent an error.  \f$X=\sum\limits_{i}
  /// i\nu_i \alpha_i\f$.
  typedef std::vector<int>
      alpha_t; // ASK JAMES: why not but these on GPU to use matrix
               // multiplication in eval_prod as well + avoid loops

  /// Indicator for multi-species Green's function or energy dispersions
  /// (spin-up vs spin-dn, multiband, etc).  This indicator should propagate
  /// through the Matsubara sums to the final expression, and might have utility
  /// for evaluating energies.  Technically this is not necessary,
  /// but may be useful.
  typedef int species_t;

  /// Indicator for statistics. A future version might use this more frequently.
  /// Current version presumes all integration frequencies are Fermionic.
  typedef enum { Bose, Fermi } stat_type;

  // Eventually these types will not appear in the ami_base class
  /// Graph types will likely be removed/replaced in a future release.  Current
  /// support is limited to Sigma and Pi_phuu graph types.  set graph_type=0 for
  /// Fermionic external line, and =1 for Bosonic.
  typedef enum {
    Sigma,
    Pi_phuu,
    Pi_phud,
    Hartree,
    Bare,
    Greens,
    density,
    doubleocc,
    Pi_ppuu,
    Pi_ppud,
    DOS,
    ENERGY,
    FORCE
  } graph_type;

  /// To be removed in a future release
  typedef enum { hubbard, coulomb } int_type;
  /// To be removed in a future release
  typedef enum { tb, fp, hf } disp_type;
  /// To be removed in a future release
  typedef enum { matsubara, real } ext_type;

  /// currently not used - to be removed
  ext_type ext_freq_type = matsubara;

  /// The `ami_vars` struct is the basic information required for the evaluation
  /// stage of AMI result.  It contains the variable internal/external
  /// quantities.  Specifically it is a list of numerical values for energies of
  /// each line and values for each frequency.  Also stored is the possibility
  /// of an overall prefactor. Also required is a value of \f$\beta=\frac{1}{k_B
  /// T}\f$ needed for evaluation of Fermi/Bose distributions.
  struct ami_vars {
    /// Numerical values of energies. - at::Tensor = {e1_vector, e2_vector, ...
    /// } columns of length batch_size to be evaluated. Allows multiple energies
    ///  to be evaluated simultaneously to then numerically evaluate spatial integrals.
    energy_t energy_;
    /// Numerical Values of frequencies stored in a at::Tensor. Follows convention that
    /// external frequency is stored in last element of rows. Then frequency_t object is
    /// a stack of objects of the form \f$(0, 0, \cdot\cdot\cdot, \nu_{ext})\f$.
    frequency_t frequency_;
    /// Overall prefactor - default(1).
    double prefactor = 1.0;
    /// Required value of inverse temperature, \f$\beta\f$.
    double BETA_ = 0.0;

    /// Experimental parameter for spectral representation.
    double gamma_ = 0;

    ami_vars(energy_t eps, frequency_t freq) {
      energy_ = eps;
      frequency_ = freq;
      prefactor = 1.0;
    }

    ami_vars(energy_t eps, frequency_t freq, double Bta) {
      energy_ = eps;
      frequency_ = freq;
      BETA_ = Bta;
    }

    ami_vars(energy_t eps, frequency_t freq, double Bta, double pf) {
      energy_ = eps;
      frequency_ = freq;
      prefactor = pf;
      BETA_ = Bta;
    }

    ami_vars() { prefactor = 1.0; }
  };

  /// Parameters for AMI construction/evaluation.
  struct ami_parms {
    /// Number of integrations to perform.
    int N_INT_;
    /// Hardcoded as (1) in this version - represents number of external
    /// variables.
    int N_EXT_ = 1;
    /// Possible energy regulator for evaluation of Fermi/Bose functions to
    /// control divergences.  Should be zero by default and switched on if
    /// needed.
    double E_REG_ = 0;

    /// Tolerence to determine equality
    double tol_ = 1e-12;

    ami_parms(int N_INT, double E_REG) {
      N_INT_ = N_INT;
      E_REG_ = E_REG;
      N_EXT_ = 1;
      TYPE_ = static_cast<TamiBase::graph_type>(0);   /// by default sigma.
      int_type_ = static_cast<TamiBase::int_type>(0); /// by default is hubbard.
      dispersion_ =
          static_cast<TamiBase::disp_type>(0); /// by default is tight-binding.
    }

    ami_parms(int N_INT, double E_REG, graph_type TYPE) {
      N_INT_ = N_INT;
      E_REG_ = E_REG;
      TYPE_ = TYPE;
      int_type_ = static_cast<TamiBase::int_type>(0);
      N_EXT_ = 1;
    }

    ami_parms(int N_INT, double E_REG, graph_type TYPE, int_type inter,
              disp_type disp) {
      N_INT_ = N_INT;
      E_REG_ = E_REG;
      TYPE_ = TYPE;
      N_EXT_ = 1;
      int_type_ = static_cast<TamiBase::int_type>(inter);
      dispersion_ = static_cast<TamiBase::disp_type>(disp);
    }

    ami_parms() {}

    graph_type TYPE_;
    int_type int_type_;
    disp_type dispersion_;
  };

  /// A Green's function structure.  This is a symbolic vector of `epsilon_t`
  /// and vector of `alpha_t`.  Also needed is what the statistics of the line
  /// are.  While this could be determined from alpha - it is better to store
  /// it.  For multistate systems the species_ index might be useful.
  struct g_struct {
    /// Constructor with state_type specified.
    g_struct(epsilon_t eps, alpha_t alpha, stat_type stat) {
      eps_ = eps;
      alpha_ = alpha;
      stat_ = stat;
      species_ = 0;
      eff_stat_ = stat_;
    }

    /// Constructor assumes fermi statistics if not specified for partially
    /// initialized structure.
    g_struct(epsilon_t eps, alpha_t alpha) {
      eps_ = eps;
      alpha_ = alpha;
      stat_ = Fermi;
      species_ = 0;
      eff_stat_ = stat_;
    }

    /// Uninitialized variant.
    g_struct() {
      stat_ = Fermi;
      species_ = 0;
      eff_stat_ = stat_;
    }

    /// Symbolic coefficients for a linear combination of epsilons.
    epsilon_t eps_;

    /// Symbolic coefficients for a linear combination of frequencies.
    alpha_t alpha_;
    /// Mark Fermi/Bose stats - Fermi is required in current version.
    stat_type stat_;
    species_t species_;

    stat_type eff_stat_;
    // Experimental Spectral representation.  Implementation incomplete.
    int pp = -1; // pp=0 means this G represents a principle part integral.
                 // pp=1 it is a delta function. else it is inert
  };

  /// Pole structure. Nearly equivalent to `g_struct`, but given its own struct
  /// for clarity. Tracks multiplicity, and which Green's function it is
  /// attached to. Also it tracks how many derivatives to take when evaluated at
  /// a fermi function.
  struct pole_struct {
    pole_struct(epsilon_t eps, alpha_t alpha) {
      eps_ = eps;
      alpha_ = alpha;
    }

    pole_struct() {}

    epsilon_t eps_;
    alpha_t alpha_;
    /// Index that specifies which frequency it is a pole with respect to.
    int index_;
    /// The multiplicity of the pole, starts at 1 and increments as needed.
    int multiplicity_ = 1;
    int der_ = 0;              /**< Counter for derivatives. */
    std::vector<int> which_g_; /**< Index to identify which `g_struct` a pole
                                originated from.*/

    /// Experimental component of Spectral evaluation.
    alpha_t x_alpha_;
  };

  // MB: Maybe need more of the types -- im pretty sure this is all that is
  // needed for terms
  typedef std::vector<g_struct> g_prod_t;
  typedef std::vector<pole_struct> pole_array_t;
  typedef std::vector<g_prod_t> Ri_t;

  /** Term Structure for term-by-term evaluation.  Conceptually simpler than SPR
   * construction. Storage translates to \f$ \prod{f(p_i)}\prod{G_j}\times sign
   * \f$.
   *
   */
  struct term {
    term() {}

    term(double s, pole_array_t p, g_prod_t g) {
      sign = s;
      p_list = p;
      g_list = g;
    }

    /// Sign prefactor
    double sign = 1;
    /// List of poles, \f$ \prod{f(p_i)}\f$.
    pole_array_t p_list;
    /// List of Green's functions, \f$ \prod{G_j}\f$.
    g_prod_t g_list;
  };

  /// The storage for the term-by-term construction.  Each `term` struct is an
  /// element of the `terms` vector.
  typedef std::vector<term> terms;

  typedef std::pair<int, int> ref_t;
  typedef std::vector<ref_t> ref_v_t;
  typedef std::vector<ref_v_t> R_ref_t;
  typedef R_ref_t ref_eval_t;

  // START OF THE POLE TREE STUFF

  /// @brief  The FermiTree class is a graph representation of sums and products
  /// of Fermi and Bose functions. It allows for a more natural factorization of
  /// the AMI integrand after each integration step.  It provides additional
  /// numerical stability compared to the fully expanded product of Fermi/Bose
  /// functions as used in the `libami` implementation. It is not expected that
  /// a user should interact with the FermiTree, since it is not a pretty
  /// structure to look at and the detailed structure will depend on the
  /// labelling of a diagram.
  class FermiTree {
  public:
    enum operation { add, mult, end };

    // Vertex info

    struct vertex_info {

      vertex_info() {
        visited = 0;
        depth_ = -1;     // default value so easier to catch.
        operation_ = -1; // set to nonsense value for catching bugs.
        prefactor_ = 1;
      }

      vertex_info(int op, double depth) {
        operation_ = op;
        depth_ = depth;
        prefactor_ = 1;
      }

      vertex_info(TamiBase::pole_struct pole, int op) {
        operation_ = op;
        pole_ = pole;
        prefactor_ = 1;
      }

      vertex_info(TamiBase::pole_struct pole, double prefactor, int op) {
        operation_ = op;
        pole_ = pole;
        prefactor_ = prefactor;
      }

      vertex_info(TamiBase::pole_struct pole, at::Tensor value, int op,
                  double depth) {
        operation_ = op;
        depth_ = depth;
        pole_ = pole;
        value_ = value;
        prefactor_ = 1;
      }

      vertex_info(TamiBase::pole_struct pole, at::Tensor value, int op,
                  double depth, double prefactor) {
        operation_ = op;
        depth_ = depth;
        pole_ = pole;
        value_ = value;
        prefactor_ = prefactor;
      }

      at::Tensor value_;
      int operation_;

      int visited;
      boost::property<boost::vertex_index_t, int> vertex_index;
      int index_;
      int depth_;

      TamiBase::pole_struct pole_;
      double prefactor_;
    };

    // Graph global info

    struct graph_info {

      int max_depth_ = 0;
      // double overall_prefactor=1.0;
    };

    /////////////////////////

    ////////////////////////
    // Edge info structure
    ////////////////////////

    struct edge_info {

      int depth_; // say this is the depth of either the target or the source

      edge_info() { depth_ = -1; }

      // todo fix this - need to remove loop id
      edge_info(int depth) { depth_ = depth; }

    }; // end edge_struct bracket

    /// @brief Boost Adjacency list used to manage fermi tree structure
    typedef boost::adjacency_list<boost::vecS, boost::listS,
                                  boost::bidirectionalS, vertex_info, edge_info,
                                  graph_info>
        fermi_tree_t;
    typedef boost::graph_traits<fermi_tree_t>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<fermi_tree_t>::edge_descriptor edge_t;

    

    // Initialize the tree to have one
    void initialize_ft(fermi_tree_t &ft);
    void initialize_ft(fermi_tree_t &ft, operation op);
    void initialize_ft(fermi_tree_t &ft, TamiBase::pole_struct &pole);
    void initialize_ft(fermi_tree_t &ft, TamiBase::pole_struct &pole,
                       double &prefactor);

    void plist_to_ft(TamiBase::pole_array_t &plist,
                     TamiBase::FermiTree::fermi_tree_t &ft);
    void plist_to_ft(TamiBase::pole_array_t &plist, double sign,
                     TamiBase::FermiTree::fermi_tree_t &ft);

    // this changes every prefactor and is not in general the same as an overall
    // prefactor.
    void update_prefactors(fermi_tree_t &ft, double sign);
    vertex_t get_root(fermi_tree_t &ft);
    bool is_empty_ft(fermi_tree_t &ft);

    void mult_prefactor(fermi_tree_t &ft, double sign);

    void get_roots(fermi_tree_t &ft, std::vector<vertex_t> &vv);

    // TODO: Do I want it to destroy the level info on call?
    void get_next_level(fermi_tree_t &ft1, vertex_t &root1,
                        std::vector<vertex_t> &next_level);

    void print_vertex(vertex_t &v, fermi_tree_t &ft);

    void print_graph(fermi_tree_t &ft);

    std::string pretty_print(fermi_tree_t &ft);
    
    // could either place
    fermi_tree_t add_ft(fermi_tree_t ft1, fermi_tree_t ft2);

    fermi_tree_t mult_ft(fermi_tree_t ft1, fermi_tree_t ft2);

    FermiTree();

  private:

    // helper functions
    void number_vertices(fermi_tree_t &g);

    // internal print functions
    std::string pretty_print_pole(TamiBase::pole_struct &pole);
    std::string pretty_print(fermi_tree_t &ft, vertex_t &v);
    std::string pretty_print_ft(fermi_tree_t &ft1, vertex_t &v);
    std::string pretty_print(fermi_tree_t &ft, std::vector<vertex_t> &vv,
                             int op);
    std::string pretty_print_level(fermi_tree_t &ft, std::vector<vertex_t> &vv,
                                   int source_op);

    // this should copy root1 from ft1 to ft2 below root2
    // will fill the pair with the vertex_t of the new vertex on ft1 and ft2
    void copy_vertex(vertex_t &root1, vertex_t &root2, fermi_tree_t &ft1,
                     fermi_tree_t &ft2, std::pair<vertex_t, vertex_t> &map);

    // copy the level below root1 to below root2
    // will populate a map relating the next level of root1 to the new vertices
    // created on ft2 below root2
    void copy_level(vertex_t &root1, vertex_t &root2, fermi_tree_t &ft1,
                    fermi_tree_t &ft2,
                    std::vector<std::pair<vertex_t, vertex_t>> &map_vec);

    // copy entire tree from ft1 starting at root1 below root2
    //  either include root1 or don't.  which is it?
    void copy_tree(fermi_tree_t &ft1, fermi_tree_t &ft2, vertex_t &root1,
                   vertex_t &root2);



  };

  /// Both the initial integrand and intermediate steps can be  written as a
  /// fermi-tree_term `ft_term`.  Non-factorizable terms will spawn additional
  /// terms. Includes a prefactor, a `fermi_tree_t` and a `g_prod_t`.
  struct ft_term {
    ft_term() {}

    ft_term(double s, TamiBase::FermiTree::fermi_tree_t ft,
            TamiBase::g_prod_t g) {
      sign_ = s;
      ft_ = ft;
      g_prod_ = g;
    }

    /// Sign prefactor
    double sign_ = 1;
    /// List of poles, \f$ \prod{f(p_i)}\f$.
    TamiBase::FermiTree::fermi_tree_t ft_;
    /// List of Green's functions, \f$ \prod{G_j}\f$.
    TamiBase::g_prod_t g_prod_;
  };

  /// The result of an integrand is this vector of `ft_term` structs. This is
  /// the main output of the `construct` function and the main object to pass to
  /// the `evaluate` function. See examples.
  typedef std::vector<ft_term> ft_terms;
  typedef std::vector<TamiBase::FermiTree::fermi_tree_t> ft_list;

  /// @brief Following the notation from `libami` the construct function
  /// performs a number of Matsubara sums, and places the result into a
  /// `ft_terms` type.
  /// @param N_INT - the number of summations to perform (beginning from the
  /// zero index frequency and going to N_INT-1)
  /// @param R0  - The starting integrand.  R0 is a product of green's functions
  /// @param terms_out  -  The resulting analytic integrand.
  void construct(int N_INT, TamiBase::g_prod_t &R0, ft_terms &terms_out);

  // TODO: not sure that the factorize=false survived - might need to
  // reimplement that.
  /// @brief Internal function to the construction. After each step attempts to
  /// combine terms into factorized form.
  ///  Can be slow for high order integrands, and can be disabled by setting
  ///  factorize=false.
  /// @param in_terms - terms to be factorized.
  /// @param out_terms - factorized terms.
  void factorize(ft_terms &in_terms, ft_terms &out_terms);

  void terms_to_ftterms(TamiBase::terms &in_terms, ft_terms &out_terms);

  void plist_to_ft(TamiBase::pole_array_t &plist,
                   TamiBase::FermiTree::fermi_tree_t &ft);

  void take_term_derivative(ft_term &in_term, TamiBase::pole_struct &pole,
                            ft_terms &out_terms);

  /// Screen IO for debugging.
  void print_terms(ft_terms &t);
  /// Screen IO for debugging.
  void print_term(ft_term &t);

  /// Experimental attempt for visualizing the post AMI integrand. Outputs
  /// string representing LaTeX formatted equation.
  std::string pretty_print_ft_term(ft_term &ft);
  std::string pretty_print_ft_terms(ft_terms &fts);

  /// @brief Evaluate just the fermi-tree
  at::Tensor eval_ft(TamiBase::ami_parms &parms,
                     TamiBase::FermiTree::fermi_tree_t &ft1,
                     TamiBase::FermiTree::vertex_t &v,
                     TamiBase::ami_vars &external);

  /// Evaluation of single term.
  at::Tensor evaluate_term(TamiBase::ami_parms &parms, ft_term &ft_term,
                           TamiBase::ami_vars &external);

  /// @brief Standard evaluate function. See examples.
  /// @return at::Tensor of shape \f$N_f \times N_e\f$, so that in practice
  /// one can sum the columns in a Monte Carlo integration
  at::Tensor evaluate(TamiBase::ami_parms &parms, ft_terms &ft_terms,
                      TamiBase::ami_vars &external);

  FermiTree FT;

  // C'tors
  /// Default Constructor.  Constructor is empty.  Currently no initialization
  /// is required in most cases.
  TamiBase() {}
  // torch device c'tor otherwise it intializes to at::kCPU
  TamiBase(at::Device &dev) {
    options = at::TensorOptions().dtype(at::kComplexDouble).device(dev);
  }

  // FUNCTIONS THAT HAVE BEEN MODIFIED FROM LIBAMI

  /*Given an array of Green's functions, finds all poles with respect to
  frequency index, and checks for multiplicities. Stores multiplicities of the
  poles as well as `which_g_`, an identifier that specifies which Green's
  function it was attached to.
  */
  pole_array_t find_poles(int index, g_prod_t &R);

  // Derivative for term construction
  void take_term_derivative(term &in_term, pole_struct &pole, terms &out_terms);

  // Functions for fermi_pole

  /// This is the central evaluation of the fermi and bose functions.  It also
  /// includes evaluating arbitrary derivatives of the functions.  See frk
  /// function that is rather complicated .  This function is also the MOST
  /// challenging function for numerical evaluation.
  at::Tensor fermi_pole(ami_parms &parms, pole_struct pole, ami_vars external);

  /// Recursive construction of fermi_bose derivatives.
  // non-private because it is used in tests
  at::Tensor fermi_bose(int m, double sigma, double beta, at::Tensor E);
  // Extensions for overflow safe primitives
  at::Tensor fermi_bose_mlte2(int m, double sigma, double beta, at::Tensor E);
  at::Tensor fermi_bose_agrad(int m, double sigma, double beta, at::Tensor E); 

  /// Given a set of tensor energies, beta, and tensor frequencies, will
  /// evaluate the energies of a pole_struct.
  at::Tensor get_energy_from_pole(pole_struct pole, ami_vars external);

  // Residue function
  g_struct update_G_pole(g_struct g_in, pole_struct pole);

  /// @brief Primarily debugging function
  /// @param `pole_struct`
  void print_pole_struct_info(pole_struct g);

  /// @brief Primarily debugging function
  /// @param `alpha_t`
  void print_alpha_info(alpha_t alpha);

  /// @brief Primarily debugging function
  /// @param `epsilon_t`
  void print_epsilon_info(epsilon_t eps);

  // Functions for Factorize_Rn

  /// Optimize function for old SPR notation.  Possibly Depricated. 
  void factorize_Rn(Ri_t &Rn, g_prod_t &unique_g, R_ref_t &Rref,
                    ref_eval_t &Eval_list);

  bool g_struct_equiv(g_struct &g1, g_struct &g2, int &sign);

  void reduce_rref(R_ref_t &Rref, ref_eval_t &Eval_list);

  bool pair_v_equiv(ref_v_t &r1, ref_v_t &r2, int &r1sign, int &r2sign);

  bool pole_equiv(pole_struct pole1, pole_struct pole2);

private:

  bool g_prod_equiv(TamiBase::g_prod_t &gp1, TamiBase::g_prod_t &gp2,
                    int &sign);

  /// @brief Integrates a single Matsubara index for every term provided using
  /// the `term_integrate_step` function.
  void integrate_step(int index, ft_terms &in_terms, ft_terms &out_terms);
  /// @brief Integrates a single step for a single term
  void term_integrate_step(int index, ft_term &in_term, ft_terms &out_terms);
  void terms_general_residue(ft_term &this_term,
                             TamiBase::pole_struct this_pole,
                             ft_terms &out_terms);

  std::string pretty_print_gprod(TamiBase::g_prod_t &gp);
  std::string pretty_print_g(TamiBase::g_struct &g);

  // LIBAMI converted functions that don't need user access
  
  /// Evaluates a product of Green's functions for all energies provided
  at::Tensor eval_gprod(ami_parms &parms, g_prod_t g_prod, ami_vars external);

  /// Convert terms to Ri structure for optimization.  This does not create a
  /// usable `terms` object. It is purely an intermediate step for the
  /// optimization functions. Note: this is largely superceded by the
  /// factorization process.
  void convert_terms_to_ri(terms &ami_terms, Ri_t &Ri);

  // Functions for amiBase::integrate_step -> integrate_Mat_ind_step

  /// Integrates a single Matsubara index for the non-fermi-tree `terms` object
  /// from previous implementation. Using this intermediate representation is
  /// useful since it is directly from `libami` it is well tested.  (renamed
  /// from the original function name (TamiBase::integrate_step()) The output
  /// terms are then factorized into a fermi tree.
  void integrate_Mat_ind_step(int index, terms &in_terms, terms &out_terms);


  // Residue Functions
  g_prod_t simple_residue(g_prod_t G_in, pole_struct pole);

  // sign function
  double get_simple_sign(int index, g_prod_t &R, pole_struct pole);

  /// Primary residue function for term construction.
  void terms_general_residue(term &this_term, pole_struct this_pole,
                             terms &out_terms);



  double get_starting_sign(g_prod_t G_in, pole_struct pole);

  // This function removes the inert parts of the gprod in the context of taking
  // derivatives
  g_prod_t reduce_gprod(g_prod_t G_in, pole_struct pole);


};

namespace mathUtils{

  /// Returns the sign of a value - or zero if it is uniquely zero.
  template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

  // template for checking equality of vectors
  template <typename T>
  bool isEqual(std::vector<T> const &v1, std::vector<T> const &v2) {
    return (v1.size() == v2.size() &&
            std::equal(v1.begin(), v1.end(), v2.begin()));
  }

  /// Simple factorial function. It is rarely called except for multipole
  /// problems, and typically for only small arguments.
  int factorial(int n);

  /*The frk function itself returns
  \f[
  f_{rk}(m,k)=\sum\limits_{m=0}^{k+1} binomialCoeff(k,m) m^r (-1)^{k-m}.
  \f]
  */
  double frk(int r, int k);

  /*
  Recursive binomial coefficient function
  */
  int binomialCoeff(int n, int k);
}


std::ostream & operator<<(std::ostream &os, TamiBase::ft_term &ft);
std::ostream & operator<<(std::ostream &os, TamiBase::ft_terms &ft);

//  std::ostream & operator<<(std::ostream &os, TamiBase::ft_term &ft){

//   TamiBase ami;
//       return os<<ami.pretty_print_ft_term(ft);
//     }
