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

class AmiBase{

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

        // Maximum argument allowed for fermi and bose functions - to prevent inf due
        // to double precision numbers
        // Note: Possible that this prevents catching overflow issues 
        double exp_max_arg = 500.0;

        /// Simple factorial function. It is rarely called except for multipole
        /// problems, and typically for only small arguments.
        int factorial(int n);

        bool drop_bosonic_diverge =
            false; // if set to true, then E_REG not needed because bosonic
                    // divergences will simply be set to zero.  Rigorously this may
                    // not be correct.
        bool drop_der = false; // if set to true, then all fermi/bose derivative terms will be dropped
        bool drop_matsubara_poles = true; // if set to true, ignores Matsubara poles with zero energy
        // bool is_real_external=false;
        bool zero_external_w = false;
        bool overflow_detected = false;
        bool verbose = false; // flag for verbose output - primarily for debugging

        double precision_cutoff =
            1e15; // By default this is set to roughly the precision of double.  If
                    // values exceed this then numerical overflow is virtually
                    // guaranteed
                    
        bool bosonic_external=0;

        // External list of energies and frequencies
        /// The energy of each denominator will always appear as a linear combination
        /// of these initial (pre integration) energies, \f$\epsilon_1, \epsilon_2\f$
        /// ..etc By convention, the energy_t contains the NEGATIVE of the energy of a
        /// given Green's function line, \f$ 1/(X+E) \f$ where \f$ E=-\epsilon \f$.
        typedef std::vector<std::complex<double>> energy_t;

        /// This is the list of internal and external frequencies values.  Typically
        /// only the last elements for external frequencies are non-zero - but one can
        /// evaluate intermediate steps where multiple external frequencies are
        /// non-zero.
        ///
        typedef std::vector<std::complex<double>> frequency_t;

        // Fundamental objects

        // the symbolic epsilon

        /// Vector of type `int` with elements \f$ a_i\f$.  This is the symbolic
        /// representation of the energy \f$E=-\sum\limits_{i}a_i\epsilon_i\f$
        /// described in AMI paper (https://doi.org/10.1103/PhysRevB.99.035120).  We
        /// use the convention \f$G=\frac{1}{X+E}\f$.  It is the coefficients for a
        /// linear combination of a set of possible values.
        typedef std::vector<int> epsilon_t;

        /// Vector of type `int` with elements \f$ \alpha_i\f$.  This is the symbolic
        /// representation of the frequency, as a linear combination of possible
        /// entries.  Typically contains only values of 0, -1 and +1. Other values at
        /// intermediate steps typically represent an error.  \f$X=\sum\limits_{i}
        /// i\nu_i \alpha_i\f$.
        typedef std::vector<int> alpha_t;

        /// Indicator for multi-species Green's function or energy dispersions
        /// (spin-up vs spin-dn, multiband, etc).  This indicator should propagate
        /// through the Matsubara sums to the final expression, and might have utility
        /// for evaluating energies.  Technically this is not necessary for libami,
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

        ext_type ext_freq_type = matsubara;

        /// The `ami_vars` struct is the basic information required for the evaluation
        /// stage of AMI result.  It contains the variable internal/external
        /// quantities.  Specifically it is a list of numerical values for energies of
        /// each line and values for each frequency.  Also stored is the possibility
        /// of an overall prefactor. Also required is a value of \f$\beta=\frac{1}{k_B
        /// T}\f$ needed for evaluation of Fermi/Bose distributions.
        struct ami_vars {
            /// Numerical values of energies.
            energy_t energy_;
            /// Numerical Values of frequencies.
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
            TYPE_ = static_cast<AmiBase::graph_type>(0);   /// by default sigma.
            int_type_ = static_cast<AmiBase::int_type>(0); /// by default is hubbard.
            dispersion_ =
                static_cast<AmiBase::disp_type>(0); /// by default is tight-binding.
            }

            ami_parms(int N_INT, double E_REG, graph_type TYPE) {
            N_INT_ = N_INT;
            E_REG_ = E_REG;
            TYPE_ = TYPE;
            int_type_ = static_cast<AmiBase::int_type>(0);
            N_EXT_ = 1;
            }

            ami_parms(int N_INT, double E_REG, graph_type TYPE, int_type inter,
                    disp_type disp) {
            N_INT_ = N_INT;
            E_REG_ = E_REG;
            TYPE_ = TYPE;
            N_EXT_ = 1;
            int_type_ = static_cast<AmiBase::int_type>(inter);
            dispersion_ = static_cast<AmiBase::disp_type>(disp);
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
            eff_stat_=stat_;
            }

            /// Constructor assumes fermi statistics if not specified for partially
            /// initialized structure.
            g_struct(epsilon_t eps, alpha_t alpha) {
            eps_ = eps;
            alpha_ = alpha;
            stat_ = Fermi;
            species_ = 0;
            eff_stat_=stat_;
            }

            /// Uninitialized variant.
            g_struct() {
            stat_ = Fermi;
            species_ = 0;
            eff_stat_=stat_;
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

        /// Pole structure. Equivalent to `g_struct`, but kept separate. Tracks
        /// multiplicity, and which Green's function it is attached to. Also it tracks
        /// how many derivatives to take when evaluated at a fermi function.
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

        // MB: Maybe need more of the types -- im pretty sure this is all that is needed for terms
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

        FermiTree FT;

        // FUNCTIONS THAT HAVE BEEN PULLED FROM LIBAMI

        /// Evaluates a product of Green's functions.
        std::complex<double> eval_gprod(ami_parms &parms, g_prod_t g_prod,
                                  ami_vars external);

        /// Convert terms to Ri structure for optimization.  This does not create a
        /// usable `terms` object. It is purely an intermediate step for the
        /// optimization functions.
        void convert_terms_to_ri(terms &ami_terms, Ri_t &Ri);

        // Functions for amiBase::integrate_step -> integrate_Mat_ind_step

        /// Integrates a single Matsubara index. --  renamed from the original function name (AmiBase::integrate_step())
        void integrate_Mat_ind_step(int index, terms &in_terms, terms &out_terms);

        /*Given an array of Green's functions, finds all poles with respect to frequency
        index, and checks for multiplicities. Stores multiplicities of the poles as
        well as `which_g_`, an identifier that specifies which Green's function it was
        attached to.
        */
        pole_array_t find_poles(int index, g_prod_t &R);

        // Residue Functions
        g_prod_t simple_residue(g_prod_t G_in, pole_struct pole);

        // sign function
        double get_simple_sign(int index, g_prod_t &R, pole_struct pole);

        /// Primary residue function for term construction.
        void terms_general_residue(term &this_term, pole_struct this_pole,
                             terms &out_terms);

        // Residue function
        g_struct update_G_pole(g_struct g_in, pole_struct pole);

        double get_starting_sign(g_prod_t G_in, pole_struct pole);

        // Derivative for term construction
        void take_term_derivative(term &in_term, pole_struct &pole, terms &out_terms);

        /// Simple factorial function. It is rarely called except for multipole
        /// problems, and typically for only small arguments.
        int factorial(int n);

        // This function removes the inert parts of the gprod in the context of taking
        // derivatives
        g_prod_t reduce_gprod(g_prod_t G_in, pole_struct pole); 

        // Functions for fermi_pole

        // This is the central evaluation of the fermi and bose functions.  It also
        // includes evaluating arbitrary derivatives of the functions.  See frk
        // function that is rather complicated .  This function is also the MOST
        // challenging function for numerical evaluation.  It is the most likely
        // source of issue or errors and some thought should go into testing this
        // carefully
        std::complex<double> fermi_pole(ami_parms &parms, pole_struct pole,
                                        ami_vars external);

        void print_pole_struct_info(pole_struct g);

        /// Recursive construction of fermi_bose derivatives.
        std::complex<double> fermi_bose(int m, double sigma, double beta,
                                  std::complex<double> E);

        void print_alpha_info(alpha_t alpha);

        /// Given a set of external energies, beta, and frequencies, will evaluate the
        /// energy of a pole_struct.
        std::complex<double> get_energy_from_pole(pole_struct pole,
                                                    ami_vars external);

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

        void print_epsilon_info(epsilon_t eps);

        // Functions for Factorize_Rn

        /// Optimize function for SPR notation.
        void factorize_Rn(Ri_t &Rn, g_prod_t &unique_g, R_ref_t &Rref,
                        ref_eval_t &Eval_list);

        bool g_struct_equiv(g_struct &g1, g_struct &g2, int &sign);

        void reduce_rref(R_ref_t &Rref, ref_eval_t &Eval_list);

        bool pair_v_equiv(ref_v_t &r1, ref_v_t &r2, int &r1sign, int &r2sign);

        




    private:


};