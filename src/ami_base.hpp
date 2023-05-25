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

namespace AmiBase{

    class place_holder_base_class{

        public:



        private:
            // Maximum argument allowed for fermi and bose functions - to prevent inf due
            // to double precision numbers
            double exp_max_arg = 500.0;


            bool drop_bosonic_diverge = false; // if set to true, then E_REG not needed because bosonic
             // divergences will simply be set to zero.  Rigorously this may
             // not be correct.

            bool drop_der = false; // if set to true, then E_REG not needed because
                         // bosonic divergences will simply be set to zero.
                         // Rigorously this may not be correct.

            bool drop_matsubara_poles = false; // if set to true, ignores Matsubara 
                        // poles with zero energy 
            
            // bool is_real_external=false;
            bool zero_external_w = false;
            bool overflow_detected = false;
            bool verbose = false; // flag for verbose output - primarily for debugging

            double precision_cutoff =
                1e15; // By default this is set to roughly the precision of double.  If
                      // values exceed this then numerical overflow is virtually
                      // guaranteed

            ext_type ext_freq_type = matsubara;

            


        
    };

    /// Returns the sign of a value - or zero if it is uniquely zero.
    template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

    /// Simple factorial function. It is rarely called except for multipole
    /// problems, and typically for only small arguments.
    int factorial(int n);

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




    








}