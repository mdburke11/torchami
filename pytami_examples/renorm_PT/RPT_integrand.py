# This is a sample integrand class that could be used inconjunction with an external integration
# library to perform the remaining spatial integrations after the AMI process 

import torch
import pytami
from typing import Callable
from collections import Counter 
import external as ext

def z_const(k : torch.Tensor):
    const = 0.1 * 1j
    return const


    
class RPT_integrand:

    # takes all the parameters needed to evaluate the momentum integrand
    # as well as a python function/lambda for the particle dispersion of the physical
    # problem at hand.

    # the __call__ method is defined to evaluate the integrand for the current
    # external physical parameters with a provided n by m momentum torch tensor. 
    # Where n is the batch size of integrands to be simulateously evaluated and
    # m is the number of independent momentum variables (dim * (2*order - 1))


    def __init__(self, tami: pytami.TamiBase, R0: pytami.TamiBase.g_prod_t, avars: pytami.TamiBase.ami_vars, 
                ft: pytami.TamiBase.ft_terms, parms: pytami.TamiBase.ami_parms, eps: Callable[[torch.Tensor], torch.Tensor],
                z: Callable[[torch.Tensor], torch.Tensor], evalReal: bool, extern_vars: ext.ext_vars) -> None:
        
        self.tami = tami
        self.R0 = R0
        self.avars = avars
        self.ft = ft
        self.parms = parms

        self.eps = eps # free particle dispersion
        self.z = z # this is the momentum dep renormalized PT shift - this will be a callable function 

        self.evalReal = evalReal
        self.external_vars = extern_vars

        # insert the parameters from extern_vars into the correct internal location
        self.update_ext_vars(self.external_vars)

        # useful things to extract
        self.dim = len(self.external_vars.k) # 2D or 3D
        self.device = self.tami.getDevice()
        I = torch.eye(self.dim, device=self.device)
        alpha: torch.Tensor = torch.tensor([self.R0[i].alpha_ for i in range(len(self.R0))], device=self.device)
        self.full_alpha = torch.kron(alpha, I).T # tensor that will be multiplied when updating the energies
        self.order = self.parms.N_INT_ # number of integrals (order) and 2 * order - 1 = len(R0) (if not Renorm PT)
        self.num_prop = len(self.R0) # number of propagators, equal to 2*self.order-1 unless renorm PT problem

        # Count the power to raise z to in each call. First count the number of occurances for all alpha vectors. Then
        # fill a list with the powers - 1 (power to be raised) at the correct index of R0.

        alps : list[tuple[int]] = [tuple(self.R0[i].alpha_) for i in range(len(R0))] # extract tuples to compare
        indices : list[int] = [alps.index(i) for i in set(alps)] # indices of the unique g_prod_t in R0
        count : dict[tuple[int], int] = Counter(alps) # get counts of each propagator
        self.powers : list[int] = [0 for i in range(len(R0))]  
        for i in indices:
            self.powers[i] = count[alps[i]] - 1 # each power is the number of greens functions present - 1



    def eff_eps(self, k: torch.Tensor) -> torch.Tensor:
        # combines all the extra components to the particle dispersion together - chemical potential, renorm PT

        return self.external_vars.mu - torch.hstack([self.eps(k.split(self.dim, dim=1)[i]) for i in range(self.num_prop)])



    def update_ext_vars(self, extern_vars: ext.ext_vars) -> None:
        self.external_vars = extern_vars
        self.avars.frequency_[-1] = extern_vars.reW + 1j * extern_vars.imW # needs to be internally updated
        self.avars.BETA_ = extern_vars.beta # needs to be internally updated for fermi function evaluations

    def get_ext_vars(self) -> list[float | complex | list[float] | float]:
        return self.external_vars.get_ext_vars()

    def update_integrand(self, k: torch.Tensor) -> None:
        # takes new set of internal momenta (eg. k = rand(0, 2pi)) with external and gets correct lin. comb.
        # then evaluates the dispersion for each propagator in the diagram and inserts it 
        # into the avars object to then evaluate thediagram

        # append on the external momentum to the momentum tensor
        K: torch.Tensor = torch.hstack([k] + [torch.full((len(k), 1), self.external_vars.k[i], device=self.device) for i in range(self.dim)])

        # get linear combinations to get energies, now in form (k1_x, k1_y, k2x, k2y, ...)
        combined = torch.matmul(K, self.full_alpha)

        #apply dispersion to dim columns at a time and update AMI integrand
        self.avars.energy_ = self.eff_eps(combined)

    def __call__(self, x: torch.Tensor) -> torch.Tensor:
        # function which is called in integration: used as integrand = AMI_integrand(...); integrand(k: torch.tensor) # returns tensor of evaluated integrand
        self.update_integrand(x)
        value: torch.Tensor = self.tami.evaluate(self.parms, self.ft, self.avars)
        if self.evalReal:
            return value.real
        return value.imag

