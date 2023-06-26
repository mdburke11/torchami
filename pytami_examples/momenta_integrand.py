# This is a sample integrand class that could be used inconjunction with an external integration
# library to perform the remaining spatial integrations after the AMI process 

import torch
import pytami
from typing import Callable

import time # TODO remove this and timing info

class ext_vars:

    def __init__(self, beta: float, mu: complex, k: list[float], reW: float, imW: float):
        self.beta = beta
        self.mu = mu
        self.k = k
        self.reW = reW
        self.imW = imW

    def get_ext_vars(self):
        return [self.beta, self.mu, self.k, self.reW, self.imW]

    
def epsilon_2D(k: torch.tensor) -> torch.tensor: # must return column vector of energies (use .unsqueeze(1))
    #return -2 * (torch.cos(k[:,0]) + torch.cos(k[:,1])).unsqueeze(1)
    return -2 * torch.cos(k).sum(dim=1, keepdim=True)


    
class AMI_integrand:

    # takes all the parameters needed to evaluate the momentum integrand
    # as well as a python function/lambda for the particle dispersion of the physical
    # problem at hand.

    # the __call__ method is defined to evaluate the integrand for the current
    # external physical parameters with a provided n by m momentum torch tensor. 
    # Where n is the batch size of integrands to be simulateously evaluated and
    # m is the number of independent momentum variables (dim * len(alpha))


    def __init__(self, tami: pytami.TamiBase, R0: pytami.g_prod_t, avars: pytami.TamiBase.ami_vars, 
                ft: pytami.ft_terms, parms: pytami.TamiBase.ami_parms, eps: Callable[[torch.tensor], torch.tensor],
                RenormPT: bool, evalReal: bool, extern_vars: ext_vars) -> None:
        
        self.tami = tami
        self.R0 = R0
        self.avars = avars
        self.ft = ft
        self.parms = parms

        self.eps = eps # free particle dispersion
        self.evalReal = evalReal
        self.RenormPT = RenormPT # flag will modify the integrand to have renormalized PT correction on integrand TODO: implement this
        self.external_vars = extern_vars

        # insert the parameters from extern_vars into the correct internal location
        self.update_ext_vars(self.external_vars)

        # useful things to extract
        self.dim = len(self.external_vars.k) # 2D or 3D
        self.device = self.tami.getDevice()
        I = torch.eye(self.dim, device=self.device)
        self.alpha: torch.tensor = torch.tensor([self.R0[i].alpha_ for i in range(len(self.R0))], device=self.device)
        self.full_alpha = torch.kron(self.alpha, I).T # tensor that will be multiplied when updating the energies
        self.order = self.parms.N_INT_ # number of integrals (order) and 2 * order - 1 = len(R0) (if not Renorm PT)



    def update_ext_vars(self, extern_vars: ext_vars) -> None:
        self.external_vars = extern_vars
        self.avars.frequency_[-1] = extern_vars.reW + 1j * extern_vars.imW # needs to be internally updated
        self.avars.BETA_ = extern_vars.beta # needs to be internally updated for fermi function evaluations

    def get_ext_vars(self) -> ext_vars:
        return self.external_vars.get_ext_vars()

    def update_integrand(self, k: torch.tensor) -> None:
        # takes new set of internal momenta (eg. k = rand(0, 2pi)) with external and gets correct lin. comb.
        # then evaluates the dispersion for each propagator in the diagram and inserts it 
        # into the avars object to then evaluate thediagram

        # append on the external momentum to the momentum tensor
        t1 = time.time()
        K: torch.tensor = torch.hstack([k] + [torch.full((len(k), 1), self.external_vars.k[i], device=self.device) for i in range(self.dim)])
        #K: torch.tensor = torch.hstack([k]+ [torch.tensor([torch.full((len(k), 1), self.external_vars.k[i], device=self.device) for i in range(self.dim)])])
        #K: torch.tensor = torch.cat((k, torch.tensor(self.external_vars.k, device=self.device).repeat(len(k), 1)), 1)
        t2 = time.time()
        # get linear combinations to get energies, now in form (k1_x, k1_y, k2x, k2y, ...)
        combined = torch.matmul(K, self.full_alpha)#.to(torch.complex64) # convert to complex values here to combine with mu
        t3 = time.time()
        # get the energies into the correct tensor format [e1, e2, e3] where e_i are column vectors - note minus sign since Greens functions are missing one
        self.avars.energy_ = self.external_vars.mu - torch.hstack([self.eps(combined[:, range(i*self.dim, self.dim*(i+1))]) for i in range(2*self.order -1)])
        t4 = time.time()

        print("\n\nupdate_integrand_interior:")
        print(f"get K: {t2 - t1}")
        print(f"product: {t3 - t2}")
        print(f"eval eps: {t4 - t3}")
        print()

    def __call__(self, x: torch.tensor) -> torch.torch:
        # function which is called in integration: used as integrand = AMI_integrand(...); integrand(k: torch.tensor) # returns tensor of evaluated integrand
        t1 = time.time()
        self.update_integrand(x)
        t2 = time.time()
        value: torch.tensor = self.tami.evaluate(self.parms, self.ft, self.avars)
        t3 = time.time()
        print(f"update integrand: {t2 - t1}")
        print(f"evaluate_val: {t3 - t2}")
        if self.evalReal:
            return value.real

        return value.imag

