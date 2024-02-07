# This is a sample integrand class that could be used inconjunction with an external integration
# library to perform the remaining spatial integrations after the AMI process

# this code will integrate a diagram for a spectrum of frequencies provided rather than iterating
# through a list. All the frequencies will be provided in the avars on initialization

import torch
import pytami
from typing import Callable, Union
import external as ext


class AMI_integrand:

    # takes all the parameters needed to evaluate the momentum integrand
    # as well as a python function/lambda for the particle dispersion of the physical
    # problem at hand.

    # the __call__ method is defined to evaluate the integrand for the current
    # external physical parameters with a provided n by m momentum torch tensor.
    # Where n is the batch size of integrands to be simulateously evaluated and
    # m is the number of independent momentum variables (dim * (2*order - 1))

    def __init__(self, tami: pytami.TamiBase, R0: pytami.TamiBase.g_prod_t,
                 avars: pytami.TamiBase.ami_vars, ft: pytami.TamiBase.ft_terms,
                 parms: pytami.TamiBase.ami_parms,
                 eps: Callable[[torch.tensor], torch.tensor], evalReal: bool,
                 extern_vars: ext.ext_vars) -> None:

        self.tami = tami
        self.R0 = R0
        self.avars = avars
        self.ft = ft
        self.parms = parms

        self.eps = eps  # free particle dispersion
        self.evalReal = evalReal
        self.external_vars = extern_vars

        # insert the parameters from extern_vars into the correct internal location
        self.update_ext_vars(self.external_vars)

        # useful things to extract
        self.dim = len(self.external_vars.k)  # 2D or 3D
        self.device = self.tami.getDevice()
        I = torch.eye(self.dim, device=self.device)
        alpha: torch.tensor = torch.tensor(
            [self.R0[i].alpha_ for i in range(len(self.R0))],
            device=self.device)
        self.full_alpha = torch.kron(
            alpha,
            I).T  # tensor that will be multiplied when updating the energies
        self.order = self.parms.N_INT_  # number of integrals (order) and 2 * order - 1 = len(R0) (if not Renorm PT)

    def update_ext_vars(self, extern_vars: ext.ext_vars) -> None:
        self.external_vars = extern_vars
        self.avars.BETA_ = extern_vars.beta  # needs to be internally updated for fermi function evaluations

    def get_ext_vars(self) -> list[Union[float, complex, list[float]]]:
        return self.external_vars.get_ext_vars()

    def update_integrand(self, k: torch.tensor) -> None:
        # takes new set of internal momenta (eg. k = rand(0, 2pi)) with external and gets correct lin. comb.
        # then evaluates the dispersion for each propagator in the diagram and inserts it
        # into the avars object to then evaluate thediagram

        # append on the external momentum to the momentum tensor
        K: torch.tensor = torch.hstack([k] + [
            torch.full(
                (len(k), 1), self.external_vars.k[i], device=self.device)
            for i in range(self.dim)
        ])

        # get linear combinations to get energies, now in form (k1_x, k1_y, k2x, k2y, ...)
        combined = torch.matmul(K, self.full_alpha)

        #apply dispersion to dim columns at a time and update AMI integrand
        self.avars.energy_ = self.external_vars.mu - torch.hstack([
            self.eps(combined.split(self.dim, dim=1)[i])
            for i in range(2 * self.order - 1)
        ])

    def __call__(self, x: torch.tensor) -> torch.torch:
        # function which is called in integration: used as integrand = AMI_integrand(...); integrand(k: torch.tensor) # returns tensor of evaluated integrand
        self.update_integrand(x)
        value: torch.tensor = self.tami.evaluate(self.parms, self.ft,
                                                 self.avars)
        if self.evalReal:
            return value.real
        return value.imag
