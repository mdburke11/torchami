import torch
import pytami
import examples as ex
import momenta_integrand as mom
import numpy as np
from torchquad import MonteCarlo, Boole, set_up_backend

def main():

    # calculates a second order self energy diagram using torchquad
    
    # init torchquad
    set_up_backend("torch", data_type="float64")

    # init device
    device = torch.device("cuda")
    ami = pytami.TamiBase(device)

    # init diagram info
    R0 = ex.construct_example4()
    avars = ex.construct_ext_example4(ami)

    # setup helpers
    ftout = pytami.ft_terms()
    E_REG = 0
    N_INT = 4
    parms = pytami.TamiBase.ami_parms(N_INT, E_REG)
    ami.construct(N_INT, R0, ftout)

    # Get external values
    beta: float = 5.0
    mu: complex = 0.0
    k: list[float] = [np.pi, np.pi]
    reW: float = 0.0
    imW: float = 0.0

    evars = mom.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, mom.epsilon_2D,
                                False, False, evars)


    mc = MonteCarlo()

    integral_value = mc.integrate(integrand, dim=8, N=1_000_000, 
                                  integration_domain=[[0, 2*np.pi]] * 8,
                                  backend="torch")

    print(integral_value)
    #
    #veg = Boole()

    #integral_value = veg.integrate(integrand, dim=4, N=1_000_001, 
    #                              integration_domain=[[0, 2*torch.pi]] * 4,
    #                              backend="torch")

    #print(integral_value)




if __name__ == "__main__":
    main()