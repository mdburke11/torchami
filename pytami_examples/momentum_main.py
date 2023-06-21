import torch
import pytami
import examples as ex
import momenta_integrand as mom
import numpy as np
import torchquad

def main():

    # calculates a second order self energy diagram using torchquad
    
    # init device
    device = torch.device("cpu")
    ami = pytami.TamiBase(device)

    # init diagram info
    R0 = ex.construct_example2()
    avars = ex.construct_ext_example2(ami)

    # setup helpers
    ftout = pytami.ft_terms()
    E_REG = 0
    N_INT = 2
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

    k = torch.rand([100, 4])
    print(integrand(k))



if __name__ == "__main__":
    main()