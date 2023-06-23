import torch
import pytami
import examples as ex
import momenta_integrand as mom
import numpy as np
#from torchquad import MonteCarlo, Boole, set_up_backend
import flat_integ as flat
    
def main():

    #torchquad_ex()
    #flat_dist_ex()
    mat_freq_flat()


def mat_freq_flat():
# init device
    device = torch.device("cuda")
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
    n: int = 1 # fermionic matsubara freq

    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [np.pi, 0.0]
    reW: float = 0.0
    imW: float = (2*n + 1) * np.pi / beta

    evars = mom.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, mom.epsilon_2D,
                                False, False, evars)


    flat_mc = flat.flat_mc_integrator(device, Max_batch=10**6)
    
    with open("2ord.dat", "a") as f:

        for n in range(100):

            evars.imW = (2*n + 1) * np.pi / beta
            integrand.update_ext_vars(evars)

            print("start")
            integral_value = flat_mc.integrate(integrand, dim=4, N=999_999, 
                                    integration_domain=[[0, 2*np.pi]] * 4)
            print("stop")

            print(f"{n} {evars.imW} {integral_value.ans} {integral_value.error}")
            f.write(f"{n} {evars.imW} {integral_value.ans} {integral_value.error}\n")
            f.flush()
        f.close()

def flat_dist_ex():

    # calculates a second order self energy diagram using flat distribution sampling

    # init device
    device = torch.device("cuda")
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
    n: int = 1 # fermionic matsubara freq

    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [np.pi, 0.0]
    reW: float = 0.0
    imW: float = (2*n + 1) * np.pi / beta

    evars = mom.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, mom.epsilon_2D,
                                False, False, evars)


    flat_mc = flat.flat_mc_integrator(device)

    integral_value = flat_mc.integrate(integrand, dim=4, N=1_999_999, 
                                  integration_domain=[[0, 2*np.pi]] * 4)

    print(f"Ans: {integral_value.ans}")
    print(f"Err: {integral_value.error}")
    #
    #veg = Boole()

    #integral_value = veg.integrate(integrand, dim=4, N=1_000_001, 
    #                              integration_domain=[[0, 2*torch.pi]] * 4,
    #                              backend="torch")

    #print(integral_value)

    
def torchquad_ex():

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
    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [np.pi, 0.0]
    reW: float = 0.0
    imW: float = 0.0

    evars = mom.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, mom.epsilon_2D,
                                False, False, evars)


    mc = MonteCarlo()

    integral_value = mc.integrate(integrand, dim=8, N=1_999_999, 
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