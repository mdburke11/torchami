import torch
import pytami
import examples as ex
import momenta_integrand as mom
import numpy as np
from torchquad import MonteCarlo, Boole, set_up_backend, VEGAS, Simpson, Trapezoid
import flat_integ as flat
import time

def main():
    mat_freq_flat_2ord()
    #mat_freq_flat_4ord()
    #mat_freq_flat_6ord()

    #torchquad_comp() # Montecarlo and vegas are working but the rest are not
    # although, the conditions for Boole and Simpsons were probably not met


    

def torchquad_comp():

    methods = {"vegas" : VEGAS(), "Monte Carlo" : MonteCarlo()}#, "boole" : Boole(),
               #"Simpson" : Simpson(), "Trapezoid" : Trapezoid()}

    print("torchquad results: \n\n")

    for method in methods:
        
        print(f"Now doing method: {method}")
        t1 = time.time()
        torchquad_ex(methods[method])
        t2 = time.time()
        print(f"{method} took: {t2 - t1} seconds\n\n")

    print("Flat MC results")
    t3 = time.time()
    flat_dist_ex()
    t4 = time.time()
    print(f"flat took: {t4 - t3} seconds")


def mat_freq_flat_2ord():
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
    mu: complex = 0.1
    k: list[float] = [np.pi, 0.0]
    reW: float = 0.0
    imW: float = (2*n + 1) * np.pi / beta

    evars = mom.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, mom.epsilon_2D,
                                 False, evars)


    flat_mc = flat.flat_mc_integrator(device, Max_batch=int(1e5))
    
    with open("2ord.dat", "a") as f:

        for n in range(100):

            evars.imW = (2*n + 1) * np.pi / beta
            integrand.update_ext_vars(evars)
            print(integrand.get_ext_vars())

            integral_value = flat_mc.integrate(integrand, dim=4, N=999_999, 
                                    integration_domain=[[0, 2*np.pi]] * 4)

            print(f"{n} {evars.imW} {integral_value.ans} {integral_value.error}")
            f.write(f"{n} {evars.imW} {integral_value.ans} {integral_value.error}\n")
            f.flush()
        f.close()

def mat_freq_flat_4ord():
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
    n: int = 1 # fermionic matsubara freq

    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [np.pi, 0.0]
    reW: float = 0.0
    imW: float = (2*n + 1) * np.pi / beta

    evars = mom.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, mom.epsilon_2D,
                                False, evars)


    flat_mc = flat.flat_mc_integrator(device, Max_batch=int(4e5))
    
    with open("4ord.dat", "a") as f:

        for n in range(0, 100, 7):

            evars.imW = (2*n + 1) * np.pi / beta
            integrand.update_ext_vars(evars)
            print(integrand.get_ext_vars())

            integral_value = flat_mc.integrate(integrand, dim=8, N=10_999_999, 
                                    integration_domain=[[0, 2*np.pi]] * 8)

            print(f"{n} {evars.imW} {integral_value.ans} {integral_value.error}")
            f.write(f"{n} {evars.imW} {integral_value.ans} {integral_value.error}\n")
            f.flush()
        f.close()

def mat_freq_flat_6ord():
    # init device
    device = torch.device("cuda")
    ami = pytami.TamiBase(device)

    # init diagram info
    R0 = ex.construct_example6()
    avars = ex.construct_ext_example6(ami)

    # setup helpers
    ftout = pytami.ft_terms()
    E_REG = 0
    N_INT = 6
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
                                False, evars)


    flat_mc = flat.flat_mc_integrator(device, Max_batch=int(1e4))
    
    with open("6ord.dat", "a") as f:

        for n in range(0, 100, 8):

            evars.imW = (2*n + 1) * np.pi / beta
            integrand.update_ext_vars(evars)
            print(integrand.get_ext_vars())

            integral_value = flat_mc.integrate(integrand, dim=12, N=199_999, 
                                    integration_domain=[[0, 2*np.pi]] * 12)

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

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, pytami.epsilon_2D_cpp,
                                False, evars)


    flat_mc = flat.flat_mc_integrator(device)

    integral_value = flat_mc.integrate(integrand, dim=4, N=1_199_999, 
                                  integration_domain=[[0, 2*np.pi]] * 4)

    print(f"Ans: {integral_value.ans}")
    print(f"Err: {integral_value.error}")
    #
    #veg = Boole()

    #integral_value = veg.integrate(integrand, dim=4, N=1_000_001, 
    #                              integration_domain=[[0, 2*torch.pi]] * 4,
    #                              backend="torch")

    #print(integral_value)

    
def torchquad_ex(mc):

    #mc = VEGAS() #MonteCarlo, Boole, set_up_backend, VEGAS, Simpson, Trapezoid

    # calculates a second order self energy diagram using torchquad
    
    # init torchquad
    set_up_backend("torch", data_type="float64")

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
    n: int = 1 # matsubara index
    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [np.pi, 0.0]
    reW: float = 0.0
    imW: float = (2*n + 1) * np.pi / beta

    evars = mom.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, pytami.epsilon_2D_cpp,
                                False, evars)

    integral_value = mc.integrate(integrand, dim=4, N=1_199_999, 
                                  integration_domain=[[0, 2*np.pi]] * 4,
                                  backend="torch")

    print(f"\nResult: {integral_value / (2*torch.pi)**4}")

    if isinstance(mc, VEGAS):
        print(f"Error: {mc._get_error()/ (2*torch.pi)**4}")


if __name__ == "__main__":
    main()