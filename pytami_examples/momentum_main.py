import torch
import pytami
import examples as ex
import momenta_integrand as mom
import external as ext
import flat_integ as flat
import time

def main():
    mat_freq_flat_2ord()
    #mat_freq_flat_4ord()
    #mat_freq_flat_6ord()

def mat_freq_flat_2ord():
# init device
    device = torch.device("cuda")
    ami = pytami.TamiBase(device)

    # init diagram info
    R0 = ex.construct_example2()
    avars = ex.construct_ext_example2(ami)

    # setup helpers
    ftout = pytami.TamiBase.ft_terms()
    E_REG = 0
    N_INT = 2
    parms = pytami.TamiBase.ami_parms(N_INT, E_REG)
    ami.construct(N_INT, R0, ftout)

    # Get external values
    beta: float = 8.33
    mu: complex = 0.1
    k: list[float] = [torch.pi, 0.0]
    N_freq: int = 100
    frequencies: torch.tensor = torch.tensor([[0.0, 0.0, 1j * (2*i + 1) * torch.pi / beta] for i in range(0, N_freq)], device=device)
    avars.frequency_ = frequencies

    evars = ext.ext_vars(beta, mu, k)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, ext.epsilon_2D,
                                 False, evars)


    flat_mc = flat.flat_mc_integrator(device, N_freq ,Max_batch=int(1e5))
    
    with open("2ord.dat", "a") as f:

        t1 = time.time()
        integral_value = flat_mc.integrate(integrand, dim=4, N=10**9, 
                                integration_domain=[[0, 2*torch.pi]] * 4)
        t2 = time.time()

        if N_freq == 1:
            n=0
            print(f"{n} {frequencies[n][-1]} {integral_value.ans} {integral_value.error}")
            f.write(f"{n} {frequencies[-1]} {integral_value.ans} {integral_value.error}\n")
            f.flush()
        else:
            for n in range(len(integral_value.ans)):
                print(f"{n} {frequencies[n][-1]} {integral_value.ans[n]} {integral_value.error[n]}")
                f.write(f"{n} {frequencies[n][-1]} {integral_value.ans[n]} {integral_value.error[n]}\n")
                f.flush()
            f.close()
    print(f"time: {t2 - t1}")
"""
def mat_freq_flat_4ord():
    # init device
    device = torch.device("cuda")
    ami = pytami.TamiBase(device)

    # init diagram info
    R0 = ex.construct_example4()
    avars = ex.construct_ext_example4(ami)

    # setup helpers
    ftout = pytami.TamiBase.ft_terms()
    E_REG = 0
    N_INT = 4
    parms = pytami.TamiBase.ami_parms(N_INT, E_REG)
    ami.construct(N_INT, R0, ftout)

    # Get external values
    n: int = 1 # fermionic matsubara freq

    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [torch.pi, 0.0]
    reW: float = 0.0
    imW: float = (2*n + 1) * torch.pi / beta

    evars = ext.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, ext.epsilon_2D,
                                False, evars)


    flat_mc = flat.flat_mc_integrator(device, Max_batch=int(4e5))
    
    with open("4ord.dat", "a") as f:

        for n in range(0, 100, 7):

            evars.imW = (2*n + 1) * torch.pi / beta
            integrand.update_ext_vars(evars)
            print(integrand.get_ext_vars())

            integral_value = flat_mc.integrate(integrand, dim=8, N=10_999_999, 
                                    integration_domain=[[0, 2*torch.pi]] * 8)

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
    ftout = pytami.TamiBase.ft_terms()
    E_REG = 0
    N_INT = 6
    parms = pytami.TamiBase.ami_parms(N_INT, E_REG)
    ami.construct(N_INT, R0, ftout)

    # Get external values
    n: int = 1 # fermionic matsubara freq

    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [torch.pi, 0.0]
    reW: float = 0.0
    imW: float = (2*n + 1) * torch.pi / beta

    evars = ext.ext_vars(beta, mu, k, reW, imW)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, ext.epsilon_2D,
                                False, evars)


    flat_mc = flat.flat_mc_integrator(device, Max_batch=int(1e4))
    
    with open("6ord.dat", "a") as f:

        for n in range(0, 100, 8):

            evars.imW = (2*n + 1) * torch.pi / beta
            integrand.update_ext_vars(evars)
            print(integrand.get_ext_vars())

            integral_value = flat_mc.integrate(integrand, dim=12, N=199_999, 
                                    integration_domain=[[0, 2*torch.pi]] * 12)

            print(f"{n} {evars.imW} {integral_value.ans} {integral_value.error}")
            f.write(f"{n} {evars.imW} {integral_value.ans} {integral_value.error}\n")
            f.flush()
        f.close()




"""

if __name__ == "__main__":
    main()