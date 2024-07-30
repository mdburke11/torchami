import sys
sys.path.append('helperScripts')
import torch
import pytami
from getDevice import getDevice
import examples as ex
import momenta_integrand as mom
import external as ext
import flat_integ as flat
import time

# example of evaluating a diagram on the imaginary axis.
# default main evaluates the imaginary part of the second 
# order self energy diagram for the first 20 matsubara 
# frequencies using 10^5 monte carlo evaluations.

# 4th and 6th order diagram functions also exist.


def main():
    mat_freq_flat_2ord()

def mat_freq_flat_2ord():
    # init device
    device = getDevice()
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
    N_freq: int = 20
    frequencies: torch.tensor = torch.tensor(
        [[0.0, 0.0, 1j * (2 * i + 1) * torch.pi / beta]
         for i in range(0, N_freq)],
        device=device)
    avars.frequency_ = frequencies

    evars = ext.ext_vars(beta, mu, k)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, ext.epsilon_2D,
                                  False, evars)

    flat_mc = flat.flat_mc_integrator(device, N_freq, Max_batch=int(1e5))

    t1 = time.time()
    integral_value = flat_mc.integrate(
        integrand,
        dim=4,
        N=10**6,
        integration_domain=[[0, 2 * torch.pi]] * 4)
    t2 = time.time()

    with open("2ord.dat", "w") as f:
        if N_freq == 1:
            n = 0
            print(
                f"{n} {frequencies[n][-1]} {integral_value.ans} {integral_value.error}"
            )
            f.write(
                f"{n} {frequencies[-1]} {integral_value.ans} {integral_value.error}\n"
            )
            f.flush()
        else:
            for n in range(len(integral_value.ans)):
                print(
                    f"{n} {frequencies[n][-1]} {integral_value.ans[n]} {integral_value.error[n]}"
                )
                f.write(
                    f"{n} {frequencies[n][-1]} {integral_value.ans[n]} {integral_value.error[n]}\n"
                )
                f.flush()
    print(f"time: {t2 - t1}")


def mat_freq_flat_4ord():
    # init device
    device = getDevice()
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
    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [torch.pi, 0.0]
    N_freq: int = 1
    frequencies: torch.tensor = torch.tensor(
        [[0.0, 0.0, 0.0, 0.0, 1j * (2 * i + 1) * torch.pi / beta]
         for i in range(0, N_freq)],
        device=device)
    avars.frequency_ = frequencies

    evars = ext.ext_vars(beta, mu, k)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, ext.epsilon_2D,
                                  False, evars)

    flat_mc = flat.flat_mc_integrator(device, N_freq, Max_batch=int(1e5))

    t1 = time.time()
    integral_value = flat_mc.integrate(
        integrand,
        dim=8,
        N=10**8,
        integration_domain=[[0, 2 * torch.pi]] * 8)
    t2 = time.time()

    with open("4ord.dat", "w") as f:
        if N_freq == 1:
            n = 0
            print(
                f"{n} {frequencies[n][-1]} {integral_value.ans} {integral_value.error}"
            )
            f.write(
                f"{n} {frequencies[-1]} {integral_value.ans} {integral_value.error}\n"
            )
            f.flush()
        else:
            for n in range(len(integral_value.ans)):
                print(
                    f"{n} {frequencies[n][-1]} {integral_value.ans[n]} {integral_value.error[n]}"
                )
                f.write(
                    f"{n} {frequencies[n][-1]} {integral_value.ans[n]} {integral_value.error[n]}\n"
                )
                f.flush()
    print(f"time: {t2 - t1}")


def mat_freq_flat_6ord():
    # init device
    device = getDevice()
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

    beta: float = 8.33
    mu: complex = 0.0
    k: list[float] = [torch.pi, 0.0]
    N_freq: int = 1
    frequencies: torch.tensor = torch.tensor(
        [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1j * (2 * i + 1) * torch.pi / beta]
         for i in range(0, N_freq)],
        device=device)
    avars.frequency_ = frequencies

    evars = ext.ext_vars(beta, mu, k)

    integrand = mom.AMI_integrand(ami, R0, avars, ftout, parms, ext.epsilon_2D,
                                  False, evars)

    flat_mc = flat.flat_mc_integrator(device, N_freq, Max_batch=int(1e4))

    t1 = time.time()
    integral_value = flat_mc.integrate(
        integrand,
        dim=12,
        N=10**6,
        integration_domain=[[0, 2 * torch.pi]] * 12)
    t2 = time.time()

    with open("6ord.dat", "w") as f:
        if N_freq == 1:
            n = 0
            print(
                f"{n} {frequencies[n][-1]} {integral_value.ans} {integral_value.error}"
            )
            f.write(
                f"{n} {frequencies[-1]} {integral_value.ans} {integral_value.error}\n"
            )
            f.flush()
        else:
            for n in range(len(integral_value.ans)):
                print(
                    f"{n} {frequencies[n][-1]} {integral_value.ans[n]} {integral_value.error[n]}"
                )
                f.write(
                    f"{n} {frequencies[n][-1]} {integral_value.ans[n]} {integral_value.error[n]}\n"
                )
                f.flush()
    print(f"time: {t2 - t1}")


if __name__ == "__main__":
    main()
