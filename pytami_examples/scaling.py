import torch
import pytami
import examples as ex
import numpy as np
import time

def construct_4ord(tami, batch_size):
    template_energy = [1, 1.1, 1.2, 1.31, 1.4, 0.01, 0.1]
    energy_size = len(template_energy)

    energy = 8.0 * torch.rand([batch_size, energy_size], device=tami.getDevice()) - 4.0
    frequency = pytami.frequency_t()
    for i in range(5):
        frequency.append(0+0j)
    
    frequency.append(0+np.pi*1j)
    beta = 1.0
    external = pytami.TamiBase.ami_vars(energy, frequency, beta)

    return external

def scaling():

    device = torch.device("cpu")
    ofname = "py_scaling.dat"

    tami = pytami.TamiBase(device)
    R0 = ex.construct_example4()
    avars = construct_4ord(tami, 5)
    N_INT=4
    E_REG=0
    amiparms = pytami.TamiBase.ami_parms(0,0)

    print("Starting energy tensor: ")
    print(avars.energy_)

    ftout = pytami.ft_terms()

    tami.construct(N_INT, R0, ftout)

    data = []

    for size in range(1, 10000, 1000):#range(100000, 1000000, 100):
        avars = construct_4ord(tami, size)

        t1 = time.time()
        result = tami.evaluate(amiparms, ftout, avars)
        t2 = time.time()

        diff = (t2 - t1) * 1000000
        
        data.append((size, diff))

    np.savetxt(ofname, data)

def main():
    scaling()


if __name__ == "__main__":
    main()


    



