import torch
import pytami
import numpy as np

def construct_example2():

    g = pytami.TamiBase.g_prod_t()

    alpha1 = pytami.TamiBase.VectorInt([1, 0, 0])
    alpha2 = pytami.TamiBase.VectorInt([0, 1, 0])
    alpha3 = pytami.TamiBase.VectorInt([-1, 1, 1])

    epsilon1 = pytami.TamiBase.VectorInt([1, 0, 0])
    epsilon2 = pytami.TamiBase.VectorInt([0, 1, 0])
    epsilon3 = pytami.TamiBase.VectorInt([0, 0, 1])

    g1 = pytami.TamiBase.g_struct(epsilon1, alpha1)
    g2 = pytami.TamiBase.g_struct(epsilon2, alpha2)
    g3 = pytami.TamiBase.g_struct(epsilon3, alpha3)

    R0 = pytami.TamiBase.g_prod_t([g1, g2, g3])

    return R0

def construct_ext_example2(tami: pytami.TamiBase) -> pytami.TamiBase.ami_vars:

    ebatchsize = 10
    fbatchsize = 5

    energy = torch.tensor([-4, 0.1, -1], device=tami.getDevice()).repeat([ebatchsize, 1])
    
    frequency_vec = []
    for i in range(2):
        frequency_vec.append(0+0j)
    
    frequency_vec.append(0+np.pi*1j)
    beta = 1.0

    frequency = torch.tensor(frequency_vec, device=tami.getDevice()).repeat([fbatchsize, 1])
    external = pytami.TamiBase.ami_vars(energy, frequency, beta)

    return external


def construct_example4():
    
    g = pytami.TamiBase.g_prod_t()

    alpha1 = pytami.TamiBase.VectorInt([1, 0, 0, 1, -1])
    alpha2 = pytami.TamiBase.VectorInt([0, 0, 0, 1, 0])
    alpha3 = pytami.TamiBase.VectorInt([1, 0, 0, 0, 0])
    alpha4 = pytami.TamiBase.VectorInt([1, 0, 0, 0, 0])
    alpha5 = pytami.TamiBase.VectorInt([0, 1, 0, 0, 0])
    alpha6 = pytami.TamiBase.VectorInt([-1, 1, 1, 0, 0])
    alpha7 = pytami.TamiBase.VectorInt([0, 0, 1, 0, 0])   

    epsilon1 = pytami.TamiBase.VectorInt([0,0,0,0,0,1,0])
    epsilon2 = pytami.TamiBase.VectorInt([0,0,0,0,1,0,0])
    epsilon3 = pytami.TamiBase.VectorInt([1,0,0,0,0,0,0])
    epsilon4 = pytami.TamiBase.VectorInt([1,0,0,0,0,0,0])
    epsilon5 = pytami.TamiBase.VectorInt([0,1,0,0,0,0,0])
    epsilon6 = pytami.TamiBase.VectorInt([0,0,0,1,0,0,0])
    epsilon7 = pytami.TamiBase.VectorInt([0,0,1,0,0,0,0])



    g1 = pytami.TamiBase.g_struct(epsilon1, alpha1)
    g2 = pytami.TamiBase.g_struct(epsilon2, alpha2)
    g3 = pytami.TamiBase.g_struct(epsilon3, alpha3)
    g4 = pytami.TamiBase.g_struct(epsilon4, alpha4)
    g5 = pytami.TamiBase.g_struct(epsilon5, alpha5)
    g6 = pytami.TamiBase.g_struct(epsilon6, alpha6)
    g7 = pytami.TamiBase.g_struct(epsilon7, alpha7)


    R0 = pytami.TamiBase.g_prod_t([g1, g2, g3, g4, g5, g6, g7])

    return R0

def construct_ext_example4(tami: pytami.TamiBase) -> pytami.TamiBase.ami_vars:

    ebatchsize = 10
    fbatchsize = 5

    energy = torch.tensor([1, 1.1, 1.2, 1.31, 1.4, 0.01, 0.1], device=tami.getDevice()).repeat([ebatchsize, 1])
    frequency_vec = []
    for i in range(4):
        frequency_vec.append(0+0j)
    
    frequency_vec.append(0+np.pi*1j)
    beta = 1.0

    frequency = torch.tensor(frequency_vec, device=tami.getDevice()).repeat([fbatchsize, 1])
    external = pytami.TamiBase.ami_vars(energy, frequency, beta)

    return external


def construct_ext_example6(tami: pytami.TamiBase) -> pytami.TamiBase.ami_vars:

    ebatchsize = 10000
    fbatchsize = 100

    energy = torch.tensor([1,1.1,1.2,1.3,1.4,0, 0.1, 0.2, 0.3,0.4, 0.5], device=tami.getDevice()).repeat([ebatchsize, 1])
    frequency_vec = []
    for i in range(6):
        frequency_vec.append(0+0j)
    
    frequency_vec.append(0+np.pi*1j)
    beta = 1.0

    frequency = torch.tensor(frequency_vec, device=tami.getDevice()).repeat([fbatchsize, 1])
    external = pytami.TamiBase.ami_vars(energy, frequency, beta)

    return external

def construct_example6():
    
    g = pytami.TamiBase.g_prod_t()

    alpha1=pytami.TamiBase.VectorInt([1,0,0,0,0,0,0])
    alpha2=pytami.TamiBase.VectorInt([0,1,0,0,0,0,0])
    alpha3=pytami.TamiBase.VectorInt([0,0,1,0,0,0,0])
    alpha4=pytami.TamiBase.VectorInt([0,0,0,1,0,0,0])
    alpha5=pytami.TamiBase.VectorInt([0,0,0,0,1,0,0])
    alpha6=pytami.TamiBase.VectorInt([0,0,0,0,0,1,0])
    alpha7=pytami.TamiBase.VectorInt([1,-1,0,0,0,1,0])
    alpha8=pytami.TamiBase.VectorInt([-1,1,0,0,0,0,1])
    alpha9=pytami.TamiBase.VectorInt([-1,1,0,0,0,0,1])
    alpha10=pytami.TamiBase.VectorInt([-1,1,-1,1,0,0,1])
    alpha11=pytami.TamiBase.VectorInt([1,0,0,0,-1,1,0])

    epsilon1=pytami.TamiBase.VectorInt([1,0,0,0,0,0,0,0,0,0,0])
    epsilon2=pytami.TamiBase.VectorInt([0,1,0,0,0,0,0,0,0,0,0])
    epsilon3=pytami.TamiBase.VectorInt([0,0,1,0,0,0,0,0,0,0,0])
    epsilon4=pytami.TamiBase.VectorInt([0,0,0,1,0,0,0,0,0,0,0])
    epsilon5=pytami.TamiBase.VectorInt([0,0,0,0,1,0,0,0,0,0,0])
    epsilon6=pytami.TamiBase.VectorInt([0,0,0,0,0,1,0,0,0,0,0])
    epsilon7=pytami.TamiBase.VectorInt([0,0,0,0,0,0,1,0,0,0,0])
    epsilon8=pytami.TamiBase.VectorInt([0,0,0,0,0,0,0,1,0,0,0])
    epsilon9=pytami.TamiBase.VectorInt([0,0,0,0,0,0,0,1,0,0,0])
    epsilon10=pytami.TamiBase.VectorInt([0,0,0,0,0,0,0,0,0,1,0])
    epsilon11=pytami.TamiBase.VectorInt([0,0,0,0,0,0,0,0,0,0,1])

    g1 = pytami.TamiBase.g_struct(epsilon1, alpha1)
    g2 = pytami.TamiBase.g_struct(epsilon2, alpha2)
    g3 = pytami.TamiBase.g_struct(epsilon3, alpha3)
    g4 = pytami.TamiBase.g_struct(epsilon4, alpha4)
    g5 = pytami.TamiBase.g_struct(epsilon5, alpha5)
    g6 = pytami.TamiBase.g_struct(epsilon6, alpha6)
    g7 = pytami.TamiBase.g_struct(epsilon7, alpha7)
    g8 = pytami.TamiBase.g_struct(epsilon8, alpha8)
    g9 = pytami.TamiBase.g_struct(epsilon9, alpha9)
    g10 = pytami.TamiBase.g_struct(epsilon10, alpha10)
    g11 = pytami.TamiBase.g_struct(epsilon11, alpha11)


    R0 = pytami.TamiBase.g_prod_t([g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11])

    return R0

