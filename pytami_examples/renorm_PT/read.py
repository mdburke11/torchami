import torch
import pytami

# script that imports graph R0 files from c++ in the convention
# even -> eps, odd -> alpha into R0 files

# Also has an initializer to setup the ami_vars for each graph

def R0_from_graph(fname: str):
    #fname = "../R0_graphs/R0_o2g0n0.graph"#"test.dat"
    f: file = open(fname, "r")
    count: int = 0
    R0: list = []
    g: list = [None, None]
 
    for l in f.readlines():

        v: list[int] = [int(i) for i in l.strip(' \n').split(" ")]
        g[count%2] = pytami.VectorInt(v) # even -> eps, odd -> alpha
        count += 1

        if count%2 == 0:
            R0.append(pytami.AmiBase.g_struct(g[0], g[1]))

    return pytami.g_prod_t(R0)

def initialize_ext_vars(R0: pytami.g_prod_t, ext_freq=0.0+0.0j: complex, beta=5.00: float, device=torch.cpu: torch.device):

    print(len(R0)) # number of Greens functions
    order_diag: int = (len(R0) + 1) // 2# order of digram n has 2n-1 greens functions 
    print(order_diag)

    # give a dummy energy thats zeros batchsize 10
    energy: torch.tensor = torch.zeros([default_batchsize, order_diag], device=dev)
    frequency: pytami.frequency_t = pytami.frequency_([0.0+0.0j for i in range(order_diag)])
    frequency.append(ext_freq)
    external: pytami.TamiBase.ami_vars = pytami.TamiBase.ami_vars(energy, frequency, beta)

    return external

def main():
    fname: str = "R0_graphs/R0_o2g0n0.graph"#"test.dat"
    R0: pytami.g_prod_t = R0_from_graph(fname)

    # test if R0 is correct
    for i in R0:
        print(i.eps_)
        print(i.eps_[0], i.eps_[1], i.eps_[2], type(i.eps_[0]))
        print(i.alpha_)
    
    print("all good!?")

if __name__ == "__main__":
    main()