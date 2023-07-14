import sys
sys.path.append("../")
import torch
import pytami
import graph_util as gru
import RPT_integrand as rpt
import external as ext
import flat_integ as flat


def main() -> None:
    mat_axis()
    #second_ord_c1_a()


def constant_z(k : torch.Tensor) -> torch.Tensor:
    const : complex = 0.1j
    return torch.full([len(k), 1], const, device=k.device)

def mat_axis() -> None:

    MAX_MAT_FREQ : int = 20

    # init diagram info
    graph_type: pytami.TamiBase.graph_type = pytami.TamiBase.Sigma
    seed: int = 0 

    g: pytami.TamiGraph = pytami.TamiGraph(graph_type, seed)

    order : int = 2; group : int = 0; n : int = 0
    folder : str = "../../examples/ggm_examples/"

    second_ord_sigma, prefactor = gru.read_diagram(g, folder, order, group, n)

    R0 : pytami.TamiBase.g_prod_t = pytami.TamiBase.g_prod_t()
    g.trojan_graph_to_R0(second_ord_sigma, R0)

    print("\nThis is the second order self energy")
    gru.print_R0(R0)
    print() 

    max_ct : int = 1   
    ct_R0 : list[pytami.TamiBase.g_prod_t]; prefactors : list[float]
    ct_R0, prefactors = gru.read_ct_diagrams(g, second_ord_sigma, order, max_ct)

    All_R0 = [R0] + ct_R0 # add the base diagram as well
    All_prefactors = [prefactor] + prefactors
    
    for i in range(1, len(All_R0)):
        R0 = All_R0[i] # lets evaluate this diagram
        prefactor = All_prefactors[i]

        fname = f"2ord_{i}.dat"

        print("We will evaluate this diagram\n")
        gru.print_R0(R0)
        print()

        # init device
        device : torch.device = torch.device("cuda")
        ami : pytami.TamiBase = pytami.TamiBase(device)

        # external vars
        mat_freq : int = 1
        beta : float = 5.00
        mu : complex = 0.0
        k : list[float] = [torch.pi, 0]
        reW = 0.0
        imW = (2*mat_freq + 1) * torch.pi / beta

        evars = ext.ext_vars(beta, mu, k, reW, imW)
        avars : pytami.TamiBase.ami_vars = gru.prep_ext(order, evars, device)

        # setup helpers
        ftout = pytami.TamiBase.ft_terms()
        E_REG = 0
        N_INT = 2
        parms = pytami.TamiBase.ami_parms(N_INT, E_REG)
        ami.construct(N_INT, R0, ftout)

        re_integrand = rpt.RPT_integrand(ami, R0, prefactor, avars, ftout, parms, ext.epsilon_2D,
                                    constant_z, True, evars)
        im_integrand = rpt.RPT_integrand(ami, R0, prefactor, avars, ftout, parms, ext.epsilon_2D,
                                    constant_z, False, evars)
        #integrand = rpt.RPT_integrand(ami, R0, prefactor, avars, ftout, parms, ext.epsilon_2D,
        #                            lambda k : torch.zeros([len(k), 1], device=k.device), False, evars)

        flat_mc = flat.flat_mc_integrator(device, Max_batch=int(1e5))
        
        with open(fname, "a") as f:

            for n in range(MAX_MAT_FREQ):

                evars.imW = (2*n + 1) * torch.pi / beta
                re_integrand.update_ext_vars(evars)
                im_integrand.update_ext_vars(evars)
                #print(integrand.get_ext_vars())

                N=770_000_000

                re_integral_value = flat_mc.integrate(re_integrand, dim=4, N=N, 
                                        integration_domain=[[0, 2*torch.pi]] * 4)
                im_integral_value = flat_mc.integrate(im_integrand, dim=4, N=N, 
                                        integration_domain=[[0, 2*torch.pi]] * 4)

                print(f"{n} {evars.imW} {re_integral_value.ans} {re_integral_value.error} {im_integral_value.ans} {im_integral_value.error}")
                f.write(f"{n} {evars.imW} {re_integral_value.ans} {re_integral_value.error} {im_integral_value.ans} {im_integral_value.error}\n")
                f.flush()
            f.close()



    



if __name__ == "__main__":
    main()