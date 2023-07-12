import sys
sys.path.append("../")
import torch
import pytami
import graph_util as gru
import RPT_integrand as rpt
import external as ext
import flat_integ as flat


def main() -> None:
    second_ord_c1_a()


def second_ord_c1_a() -> None:

    #
    # graph set up
    #
    
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

    max_ct : int = 2    
    ct_R0 : list[pytami.TamiBase.g_prod_t]; prefactors : list[float]
    ct_R0, prefactors = gru.read_ct_diagrams(g, second_ord_sigma, order, max_ct)

    diag_index : int = 5
    R0 = ct_R0[diag_index] # lets evaluate this diagram
    prefactor = prefactors[diag_index]

    print("We will evaluate this diagram\n")
    gru.print_R0(R0)
    print() 

    #
    # integrand setup
    #
    for _ in range(10):
        device : torch.device = torch.device("cuda")
        ami : pytami.TamiBase = pytami.TamiBase(device)

        # external vars
        mat_freq : int = 10
        beta : float = 8.33
        mu : complex = 0
        k : list[float] = [torch.pi, 0]
        reW = 0.0
        imW = (2*mat_freq + 1) * torch.pi / beta

        evars = ext.ext_vars(beta, mu, k, reW, imW)
        avars : pytami.TamiBase.ami_vars = gru.prep_ext(order, evars, device)

        # integrand
        ftout : pytami.TamiBase.ft_terms = pytami.TamiBase.ft_terms()
        E_REG : float = 0
        N_INT : int = order
        parms : pytami.TamiBase.ami_parms = pytami.TamiBase.ami_parms(N_INT, E_REG)
        ami.construct(N_INT, R0, ftout)

        integrand : rpt.RPT_integrand = rpt.RPT_integrand(ami, R0, prefactor, avars, ftout, parms, ext.epsilon_2D, rpt.dampt_eps, True, evars)

        # integrate
        flat_mc = flat.flat_mc_integrator(device)

    
        value = flat_mc.integrate(integrand, dim=4, N=6_999_999,
                                integration_domain=[[0, 2*torch.pi]] * 4)

        print(f"Ans: {value.ans}")
        print(f"Err: {value.error}")


    



if __name__ == "__main__":
    main()