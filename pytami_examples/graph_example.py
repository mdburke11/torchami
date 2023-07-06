import torch
import pytami


def graph_eg():

    graph_dir = "../examples/ggm_examples/"

    max_ord = 6

    graph_type = pytami.TamiBase.Sigma
    print(graph_type)

    seed = 0 

    g = pytami.TamiGraph(graph_type, seed)

    ggm = pytami.TamiGraph.gg_matrix_t()

    g.read_ggmp(graph_dir, ggm, max_ord)
    g.print_ggm(ggm)
    g.ggm_label(ggm,0)

    R0 = pytami.g_prod_t()

    for o in range(2, max_ord+1, 2):

        test = ggm[o][0].graph_vec

        print(type(test))
        print(test)

        #g.graph_to_R0(ggm[2][0].graph_vec[0], R0)


        trojan_g = pytami.TamiGraph.trojan_graph(test, 0)

        g.trojan_graph_to_R0(trojan_g, R0)
        prefactor = g.trojan_get_prefactor(trojan_g, o)

        print(f"\n\nPrefactor: {prefactor}")

        print("Alpha:")
        for x in R0:
            print(x.alpha_)

        print("Epsilon:")
        for x in R0:
            print(x.eps_)


def renorm_PT_eg():
    pass





if __name__ == "__main__":
    graph_eg()