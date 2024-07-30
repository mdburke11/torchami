import sys
sys.path.append('helperScripts')
from graph_util import print_R0
import torch
import pytami
import copy

# 3 examples: graph_eg is a plain load and print the R0 object (main)
# renorm_PT_eg is generates the counter term diagrams (here as a test, not in main)
# bose_alphas_eg is generates the bosonic propagator's momenta alphas (here as a test, not in main)

def graph_eg() -> None:

    graph_dir: str = "../examples/ggm_examples/"

    min_ord: int = 4
    max_ord: int = 6

    graph_type: pytami.TamiBase.graph_type = pytami.TamiBase.Sigma
    print(graph_type)

    seed: int = 0

    g: pytami.TamiGraph = pytami.TamiGraph(graph_type, seed)

    ggm: pytami.TamiGraph.gg_matrix_t = pytami.TamiGraph.gg_matrix_t()

    g.read_ggmp(graph_dir, ggm, min_ord, max_ord)
    g.print_ggm(ggm)
    g.ggm_remove_pairs(ggm, 0)

    print("After removing pairs")
    g.print_ggm(ggm)
    g.ggm_label(ggm, 0)

    R0: pytami.TamiBase.g_prod_t = pytami.TamiBase.g_prod_t()

    o: int
    for o in range(min_ord, max_ord + 1, 2):

        test: pytami.TamiGraph.graph_group = ggm[o][0].graph_vec

        trojan_g: pytami.TamiGraph.trojan_graph = pytami.TamiGraph.trojan_graph(
            test, 0)  # Wrapped graph_t object (Boost obj in disguise)

        g.trojan_graph_to_R0(trojan_g, R0)
        prefactor: float = g.trojan_get_prefactor(trojan_g, o)

        print(f"\n\nPrefactor: {prefactor}")
        print_R0(R0)


def renorm_PT_eg() -> None:

    # example python script that gets all the second order self energy diagrams up to 2 insertions

    graph_dir: str = "../examples/ggm_examples/"

    min_ord: int = 4
    max_ord: int = 6
    o: int = 2  # only second order

    max_ct: int = 2

    graph_type: pytami.TamiBase.graph_type = pytami.TamiBase.Sigma
    seed: int = 0

    g: pytami.TamiGraph = pytami.TamiGraph(graph_type, seed)

    ggm: pytami.TamiGraph.gg_matrix_t = pytami.TamiGraph.gg_matrix_t()

    g.read_ggmp(graph_dir, ggm, min_ord, max_ord)
    g.print_ggm(ggm)
    #g.ggm_remove_pairs(ggm, 0)

    print("After removing pairs")
    g.print_ggm(ggm)
    g.ggm_label(ggm, 0)

    gg: pytami.TamiGraph.graph_group = ggm[o][0].graph_vec
    secondOrdSigma: pytami.TamiGraph.trojan_graph = pytami.TamiGraph.trojan_graph(
        gg, 0)

    ct_temp: pytami.TamiGraph.graph_vector = pytami.TamiGraph.graph_vector(
    )  # catch all ct graphs

    R0_vec: list[pytami.TamiBase.g_prod_t] = [
    ]  # place to store all final R0 objects
    prefactors: list[float] = []  # place to store all prefactors
    R0: pytami.TamiBase.g_prod_t = pytami.TamiBase.g_prod_t(
    )  # temp R0 object to unpack R0's

    g.trojan_generate_sigma_ct(secondOrdSigma, ct_temp, max_ct)

    # extract to python friendly trojan_graphs
    trojan_results: list[pytami.TamiGraph.trojan_graph] = [
        pytami.TamiGraph.trojan_graph(ct_temp, i) for i in range(len(ct_temp))
    ]

    ct_g: pytami.TamiGraph.trojan_graph
    for ct_g in trojan_results:
        g.trojan_graph_to_R0(ct_g, R0)
        R0_vec.append(copy.copy(R0))  # avoid mutability of g_prod_t objects
        pf: float = g.trojan_get_prefactor(ct_g, o)  # order = 2
        prefactors.append(pf)

    #print all the ct graphs
    print(
        f"\n\nCT diagrams Second order self energy diagram up to {max_ct} insertions\n\n"
    )

    i: int
    r: pytami.TamiBase.g_prod_t
    for i, r in enumerate(R0_vec):
        print(f"--------------------------")
        print(f"graph: {i}")
        print(f"Prefactor: {prefactors[i]}\n")
        print_R0(r)

def bose_alphas_eg():

    # testing of the trojan implementation of the extract_bose_alphas function

    # load in the diagrams from min_ord to max_ord and print them and the bose alphas

    graph_dir: str = "../examples/ggm_examples/"

    min_ord: int = 2
    max_ord: int = 6

    graph_type: pytami.TamiBase.graph_type = pytami.TamiBase.Sigma
    print(f"Graph type: {graph_type}")

    seed: int = 0

    g: pytami.TamiGraph = pytami.TamiGraph(graph_type, seed)

    ggm: pytami.TamiGraph.gg_matrix_t = pytami.TamiGraph.gg_matrix_t()

    g.read_ggmp(graph_dir, ggm, min_ord, max_ord)
    g.print_ggm(ggm)
    g.ggm_remove_pairs(ggm, 0)

    print("After removing pairs")
    g.print_ggm(ggm)
    g.ggm_label(ggm, 0)

    # alpha_t containers to fill
    R0: pytami.TamiBase.g_prod_t = pytami.TamiBase.g_prod_t()
    bose_alphas: pytami.TamiBase.bose_alphas = pytami.TamiBase.bose_alphas()

    o: int
    for o in range(min_ord, max_ord + 1, 2):

        test: pytami.TamiGraph.graph_group = ggm[o][0].graph_vec

        trojan_g: pytami.TamiGraph.trojan_graph = pytami.TamiGraph.trojan_graph(
            test, 0)  # Wrapped graph_t object (Boost obj in disguise)

        # regular R0 graph
        g.trojan_graph_to_R0(trojan_g, R0)
        prefactor: float = g.trojan_get_prefactor(trojan_g, o)

        # new bosonic alphas
        g.trojan_extract_bose_alphas(trojan_g, bose_alphas)

        print(f"\n\nPrefactor: {prefactor}")
        print("Regular R0 Object:")
        print_R0(R0)
        print("Bosonic alphas")
        for a in bose_alphas:
            for x in a:
                print(x, end=" ")
            print()


if __name__ == "__main__":
    graph_eg()
    #renorm_PT_eg()
    #bose_alphas_eg()
