import torch
import pytami
import copy
import external as ext


def main() -> None:
    graph_type: pytami.TamiBase.graph_type = pytami.TamiBase.Sigma
    seed: int = 0

    g: pytami.TamiGraph = pytami.TamiGraph(graph_type, seed)

    order: int = 2
    group: int = 0
    n: int = 0
    folder: str = "../examples/ggm_examples/"

    second_ord_sigma, pref = read_diagram(g, folder, order, group, n)

    R0: pytami.TamiBase.g_prod_t = pytami.TamiBase.g_prod_t()
    g.trojan_graph_to_R0(second_ord_sigma, R0)

    print("\nThis is the second order self energy")
    print(f"Prefactor: {pref}")
    print_R0(R0)
    print()

    # now try the CT diagrams
    max_ct = 2
    ct_R0, prefactors = read_ct_diagrams(g, second_ord_sigma, order, max_ct)

    for i, r in enumerate(ct_R0):
        print(f"\nGraph: {i+1}")
        print(f"Prefactor: {prefactors[i]}")
        print_R0(r)


def print_R0(R0: pytami.TamiBase.g_prod_t) -> None:

    x: pytami.TamiBase.g_struct

    print("Alpha:")
    for x in R0:
        a: int
        for a in x.alpha_:
            print(a, end=" ")
        print()

    print()

    print("Epsilon:")
    for x in R0:
        e: int
        for e in x.eps_:
            print(e, end=" ")
        print()


def read_ct_diagrams(
        ptg: pytami.TamiGraph, diagram: pytami.TamiGraph.trojan_graph,
        order: int,
        max_ct: int) -> tuple[list[pytami.TamiBase.g_prod_t], list[float]]:

    ct_temp: pytami.TamiGraph.graph_vector = pytami.TamiGraph.graph_vector()
    R0_vec: list[pytami.TamiBase.g_prod_t] = [
    ]  # place to store all final R0 objects
    R0: pytami.TamiBase.g_prod_t = pytami.TamiBase.g_prod_t(
    )  # temp R0 object to unpack R0's

    ptg.trojan_generate_sigma_ct(diagram, ct_temp, max_ct)

    # extract to python friendly trojan_graphs
    trojan_results: list[pytami.TamiGraph.trojan_graph] = [
        pytami.TamiGraph.trojan_graph(ct_temp, i) for i in range(len(ct_temp))
    ]
    prefactors: list[float] = []

    ct_g: pytami.TamiGraph.trojan_graph
    for ct_g in trojan_results:
        ptg.trojan_graph_to_R0(ct_g, R0)
        R0_vec.append(copy.copy(R0))  # avoid mutability of g_prod_t objects

        prefactor: float = ptg.trojan_get_prefactor(ct_g, order)
        prefactors.append(prefactor)

    return R0_vec, prefactors


def read_diagram(ptg: pytami.TamiGraph, folder: str, ord: int, group: int,
                 n: int) -> tuple[pytami.TamiGraph.trojan_graph, float]:

    # reads in the specified graph from the folder provided

    ggm: pytami.TamiGraph.gg_matrix_t = pytami.TamiGraph.gg_matrix_t()

    ptg.read_ggmp(folder, ggm, ord, ord + 1)
    ptg.print_ggm(ggm)

    print("After removing pairs")
    ptg.print_ggm(ggm)
    ptg.ggm_label(ggm, 0)

    gg: pytami.TamiGraph.graph_group = ggm[ord][group].graph_vec
    diagram: pytami.TamiGraph.trojan_graph = pytami.TamiGraph.trojan_graph(
        gg, n)

    prefactor: float = ptg.trojan_get_prefactor(diagram, ord)

    return diagram, prefactor


def prep_ext(
    order: int,
    evars: ext.ext_vars,
    device: torch.device = torch.device("cpu")
) -> ext.ext_vars:

    default_batch_size: int = 10  # this is fixed when being integrating
    fbatch_size: int = 5

    energy: torch.Tensor = torch.zeros([default_batch_size, order + 1],
                                       device=device)
    frequency_vec: list[tuple[complex]] = [
        torch.tensor((0.0 + 0.0j, ) * order +
                     (1j * (2 * i + 1) * torch.pi / evars.beta, ),
                     device=device) for i in range(fbatch_size)
    ]
    frequency: pytami.TamiBase.frequency_t = torch.vstack(frequency_vec)
    external: pytami.TamiBase.ami_vars = pytami.TamiBase.ami_vars(
        energy, frequency, evars.beta)

    return external


if __name__ == "__main__":
    main()
