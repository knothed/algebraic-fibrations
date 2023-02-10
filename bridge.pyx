

#### GRAPH FIBERING ####

def colorings(g,c):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    #cdef arr2d_fixed abc = arr2d_fixed_create_from(adj.data,n*n,1);
    #print_arrf(abc);

    # Isometries
    cdef arr2d_fixed isos = get_isometries(adj)
    #print(f"num isos: {isos.len}")

    # Colorings: max #colors
    from sage.graphs.cliquer import all_cliques
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)
    cdef arr2d_var partitions = cliquewise_vertex_partition(n, cliques)

    now = time_ms()
    cdef arr2d_fixed cols = find_all_colorings(adj,c,partitions)
    print(f"found {cols.len} cols in {time_ms()-now} ms")
    #return np_array_from_arrf(cols)

    now = time_ms()
    cdef arr2d_fixed reduced_cols = kill_permutations_and_isos(n,c,cols,isos)
    print(f"reduced to {reduced_cols.len} in {time_ms()-now} ms")
    free_arrf(cols)
    return np_array_from_arrf(reduced_cols)




#### OTHER INTERFACE ####

def state_as_list(state: int):
    res = []
    cdef int i
    for i in range(log2_int(state)+1):
        i

#### UTILITIES ####







cdef extern from "coloring.c":
    int num_colors_upper_bound(int n, arr2d_var cliques, arr2d_fixed legal_states)

    arr2d_fixed find_all_colorings(arr2d_fixed adj, int num_cols, arr2d_var partitions)
    arr2d_fixed kill_permutations_and_isos(int n, int num_colors, arr2d_fixed cols, arr2d_fixed isos)


    legal_orbits_result find_legal_orbits(int n, arr2d_fixed colorings, arr2d_fixed legal_states, int num_threads, bint stop_after_first)
    arr2d_fixed all_legal_states(arr2d_fixed adj)

cdef extern from "graph.c":
    bint is_graph_hyperbolic(arr2d_fixed adj)
    arr2d_fixed get_isometries(arr2d_fixed adj)


