from sage.all import *
from libc.stdlib cimport malloc, realloc, free

#### GRAPH FIBERING ####

def is_hyperbolic(g):
    cdef int n = g.order()
    cdef int* adj = list_to_array(g.adjacency_matrix()._list())
    res = is_graph_hyperbolic(n,adj);
    free(adj)
    return res

def graph_fiberings(g, max_cols=None, verbose=True):
    cdef int n = g.order()
    cdef int* adj = list_to_array(g.adjacency_matrix()._list())

    # Legal states
    if verbose: now = time_ms()
    cdef int legal_count = 0
    cdef int* legal_states = all_legal_states(n,adj,&legal_count)
    if verbose: print(f"Get {legal_count} legal states: took {time_ms()-now} ms")

    # Isometries
    if verbose: now = time_ms()
    cdef int isos_count = 0
    cdef int* isos = get_isometries(n, adj, &isos_count)
    if verbose: print(f"Get {isos_count} isometries: took {time_ms()-now} ms")

    # Colorings: max #colors
    from sage.graphs.cliquer import all_cliques
    if verbose: now = time_ms()
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr_of_arrs cliques = list_of_variable_length_lists_to_array(cliques_py)
    if verbose: print(f"Convert {len(cliques_py)} cliques to C: took {time_ms()-now} ms")

    if verbose: now = time_ms()
    cdef int cmax = max_possible_colors(n, cliques, legal_states, legal_count)
    if verbose: print(f"Get {cmax} max possible colors: took {time_ms()-now} ms")

    cdef int cmin = g.chromatic_number()
    if max_cols is not None: cmax = min(cmax, max_cols)
    if verbose: print(f"#colors: between {cmin} and {cmax}")

    # Colorings: all colorings of specific #colors
    cdef int i, c, reduced_count, orbit_count = 0
    cdef int* cols
    cdef int* reduced_cols
    cdef int* orbit
    fibered = False
    for c in range(cmin,cmax+1):
        # Get all colorings and reduce by color swapping and graph isometries
        from sage.graphs.graph_coloring import all_graph_colorings
        if verbose: now = time_ms()
        cols_py = list(all_graph_colorings(g,c,vertex_color_dict=True))
        cols = list_of_dicts_to_array(cols_py,n)
        reduced_count = 0
        reduced_cols = kill_permutations_and_isos(n,c,cols,len(cols_py),isos,isos_count,&reduced_count)

        # Find legal orbits
        for i in range(reduced_count):
            orbit_count = 0
            orbit = find_legal_orbits(n,reduced_cols+n*i,legal_states,legal_count,&orbit_count)
            if (orbit_count > 0):
                print(f"found legal orbit! col: {list_from_array(reduced_cols+n*i,n)}, states: {list_from_array(orbit,orbit_count)}")
                fibered = True

        if verbose: print(f"Checking colorings with {c} colors: took {time_ms()-now} ms")

        free(cols)
        free(reduced_cols)
        free(orbit)

    free(adj)
    free(legal_states)
    free(isos)
    free_arr(cliques)

    return fibered

## TOOLS ##

cdef extern from "utils.c":
    ctypedef struct arr_of_arrs:
        int* data
        int* start_indices
        int len
    void free_arr(arr_of_arrs arr)
    arr_of_arrs create_arr_of_arrs(int max_total_elems, int max_arrs)
    void print_arr(arr_of_arrs arr)

cdef extern from "coloring.c":
    int max_possible_colors(int n, arr_of_arrs cliques, int* legal_states, int legal_count)
    int* kill_permutations_and_isos(int n, int num_cols, int* cols, int cols_count, int* isos, int isos_count, int* result_count)

cdef extern from "legal.c":
    int* all_legal_states(int n, int* adj_matrix, int* result_count)
    bint is_state_legal(int n, int* adj_matrix, int state)
    int* find_legal_orbits(int n, int* coloring, int* legal_states, int num_states, int* result_count)

cdef extern from "graph.c":
    bint is_graph_hyperbolic(int n, int* adj)
    int* get_isometries(int n, int* adj, int* result_count)

def time_ms():
    import time
    return time.time_ns() // 1000000

cdef int* list_to_array(list):
    cdef int l = len(list)
    cdef int* res = <int*> malloc(l*sizeof(int))
    for i in range(l):
        res[i] = list[i]
    return res

cdef arr_of_arrs list_of_variable_length_lists_to_array(list):
    cdef int total_len = 0
    cdef int i,j
    for i in range(len(list)):
        total_len += len(list[i])

    cdef arr_of_arrs result = create_arr_of_arrs(total_len, len(list))
    result.len = len(list);

    cdef int current_count = 0
    for i in range(len(list)):
        l = list[i]
        result.start_indices[i] = current_count
        for j in range(len(l)):
            result.data[current_count+j] = l[j]
        current_count += len(l)
    result.start_indices[len(list)] = current_count

    return result

cdef int* list_of_dicts_to_array(colorings, n):
    cdef int l = len(colorings)
    cdef int* res = <int*> malloc(n*l*sizeof(int))
    for i in range(l):
        col = colorings[i]
        for k,v in col.items():
            res[n*i+k] = v
    return res

cdef int* list_of_lists_to_array(isometries, n):
    cdef int l = len(isometries)
    cdef int* res = <int*> malloc(n*l*sizeof(int))
    for i in range(l):
        for j in range(n):
            res[n*i+j] = isometries[i][j]
    return res

cdef list_of_lists_from_array(int* cols, int num_cols, int n):
    result = []
    for i in range(num_cols):
        col = []
        for j in range(n):
            col.append(cols[n*i+j])
        result.append(col)
    return result

cdef list_from_array(int* cols, int n):
    result = []
    for i in range(n):
        result.append(cols[i])
    return result