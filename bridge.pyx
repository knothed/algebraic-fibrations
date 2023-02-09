from sage.all import *
from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as np
np.import_array()
import ctypes

#### GRAPH FIBERING ####

def colorings(g,c):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cdef arr2d_fixed abc = arr2d_fixed_create_from(adj.data,n*n,1);
    print_arrf(abc);

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

def is_hyperbolic(g):
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    res = is_graph_hyperbolic(adj);
    free_arrf(adj)
    return res

def graph_fiberings(g, max_cols=None, verbose=True, num_threads=1):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())

    # Preparations

    if verbose: print(f"Preparations ...")
    if verbose: now = time_ms()

    cdef arr2d_fixed legal_states = all_legal_states(adj)
    cdef arr2d_fixed isos = get_isometries(adj)

    from sage.graphs.cliquer import all_cliques
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)
    cdef int cmax = num_colors_upper_bound(n, cliques, legal_states)
    cdef arr2d_var partitions = cliquewise_vertex_partition(n, cliques)

    cdef int cmin = 0 # g.chromatic_number() REMOVE BECAUSE SLOW!
    if max_cols is not None: cmax = min(cmax, max_cols)

    if verbose: print(f"... took {time_ms()-now} ms: {legal_states.len} legal states, {isos.len} isos; {cmin} <= #colors <= {cmax}")

    # Colorings: all colorings of specific #colors

    cdef int i, c = 0
    cdef arr2d_fixed cols
    cdef arr2d_fixed orbits
    fibers = False

    for c in range(cmin,cmax+1):
        # Get all colorings and reduce by color swapping and graph isometries
        if verbose: now = time_ms()
        cols = find_all_colorings(adj,c,partitions)
        if verbose: now2 = time_ms()
        #cols = kill_permutations_and_isos(n,c,cols,isos)
        if verbose: now3 = time_ms()

        # Find legal orbits
        if verbose:
            print(f"{c} colors: Testing {cols.len} colorings ...")
        orbits = find_legal_orbits(n,cols,legal_states,num_threads)
        if (orbits.len > 0):
            #np_array_from_arrf(orbits)
            print(f"found legal orbit! col: todo, states: {orbits.len}")
            fibers = True

        if verbose:
            print(f"... took {time_ms()-now} ms (all_cols): {now2-now}, reduce: {now3-now2}, orbits: {time_ms()-now3}")

        free_arrf(cols)
        free_arrf(orbits)

    free_arrf(adj)
    free_arrf(legal_states)
    free_arrf(isos)
    free_arrv(cliques)
    free_arrv(partitions)

    return fibers


def time_ms():
    import time
    return time.time_ns() // 1000000


#### C IMPORTS ####

cdef extern from "utils.c":
    ctypedef struct arr2d_fixed:
        int* data
        int row_len
        int len
    void free_arrf(arr2d_fixed arr)
    void print_arrf(arr2d_fixed arr)
    arr2d_fixed arr2d_fixed_create_from(int* data, int row_len, int len)

    ctypedef struct arr2d_var:
        int* data
        int* end_indices
        int len
    void free_arrv(arr2d_var arr)
    void print_arrv(arr2d_var arr)
    arr2d_var arr2d_var_create_empty(int total_capacity, int num_rows_capacity)
    arr2d_var append_arrv(arr2d_var arr, int* src, int n)

cdef extern from "coloring.c":
    int num_colors_upper_bound(int n, arr2d_var cliques, arr2d_fixed legal_states)
    arr2d_var cliquewise_vertex_partition(int n, arr2d_var cliques)
    arr2d_fixed find_all_colorings(arr2d_fixed adj, int num_cols, arr2d_var partitions)
    arr2d_fixed kill_permutations_and_isos(int n, int num_colors, arr2d_fixed cols, arr2d_fixed isos)

cdef extern from "legal.c":
    ctypedef struct legal_orbits_result:
        arr2d_fixed colorings
        arr2d_var states
    legal_orbits_result find_legal_orbits(int n, arr2d_fixed colorings, arr2d_fixed legal_states, int num_threads, bint stop_after_first)
    arr2d_fixed all_legal_states(arr2d_fixed adj)

cdef extern from "graph.c":
    bint is_graph_hyperbolic(arr2d_fixed adj)
    arr2d_fixed get_isometries(arr2d_fixed adj)


#### ARRAY CONVERSION ####

# todo: test speed of list_to_array etc.

# Convert a graph adjacency matrix into C form.
cdef arr2d_fixed arrf_from_adj_matrix(adj):
    cdef int n = adj.dimensions()[0]
    cdef int* data = list_to_array(adj._list())
    return arr2d_fixed_create_from(data,n,n)

# Convert a graph adjacency matrix into C form.
cdef arr2d_var arrv_from_nested_list(list):
    cdef int total_len = sum([len(sub) for sub in list])
    cdef arr2d_var res = arr2d_var_create_empty(total_len, len(list));
    cdef int* p
    for sub in list:
        p = list_to_array(sub)
        res = append_arrv(res,p,len(sub))
        free(p)
    return res

# Convert a C array to a Python numpy array.
cdef object np_array_from_arrf(arr2d_fixed arr):
    cdef np.npy_intp shape[2]
    shape[0] = <np.npy_intp> arr.len
    shape[1] = <np.npy_intp> arr.row_len
    ndarray = np.PyArray_SimpleNewFromData(2, shape, np.NPY_INT, <void*>arr.data)
    np.PyArray_UpdateFlags(ndarray, ndarray.flags.num | np.NPY_OWNDATA)
    return ndarray

# Convert a list into a C array.
cdef int* list_to_array(list):
    cdef int* arr = <int*> malloc(len(list)*sizeof(int))
    cdef int i=0
    for i in range(len(list)):
        arr[i]=list[i]
    return arr