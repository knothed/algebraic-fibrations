###
### A Python interface around the C code.
### Contains all exposed functionality regarding legal orbit search and related functionality.
###

from sage.all import *
from sage.graphs.cliquer import all_cliques
from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as np
np.import_array()
import ctypes

#### VIRTUAL ALGEBRAIC FIBERING: LEGAL ORBIT SEARCH ####

# Determine (by brute force) whether a graph virtually algebraically fibers. Only consider colorings with up to max_cols colors, if desired.
# Split up the orbit search work onto multiple threads, if desired.
def has_legal_orbit(g, max_cols=None, verbose=True, num_threads=1):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)

    if max_cols is None: max_cols = 0
    cdef legal_orbits_result orbits = graph_fiberings(adj, cliques, max_cols, verbose, max(1,num_threads), true)
    fibers = orbits.colorings.len > 0

    free_arrf(adj, orbits.colorings)
    free_arrv(cliques, orbits.states)
    return fibers

# Get a single legal orbit (which has at most max_cols colors, if desired).
# Split up the orbit search work onto multiple threads, if desired.
# If such an orbit exists, the result will be of the form {'coloring': list, 'state': list} where coloring is a list representing the coloring,
# and state represents a SINGLE state: it contains all vertices which are part of the state.
# Else, None is returned.
def one_legal_orbit(g, max_cols=None, verbose=True, num_threads=1):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)

    if max_cols is None: max_cols = 0
    cdef legal_orbits_result orbits = graph_fiberings(adj, cliques, max_cols, verbose, max(1,num_threads), true)

    result = None
    if orbits.colorings.len > 0:
        col = list_from_arrf_row(orbits.colorings, 0)
        state = vertices_in_state(get_arrv(orbits.states, 0, 0))
        result = {'coloring': col, 'state': state}

    free_arrf(adj, orbits.colorings)
    free_arrv(cliques, orbits.states)
    return result

# Get all legal orbits that exist for a graph. Stop at colorings with more than max_cols colors, if desired.
# Split up the orbit search work onto multiple threads, if desired.
# The result will be an array of dictionaries {'coloring': list, 'states': list} where coloring is a list representing the coloring,
# and states is a list representing all the different legal orbits that exist: it contains one state per legal orbit.
def all_legal_orbits(g, max_cols=None, verbose=True, num_threads=1):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)

    if max_cols is None: max_cols = 0
    cdef legal_orbits_result orbits = graph_fiberings(adj, cliques, max_cols, verbose, max(1,num_threads), false)

    # Convert to list
    now = time_ms()
    result = [0] * orbits.colorings.len
    cdef int i
    for i in range(orbits.colorings.len):
        result[i] = {}
        result[i]['coloring'] = list_from_arrf_row(orbits.colorings, i)
        result[i]['states'] = list_from_arrv_row(orbits.states, i)

    free_arrf(adj, orbits.colorings)
    free_arrv(cliques, orbits.states)
    return result


#### RELATED FUNCTIONALITY #####

# Determine whether the RACG given by a graph is a word-hyperbolic group.
def is_hyperbolic(g):
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    res = is_graph_hyperbolic(adj);
    free_arrf(adj)
    return res

# Given a state as an integer, get the list of vertices in the state.
def vertices_in_state(state):
    res = []
    cdef int i
    for i in range(32):
        if (state >> i) & 1:
            res.append(i)
    return res

# Find all colorings for a graph with exactly num_cols colors.
# Reduce the colorings modulo color swapping and graph isometries.
# The result is a 2D numpy array where each row represents a coloring.
def all_reduced_colorings(g, num_cols, verbose=False):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cdef arr2d_fixed isos = get_isometries(adj)

    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)
    cdef arr2d_var partitions = cliquewise_vertex_partition(n, cliques)

    if verbose: t1 = time_ms()
    cdef arr2d_fixed cols = find_all_colorings(adj,num_cols,partitions)
    if verbose: t2 = time_ms()
    cdef arr2d_fixed reduced = reduce_colorings(n,num_cols,cols,isos)
    if verbose: print(f"Found {cols.len} colorings in {t2-t1}ms; reduced to {reduced.len} unique ones in {time_ms()-t2} ms.")

    free_arrf(adj, isos, cols)
    free_arrv(cliques, partitions)
    return np_array_from_arrf(reduced)


# ...


#### C IMPORTS ####

cdef extern from "impl/coloring.c":
    arr2d_var cliquewise_vertex_partition(int n, arr2d_var cliques); # todo: remove this function entirely and just focus on a single largest clique
    arr2d_fixed find_all_colorings(arr2d_fixed adj, int num_cols, arr2d_var partitions)
    arr2d_fixed reduce_colorings(int n, int num_colors, arr2d_fixed cols, arr2d_fixed isos)

cdef extern from "impl/fibering.c":
    legal_orbits_result graph_fiberings(arr2d_fixed adj, arr2d_var cliques, int max_cols, bint verbose, int num_threads, bint single_orbit)

cdef extern from "impl/graph.c":
    bint is_graph_hyperbolic(arr2d_fixed adj)
    arr2d_fixed get_isometries(arr2d_fixed adj)

cdef extern from "impl/utils.c":
    ctypedef struct arr2d_fixed:
        int* data
        int row_len
        int len
    int get_arrf(arr2d_fixed arr, int i, int j)
    void free_arrf(...)
    void print_arrf(arr2d_fixed arr)
    arr2d_fixed arr2d_fixed_create_from(int* data, int row_len, int len)

    ctypedef struct arr2d_var:
        pass
    int get_arrv(arr2d_var arr, int i, int j)
    int size_arrv(arr2d_var arr, int i)
    void free_arrv(...)
    void print_arrv(arr2d_var arr)
    arr2d_var arr2d_var_create_empty(int total_capacity, int num_rows_capacity)
    arr2d_var append_arrv(arr2d_var arr, int* src, int n)

cdef extern from "impl/legal.c":
    ctypedef struct legal_orbits_result:
        arr2d_fixed colorings
        arr2d_var states

#### ARRAY CONVERSION ####

# todo: test speed of list_to_array etc.
def time_ms():
    import time
    return time.time_ns() // 1000000

## PYTHON --> C ##

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

# Convert a list into a C array.
cdef int* list_to_array(list):
    cdef int* arr = <int*> malloc(len(list)*sizeof(int))
    cdef int i=0
    for i in range(len(list)):
        arr[i]=list[i]
    return arr

## C --> PYTHON ##

# Convert an arr2d_fixed to a Python numpy array.
# Do not free the original array after calling this!
cdef object np_array_from_arrf(arr2d_fixed arr):
    cdef np.npy_intp shape[2]
    shape[0] = <np.npy_intp> arr.len
    shape[1] = <np.npy_intp> arr.row_len
    ndarray = np.PyArray_SimpleNewFromData(2, shape, np.NPY_INT, <void*>arr.data)
    np.PyArray_UpdateFlags(ndarray, ndarray.flags.num | np.NPY_OWNDATA)
    return ndarray

cdef list list_from_arrf_row(arr2d_fixed arr, int i):
    res: list = [0] * arr.row_len
    cdef int j
    for j in range(arr.row_len):
        res[j] = get_arrf(arr,i,j)
    return res

cdef list list_from_arrv_row(arr2d_var arr, int i):
    res: list = [0] * size_arrv(arr,i)
    cdef int j
    for j in range(size_arrv(arr,i)):
        res[j] = get_arrv(arr,i,j)
    return res