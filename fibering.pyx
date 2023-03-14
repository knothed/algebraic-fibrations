#cython: language_level=3
###
### A Python interface around the C code.
### Contains all exposed functionality regarding legal orbit search and related functionality.
###

from sage.all import *
from sage.graphs.cliquer import all_cliques
from libc.stdlib cimport malloc, free
from libc.stdint cimport uint64_t, int64_t, intptr_t

import numpy as np
cimport numpy as np
np.import_array()
import ctypes

from sage.matrix.matrix_integer_dense import Matrix_integer_dense
import subprocess

#### VIRTUAL ALGEBRAIC FIBERING: LEGAL ORBIT SEARCH ####

# Determine (by brute force) whether a graph virtually algebraically fibers. Only consider colorings with up to max_cols colors, if desired.
# Split up the parallelizable work into multiple threads, if desired.
def has_legal_orbit(g, min_cols=0, max_cols=0, verbose=True, total_progress_bar=True, num_threads=1):
    cdef int n = g.order()

    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)

    cdef legal_orbits_result orbits = graph_fiberings(adj, cliques, min_cols, max_cols, verbose, total_progress_bar, num_threads, true)
    fibers = orbits.colorings.len > 0

    free_arrf(adj, orbits.colorings)
    free_arrv(cliques, orbits.states)
    return fibers

# Get a single legal orbit (which has at most max_cols colors, if desired).
# Split up the parallelizable work into multiple threads, if desired.
# If such an orbit exists, the result will be of the form {'coloring': list, 'state': list} where coloring is a list representing the coloring,
# and state represents a SINGLE state: it contains all vertices which are part of the state.
# Else, None is returned.
def one_legal_orbit(g, min_cols=0, max_cols=0, verbose=True, total_progress_bar=False, num_threads=1):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)

    cdef legal_orbits_result orbits = graph_fiberings(adj, cliques, min_cols, max_cols, verbose, total_progress_bar, num_threads, true)

    result = None
    if orbits.colorings.len > 0:
        col = list_from_arrf_row(orbits.colorings, 0)
        state = vertices_in_state(get_arrv(orbits.states, 0, 0))
        result = {'coloring': col, 'state': state}

    free_arrf(adj, orbits.colorings)
    free_arrv(cliques, orbits.states)
    return result

# Get all legal orbits that exist for a graph. Stop at colorings with more than max_cols colors, if desired.
# Split up the parallelizable work into multiple threads, if desired.
# The result will be an array of dictionaries {'coloring': list, 'states': list} where coloring is a list representing the coloring,
# and states is a list representing all the different legal orbits that exist: it contains one state per legal orbit.
def all_legal_orbits(g, min_cols=0, max_cols=0, verbose=True, total_progress_bar=True, num_threads=1, states_as_vertex_lists=False):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)

    cdef legal_orbits_result orbits = graph_fiberings(adj, cliques, min_cols, max_cols, verbose, total_progress_bar, num_threads, false)

    # Convert to list
    result = [0] * orbits.colorings.len
    cdef int i
    for i in range(orbits.colorings.len):
        result[i] = {}
        result[i]['coloring'] = list_from_arrf_row(orbits.colorings, i)
        if states_as_vertex_lists:
            result[i]['states'] = list(map(vertices_in_state,list_from_arrv_row(orbits.states, i)))
        else:
            result[i]['states'] = list_from_arrv_row(orbits.states, i)

    free_arrf(adj, orbits.colorings)
    free_arrv(cliques, orbits.states)
    return result

# Find all legal orbits for a single coloring.
def legal_orbits_for_coloring(g, coloring_py):
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cdef int* coloring = list_to_array(coloring_py)

    cdef legal_orbits_result orbits = graph_fiberings_single_coloring(adj, coloring)
    result = [] if orbits.states.len == 0 else list_from_arrv_row(orbits.states, 0)

    free(coloring)
    free_arrf(adj, orbits.colorings)
    free_arrv(orbits.states)
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
def all_reduced_colorings(g, num_cols, verbose=False, num_threads=1):
    cdef int n = g.order()
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cdef arr2d_fixed isos = get_isometries(adj)

    cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
    cdef arr2d_var cliques = arrv_from_nested_list(cliques_py)
    cdef arr2d_var partitions = cliquewise_vertex_partition(n, cliques)

    if verbose: t1 = time_ms()
    cdef arr2d_fixed cols = find_all_colorings(adj, num_cols, partitions)
    if verbose: t2 = time_ms()
    cdef arr2d_fixed reduced = reduce_colorings(n, num_cols, cols, isos, num_threads)
    if verbose: print(f"Found { pretty_int(cols.len).decode('utf-8') } colorings in { pretty_ms(t2-t1,True).decode('utf-8') }; reduced to { pretty_int(reduced.len).decode('utf-8') } unique ones in { pretty_ms(time_ms()-t2,True).decode('utf-8') }.")

    free_arrf(adj, isos, cols)
    free_arrv(cliques, partitions)
    return np_array_from_arrf(reduced)


#### GENG: MULTIPLE GRAPHS ####

# Use geng to produce a stream of graphs, each of which is checked for fibering.
# When hyp_check, only hyperbolic graphs are checked.
# Returns all the graphs which fiber. Optionally, writes them to a file in graph6 format.
# Sensible args are: -c (connected), '{2*n-4}:0' (minimum edges)
def geng_fibering(n: int, geng: str, args: str, num_queues: int, queue_capacity: int, threads_per_queue: int, hyp_check: bint = 1, total: int = 0, write_to_file: str = ""):
    # Create fifo and start geng > fifo
    fifo = "/tmp/geng"
    try: os.remove(fifo)
    except OSError: pass
    os.mkfifo(fifo)
    subprocess.Popen([f'{geng} {n} -c {args} > {fifo}'], shell=True)

    # Preparations
    cdef arr2d_fixed adj;
    cdef arr2d_var cliques;
    cdef int i = 0
    cdef int step = max(1, total/100)

    cdef fibering_scheduler scheduler = make_scheduler(n, num_queues, queue_capacity, threads_per_queue, write_to_file.encode('UTF-8'))

    mat_space = Mat(ZZ,n)
    pytext = (f"Checking {pretty_int(total).decode('utf8')} graphs: ").encode('UTF-8')
    cdef char* text = pytext

    with open(fifo, 'r') as f:
        for line in iter(f.readline, ''):
            i += 1
            if (total > 0 and i % step == 0):
                print_progress(text, (<double>i)/(<double>total), -1)

            # Read matrix
            adj = arr2d_fixed_create_empty(n, n*n)
            adj.len = n
            read_adj_matrix_graph6(line.encode('ascii'), adj.data)
            if (hyp_check and not is_graph_hyperbolic(adj)): continue

            # Create Graph - required for all_cliques call.
            # This is quite slow, better: make own all_cliques algorithm
            g = graph_from_arrf(adj, mat_space)

            # C stuff is done multithreaded: add to scheduler
            cliques_py = sorted(list(all_cliques(g,min_size=2)), key=len, reverse=True)
            cliques = arrv_from_nested_list(cliques_py)

            add_to_scheduler(scheduler, adj, cliques)

    # Convert to graphs
    graphs = []
    cdef stream_result res = scheduler_finish(scheduler)
    res.results.len = n # tweak so that each graph_from_arrf call returns a single graph
    cdef int offset = n*n
    for i in range(res.num_fiber):
        g = graph_from_arrf(res.results,mat_space)
        # g = graph_from_graph6(graph6_from_adj_matrix(res.results).decode('ascii'), n)
        graphs.append(g)
        res.results.data = res.results.data + <intptr_t> offset
    return graphs

# Do not free the original adjacency matrix after calling this.
cdef graph_from_arrf(adj: arr2d_fixed, mat_space):
    l: list = [list(row) for row in np_array_from_arrf(adj)]
    mat = Matrix_integer_dense.__new__(Matrix_integer_dense, mat_space)
    mat.__init__(mat_space,l)
    return Graph(mat, format='adjacency_matrix')

# Convert Graph -> graph6
def graph6_from_graph(g):
    cdef arr2d_fixed adj = arrf_from_adj_matrix(g.adjacency_matrix())
    cdef char* res = graph6_from_adj_matrix(adj);
    free_arrf(adj)
    return res.decode('utf-8');

# Convert graph6 -> Graph
def graph_from_graph6(g6,n):
    cdef arr2d_fixed adj = arr2d_fixed_create_empty(n, n*n)
    adj.len = n
    read_adj_matrix_graph6(g6.encode('ascii'), adj.data)
    mat_space = Mat(ZZ,n)
    return graph_from_arrf(adj, mat_space)


#### C IMPORTS ####

cdef extern from "impl/utils.c":
    ctypedef struct arr2d_fixed:
        int* data
        int row_len
        int len
    int get_arrf(arr2d_fixed arr, int i, int j)
    void free_arrf(...)
    void print_arrf(arr2d_fixed arr)
    arr2d_fixed arr2d_fixed_create_from(int* data, int row_len, int len)
    arr2d_fixed arr2d_fixed_create_empty(int row_len, int capacity);

    ctypedef struct arr2d_var:
        int len
    int get_arrv(arr2d_var arr, int i, int j)
    int size_arrv(arr2d_var arr, int i)
    void free_arrv(...)
    void print_arrv(arr2d_var arr)
    arr2d_var arr2d_var_create_empty(int total_capacity, int num_rows_capacity)
    arr2d_var append_arrv(arr2d_var arr, int* src, int n)

    char* pretty_ms(uint64_t ms, bint subsecond_precision)
    char* pretty_int(int num)
    void print_progress(char* prefix, double progress, int64_t estimated_ms);

cdef extern from "impl/coloring.c":
    arr2d_var cliquewise_vertex_partition(int n, arr2d_var cliques); # todo: remove this function entirely and just focus on a single largest clique
    arr2d_fixed find_all_colorings(arr2d_fixed adj, int num_cols, arr2d_var partitions)
    arr2d_fixed reduce_colorings(int n, int num_colors, arr2d_fixed cols, arr2d_fixed isos, int num_threads)

cdef extern from "impl/fibering_single.c":
    legal_orbits_result graph_fiberings(arr2d_fixed adj, arr2d_var cliques, int min_cols, int max_cols, bint verbose, bint total_progress_bar, int num_threads, bint single_orbit)
    legal_orbits_result graph_fiberings_single_coloring(arr2d_fixed adj, int* coloring);

cdef extern from "impl/fibering_multi.c":
    void read_adj_matrix_graph6(char* geng, int* adj)
    char* graph6_from_adj_matrix(arr2d_fixed adj);
    ctypedef struct fibering_scheduler:
        pass
    ctypedef struct stream_result:
        arr2d_fixed results
        int num_checked
        int num_fiber
    fibering_scheduler make_scheduler(int n, int num_queues, int capacity_per_queue, int threads_per_queue, char* results_file_str)
    void add_to_scheduler(fibering_scheduler scheduler, arr2d_fixed adj, arr2d_var cliques)
    stream_result scheduler_finish(fibering_scheduler scheduler)

cdef extern from "impl/graph.c":
    bint graph_can_fiber(arr2d_fixed adj)
    bint is_graph_hyperbolic(arr2d_fixed adj)
    arr2d_fixed get_isometries(arr2d_fixed adj)

cdef extern from "impl/legal.c":
    ctypedef struct legal_orbits_result:
        arr2d_fixed colorings
        arr2d_var states


#### ARRAY CONVERSION ####

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