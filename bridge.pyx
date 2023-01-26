cdef extern from "fast.c":
    int* kill_permutations_and_isos(int n, int num_cols,
                                    int* cols, int cols_count,
                                    int* isos, int isos_count,
                                    int* result_count)

    bint is_state_legal(int n, int* adj_matrix, int state)

    int* get_isometries(int n, int* adj, int* result_count)

from sage.all import *
# from cpython cimport array
# import array

from libc.stdlib cimport malloc, free

def reduced_colorings(g,c):
    from sage.graphs.graph_coloring import all_graph_colorings
    pycols = list(all_graph_colorings(g,c,vertex_color_dict=True))
    pyisos = get_isometries_(g.adjacency_matrix())

    n = g.order()
    cdef int* ccols = list_of_dicts_to_array(pycols,n)
    cdef int* cisos = list_of_lists_to_array(pyisos,n)

    cdef int result_count = 0
    cdef int* cresult = kill_permutations_and_isos(n,c,ccols,len(pycols),cisos,len(pyisos),&result_count)
    pyresult = list_of_lists_from_array(cresult, result_count, n)

    free(ccols)
    free(cisos)
    return pyresult

def state_legal(adj, int state):
    cdef int n = adj.dimensions()[0]
    cdef int* cadj = list_to_array(adj._list())
    cdef int legal = is_state_legal(n,cadj,state)
    free(cadj)
    return legal == 1

def get_isometries_c(adj):
    cdef int n = adj.dimensions()[0]
    cdef int* cadj = list_to_array(adj._list())

    cdef int result_count = 0
    cdef int* cisos = get_isometries(n,cadj,&result_count)
    pyisos = list_of_lists_from_array(cisos, result_count, n)

    free(cadj)
    return pyisos

## C CONVERSION TOOLS ##

cdef int* list_to_array(list):
    cdef int l = len(list)
    cdef int* res = <int*> malloc(l*sizeof(int))
    for i in range(l):
        res[i] = list[i]
    return res

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


## FROM MARTELLI, TODO: REWRITE IN CYTHON

def get_isometries_(G, level = 0, isometries = [], isometry = [], n = 0):
    if n == 0:
        n = G.dimensions()[0]
        isometry = [0 for i in range(n)]
        isometries = []
    if level == n:
        isometries.append(isometry[:])
        # print (isometry)    # ACHTUNG
        # print (len(isometries)) # ACHTUNG
    for i in range(n):
        if isometry[:level].count(i) == 0:
            is_ok = True
            for j in range(level):
                if G[level, j] != G[i, isometry[j]]:
                    is_ok = False
                    break
            if is_ok == True:
                isometry[level] = i
                get_isometries_(G, level + 1, isometries, isometry, n)
    if level == 0:
        return isometries