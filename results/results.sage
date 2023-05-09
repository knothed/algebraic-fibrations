###
### Functionality for reading and analyzing lists of graphs.
###

def read_graphs(file):
    i = 0
    res = []
    with open(file) as f:
        for line in f:
            print(i)
            i += 1
            res.append(graph_from_graph6(line))
    return res

# for n in {7,8,9,10}, read all fibering graphs.
def read_all(n):
    return read_graphs(f"results/all{n}.txt")

# for n in {10,11,12}, read all hyperbolic fibering graphs.
def read_hyp(n):
    return read_graphs(f"results/hyp{n}.txt")

######## GRAPH ANALYSIS ########

def all_that_k_fiber(graphs, k):
    return list(filter(lambda a: graph_k_fibers(a, k=k, num_threads=3), graphs))

# Frequency analysis of graph properties.
def frequencies(graphs, func):
    res = {}
    for g in graphs:
        s = func(g)
        if s in res:
            res[s] += 1
        else:
            res[s] = 1
    return res

def clique_number_frequencies(graphs):
    return frequencies(graphs, lambda a: a.clique_number())

def isom_frequencies(graphs):
    return frequencies(graphs, num_isometries)

def edge_frequencies(graphs):
    return frequencies(graphs, lambda a: a.size())

# Return all graphs in `graphs` which are not trivial extensions of any of `other` or the singular graph.
def filter_nontrivial(graphs, other):
    z = Graph({0:[],1:[]})
    trivial = graphs_which_are_trivial_extensions(graphs, other+[z])
    return [g for g in graphs if g not in trivial]

# Return all graphs in `graphs` which are trivial extensions of any of `other`.
# Thereby, a trivial extension of G is the join of G with a complete graph (with at least 1 vertex).
def graphs_which_are_trivial_extensions(graphs, other):
    from sage.graphs.graph_generators import GraphGenerators
    def is_trivial_extension(g):
        for o in other:
            n = g.order() - o.order()
            if n > 0 and g.is_isomorphic(o.join(GraphGenerators.CompleteGraph(n))):
                return True
        return False
    return list(filter(is_trivial_extension, graphs))

# The complement graph of a graph.
def dual(g):
    n = g.order()
    mat = ones_matrix(n) - identity_matrix(n) - g.am()
    return Graph(mat, format='adjacency_matrix')