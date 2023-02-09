from sage.graphs.graph_coloring import *

def check_all_graphs_c(n,add=2,verbose=True,hyper=False):
    i = 0
    for g in graphs(n):
        i += 1
        n = time_ms()
        if verbose: print(f'{i}:')
        if not precheck(g): continue
        if hyper and not is_hyperbolic(g): continue
        if graph_fiberings(g, max_cols=g.chromatic_number()+add-1, verbose=verbose): yield g
        if verbose: print(f'{i} took {time_ms()-n} ms')

def check_all_graphs(n,add=2,verbose=True):
    i = 0
    for g in graphs(n):
        i += 1
        if verbose: print(f'{i}:')
        if not precheck(g): continue
        if graph_fibers(g,add,verbose): yield g

def precheck(g):
    if curv2(g) < 0: return False
    if g.order() > 2 and not g.is_connected(): return False
    if g.order() > 3 and g.is_tree(): return False
    # etc.
    return True

def graph_fibers(g,col_tries=2,verbose=True):
    c = g.chromatic_number()
    n = g.order()

    import time
    def now():
        return time.time_ns() // 1000000

    t1 = now()
    legals = legal_states(g.adjacency_matrix())
    t2 = now()
    if verbose: print(f'get legal states: {t2-t1}ms')

    for c in range(c,c+col_tries):
        t3 = now()
        cols = all_colorings(g,c)
        t4 = now()
        if verbose: print(f'get {len(cols)} {c}-colorings: {t4-t3}ms')

        total = len(cols)
        i = 0
        for col in cols:
            legal = legals.copy()
            while legal: #todo: length check, at least 2^(c-1)
                state = next(iter(legal))
                orbit = get_orbit_nat(n,col,c,state)
                if orbit <= legal:
                    return col, get_binary(n,state)
                legal -= orbit
            i += 1
        if verbose: print(f'checking the {c}-colorings: {now()-t4}ms')

# the (half-)orbit of a coloring on a state which is represented as a number in 0 ..< 2^(n-1).
# vertex n is never in the orbit.
def get_orbit_nat(n,col,num_colours,i):
    orbit = set()
    actors = list_of_binaries(num_colours)
    redundant_col = col[n-1] # the color of vertex n

    for actor in actors: # all 2^(c-1) subsets of 0..<col (modulo negation) acting on the state i
        if actor[redundant_col] == 1: continue
        new_state = get_binary(n,i)
        for d in range(n): # swap every bit of i which is present in actor
            if actor[col[d]] == 1: # in C: use bitwise xor
                new_state[d] = 1 - new_state[d]
        orbit.add(get_integer(new_state))

    return orbit

######## CONCRETE GRAPHS ########

# Create a löbell graph.
# layer 0: 2 outer vertices. layers 1+2: k+k inner vertices.
def löbell_graph(k):
    g = Graph()

    # convert from tuple (layer,vertex) to int
    def c(l,v):
        return v + [0,2,k+2][l]

    # layer 0 inward
    for i in range(k):
        g.add_edge(c(0,0), c(1,i))
        g.add_edge(c(0,1), c(2,i))
    # inside layers 1 and 2
        g.add_edge(c(1,i), c(1,(i+1)%k))
        g.add_edge(c(2,i), c(2,(i+1)%k))
    # between layer 1 and 2
    for i in range(k):
        g.add_edge(c(1,i), c(2,i))
        g.add_edge(c(1,i), c(2,(i+1)%k))

    return g

# Extend the löbell graph by doing "inverse" edge surgery to the hyperbolic polygon defined by it.
# This means, when g's faces are all triangles, v is a vertex of valence >= 6, v is connected to valence >=5 vertices,
# in particular v1 and v2, and v1-v-v2 isn't part of a 4-cycle, we replace v by a new edge and connect the two now vertices with v's neighbors
# such that only v1 and v2 are connected to both, the others to one of them.
# This method changes the graph by modifying the existing vertex v and adding a new vertex whose label is returned.
def add_inoue_edges(graph, v, v1, v2):
    # Get vertices incident to v
    incident = graph.edges_incident(v)
    if len(incident) < 6:
        print(f"Error: {v} has valence < 6")
        return None
    for i in range(len(incident)):
        incident[i] = incident[i][0] if incident[i][1] == v else incident[i][1]

    # Sort vertices so they form a cycle
    cycle = [incident[0]]
    incident.remove(cycle[0])
    while incident:
        neighbors = list(filter(lambda x: graph.has_edge(x, cycle[-1]), incident))
        if not neighbors:
            print(f"Error: {v} is not center of a valid cycle")
            return None
        cycle.append(neighbors[0])
        incident.remove(neighbors[0])

    # Get positions of v1 and v2 in the cycle
    try:
        idx1 = cycle.index(v1)
        idx2 = cycle.index(v2)
        idx_low = min(idx1, idx2)
        idx_hi = max(idx1, idx2)
    except:
        print(f"Error: both vertices {v1} and {v2} must be incident to {v}")
        return None

    n = len(cycle)
    if (idx1 - idx2) % n < 3 or (idx2 - idx1) % n < 3:
        print(f"Error: vertices {v1} and {v2} are too close: they must not form a 3- or 4-cycle with {v}")
        return None

    # Add new vertices to the graph
    graph.delete_vertex(v)
    graph.add_vertex(v)
    new1 = v
    new2 = graph.add_vertex()
    graph.add_edges([(new1,new2), (new1,v1), (new1,v2), (new2,v1), (new2,v2)])
    for i in range(idx_low+1, idx_hi):
        graph.add_edge(new1, cycle[i])
    for i in range(idx_hi+1, idx_low+n):
        graph.add_edge(new2, cycle[i%n])

    return new2


######## GRAPH PROPERTIES ########

def is_hyperbolic_old(g):
    cycles = g.to_directed().all_simple_cycles(max_length=4)
    tris = list(map(set, filter(lambda c: len(c)==4, cycles)))
    quads = map(set, filter(lambda c: len(c)==5, cycles))
    for quad in quads:
        for tri in tris:
            if tri <= quad: break
        else:
            return False
    return True

def curv2(g):
    return 1-g.order()/2+g.size()/4

def curv3(g):
    return curv2(g)-g.triangles_count()/8


######## COLORINGS UP TO ISOMORPHISM AND CONTAINMENT ########

# 1. color label permutation: {0: a, 1: b, ...} ~ {0: b, 1: a, ...}
# 2. colorings come from graph isomorphism: {0: a, 1: b} ~ {0: f(a), 1: f(b)}
# 3. coloring is excessive, i.e. finer than another: {0: a, 1: b, 2: c} < {0: a, 1: b+c}


# all colorings with n colors up to color permutation and graph isomorphism.
def all_colorings(g,n):
    from sage.graphs.graph_coloring import all_graph_colorings
    cols = list(all_graph_colorings(g,n,vertex_color_dict=True))
    isos = get_isometries(g.adjacency_matrix())
    return kill_permutations_and_isos(cols, isos)

def kill_permutations_and_isos(cols, isos):
    result = []

    def is_new(col):
        c2i = 0
        for c2 in result:
            fi = 0
            for f in isos:
                if is_color_permutation_iso(col,c2,f):
                    return False
                fi += 1
            c2i += 1
        return True

    i = 0
    l = len(cols)
    for col in cols:
        #print(f'{i} of {l}')
        if is_new(col): result.append(col)
        i += 1

    return result

# check whether c2 is some color permutation of c1 under the automorphism f.
def is_color_permutation_iso(c1, c2, f):
    swaps = {}
    for v,col2 in c2.items():
        col1 = c1[f[v]]
        if not col2 in swaps:
            swaps[col2] = col1
        elif swaps[col2] != col1:
            return False
    return True


######## LEGALITY TESTING ########

# find all states of g where there links are legal.
# (only consider states where the highest vertex is excluded; this suffices because of I/O symmetry)
def legal_states(m):
    result = set()
    n = m.dimensions()[0]
    for i in range(2^(n-1)):
        state = get_binary(n-1,i)
        links = adlinks2(m,state)
        for g in links:
            cc = g.clique_complex()
            if g.order() == 0 or not cc.is_connected(): break
        else:
            result.add(i)
    return result

#1-legality:
#for i in range(1000):
#     print(i)
#     g=next(gs)
#     if g.is_tree(): continue
#     if g.is_connected() and len(g.clique_complex().fundamental_group().gens())==0: a += 1


# takes states where the highest vertex is excluded, i.e. OUT.
def adlinks2(G, state):
    n = G.dimensions()[0]
    # case = 0 is ascending, case = 1 is descending
    links = []
    for case in range(2):
        columns = []
        for i in range(n-1):
            if state[i] == case:
                columns.append(i)
        if case == 0: columns.append(n-1)
        G_reduced = G[columns, columns]
        links.append(Graph(G_reduced, format='adjacency_matrix'))
    return links