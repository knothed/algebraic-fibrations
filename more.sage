###
### Additional functionality that is used in the masters thesis.
###

######## K-LEGALITY CHECKING ########

# Classify the given graphs by the maximum k for which they k-fiber.
def classify_max_k_fiberings(graphs, num_threads=1, verbose=True):
    result = {}
    for g in graphs:
        k = max_k_fibers(g, num_threads, verbose)
        if k in result:
            result[k].append(g)
        else:
            result[k] = [g]
    return result

# All orbits with are k-legal. Precondition: legal_orbits are all legal.
def k_legal_orbits(graph, k, legal_orbits):
    result = []
    for pair in legal_orbits:
        coloring = pair['coloring']
        for state in pair['states']:
            if max_k_orbit(graph, k, coloring, state) >= k:
                result.append({'coloring': coloring, 'state': state})
    return result

# Determine the maximum k such that the graph k-fibers.
# Returns -1 if the graph does not fiber, and returns infinity if the graph fibers for all k.
def max_k_fibers(graph, num_threads=1, verbose=True):
    orbits = all_legal_orbits(graph, num_threads=num_threads, verbose=verbose)
    if not orbits: return -1
    return max_k_orbits(graph, orbits)

# Determine the maximum k for which there is an orbit in legal_orbits which k-fibers.
def max_k_orbits(graph, legal_orbits):
    if curv3(graph) > 0: return 0
    k = 0
    for pair in legal_orbits:
        coloring = pair['coloring']
        for state in pair['states']:
            k = max(k, max_k_orbit(graph, k+1, coloring, state))
            if k == infinity: return k
    return k

# The maximum k for which this orbit k-fibers.
# Preconidition: orbit is 0-legal.
# If k drops below `min_k`, return 0.
def max_k_orbit(graph, min_k, coloring, state):
    k = infinity
    orbit = orbit_of(state, coloring)
    for verts in map(vertices_in_state, orbit):
        g = Graph(graph.am()[verts,verts], format='adjacency_matrix')
        if g.is_tree(): continue

        cc = g.clique_complex()
        if cc.fundamental_group().simplified().ngens() > 0: # nontrivial pi_1
            return 0

        if k <= 1: continue
        hom = cc.homology()
        m = max(hom.keys())
        for j in [2..min(m, k)]:
            if hom[j].ngens() > 0:
                k = j-1
                if k < min_k: return 0
                break
    return k


# For a graph to have a legal orbit, curv2 must be >= 0. (see JNW paper)
def curv2(g):
    return 1-g.order()/2+g.size()/4

# For a graph to have a 1-legal orbit, curv3 must be <= 0. (similar argument)
def curv3(g):
    return curv2(g)-g.triangles_count()/8



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
# This method changes the graph by modifying the existing vertex v and adding a new vertex having the highest label.
# Return the modified graph.
def inoue_edge_surgery(graph, v, v1, v2):
    cycle = get_incident_cycle(graph, v)
    if cycle is None: return None
    if len(cycle) < 6:
        print(f"Error: {v} has valence < 6")
        return None

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

    return graph

# Combine two löbell graphs into one by removing one same-valence vertex per graph and gluing its neighbors together along an antiprism.
# This method returns the new graph. It has v1+v2-2 vertices.
def inoue_antiprism(g1, g2, v1, v2):
    c1 = get_incident_cycle(g1, v1)
    c2 = get_incident_cycle(g2, v2)
    if c1 is None or c2 is None: return None
    if len(c1) != len(c2):
        print(f"Error: cycles around vertices have different lengths: {len(c1)} != {len(c2)}")
        return None

    # Create g1+g2 - (v1+v2)
    n1 = g1.num_verts()
    graph = Graph()
    graph.add_edges(g1.edges(sort=False))
    graph.delete_vertex(v1)
    graph.add_edges(map(lambda e: (e[0]+n1, e[1]+n1, e[2]), g2.edges(sort=False)))
    graph.delete_vertex(v2+n1)

    # Add antiprism
    n = len(c1)
    for i in range(n):
        graph.add_edge(c1[i],c2[i]+n1)
        graph.add_edge(c1[i],c2[(i+1)%n]+n1)

    graph.relabel()
    return graph

# Get the vertices incident to some vertex, sorted as a cycle. If they do not form a cycle, print an error and return None.
def get_incident_cycle(graph, v):
    incident = graph.edges_incident(v)
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

    return cycle

# All vertices with valence 6 or more. Use get_incidence_cycle to see which vertices to use in an inoue_edge_surgery call.
def valence_geq_6_verts(g):
    return list(filter(lambda v: len(g.edges_incident(v,sort=False)) > 5, g.vertices(sort=False)))


# The first nontrivial hyperbolic graph which fibers. It has 10 vertices.
def ten_graph():
    return wheel([2,2,2,2,2])

def prism(g, anti=False):
    p = Graph()
    n = g.order()
    p.add_edges(g.edges(sort=False))
    p.add_edges(map(lambda a: (a[0]+n,a[1]+n),g.edges(sort=False)))
    for i in range(n):
        p.add_edge(i,i+n)
        if anti:
            p.add_edge((i+1)%n,i+n)
    return p

def wheel(ns):
    g = Graph()
    ks = list(map(graphs.CompleteGraph, ns))
    sums = [sum(ns[:i]) for i, x in enumerate(ns)]

    for i in range(len(ns)):
        g.add_edges(map(lambda e: (e[0]+sums[i],e[1]+sums[i]), ks[i].edges(sort=False)))
        j = (i+1)%len(ns)
        for a in range(ns[i]):
            for b in range(ns[j]):
                g.add_edge(a+sums[i], b+sums[j])

    return g

def cycle(n):
    return graphs.CycleGraph(n)

def empty(n):
    g = Graph()
    for i in range(n):
        g.add_vertex(i)
    return g

def complete(n):
    return graphs.CompleteGraph(n)

def join(gs):
    if not gs: return graphs.EmptyGraph()
    return gs[0].join(join(gs[1:]))