######## Useful functionality ########
## TODO

def coloring_has_legal_orbit():
# todo
    return

def all_k_legal_orbits(g, min_cols=0, max_cols=0, verbose=True, total_progress_bar=True, num_threads=1):
# todo
    return

def curv2(g):
    return 1-g.order()/2+g.size()/4

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
def ten_graph(n=5):
    g = Graph()

    for i in range(n):
        j = (i+1)%n
        g.add_edge(2*i, 2*i+1)
        g.add_edge(2*i, 2*j)
        g.add_edge(2*i, 2*j+1)
        g.add_edge(2*i+1, 2*j)
        g.add_edge(2*i+1, 2*j+1)

    return g