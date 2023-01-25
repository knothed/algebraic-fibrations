# code from: https://people.dm.unipi.it/martelli/files/Diagonal.sage

from copy import deepcopy

# Routines in Sage to play the JNW game with some right-angled polytope

# The right-angled polytope P is assumed to lie in some constant curvature space and have finite volume.
# It is encoded as a matrix G with entries 0 and 1 that represents the adjacence graph between facets.
# Often the parameter d = dim P is also needed (although in principle d could be deduced from G by studying subgraphs
# that represent n-simplexes and n-cubes). So P is often encoded as a pair (G, d).

# The number of facets of P is also the size of G and is usually indicated with n.

#######

# TECHNICAL GENERAL ROUTINES:

# Returns a list of length n with entries in 0, 1 (called a "binary" below)
# that represents the number i as a binary number. Example: (4, 2) gives [0, 1, 0, 0, 0]

def get_binary(n, i):
    answer = []
    i2 = i
    for j in range(n):
        if i2 % 2:
            answer.append(1)
            i2 = i2-1
        else:
            answer.append(0)
        i2 = i2/2
    return answer

# Returns the integer represented by the given binary. Example: [0, 1, 1, 0, 0, 0, 0] gives 6

def get_integer(b):
    answer = 0
    for i in range(len(b)):
        if b[i] == 1:
            answer = answer + 2^i
    return answer

# Returns the list of 2^n binaries with length n and entries 0, 1

def list_of_binaries(n):
    if n == 1:
        return [[0], [1]]
    sublist = list_of_binaries(n-1)
    answer = []
    for i in range(2):
        for binary in sublist:
            new_binary = [i]
            new_binary.extend(binary[:])
            answer.append(new_binary[:])
    return answer

######

# QUATERNIONS AND OCTONIONS:

# Quaternion multiplication

def quat_mult(z, w):
    (a,b,c,d,e,f,g,h) = (z[0], z[1], z[2], z[3], w[0], w[1], w[2], w[3])
    A = a*e - b*f - c*g - d*h
    B = a*f + b*e + c*h - d*g
    C = a*g + c*e + d*f - b*h
    D = a*h + d*e + b*g - c*f
    return [A,B,C,D]

# Returns the 7 oriented lines in the Fano plane
# If expanded == True, returns each [a, b, c] also as [b, c, a] and [c, a, b]

def Fano_lines(expanded = False):
    lines = [[1, 2, 4], [2, 3, 5], [3, 4, 6], [4, 5, 7], [5, 6, 1], [6, 7, 2], [7, 1, 3]]
    if expanded == True:
        expanded_lines = []
        for line in lines:
            for i in range(3):
                new_line = line[:]
                for j in range(3):
                    new_line[j] = line[(i+j) % 3]
                expanded_lines.append(new_line)
        lines = expanded_lines[:]
    return lines

# Octonion multiplication

def octo_mult(v, w):
    lines = Fano_lines(True)

    answer = [0 for i in range(8)]
    for i in range(8):
        for j in range(8):
            if i == 0:
                answer[j] = answer[j] + v[i] * w[j]
            if i > 0 and j == 0:
                answer[i] = answer[i] + v[i] * w[j]
            if i > 0 and i == j:
                answer[0] = answer[0] - v[i] * w[j]
            if i > 0 and j > 0 and i != j:
                for line in lines:
                    if line[0] == i and line[1] == j:
                        k = line[2]
                        answer[k] = answer[k] + v[i] * w[j]
                    if line[0] == j and line[1] == i:
                        k = line[2]
                        answer[k] = answer[k] - v[i] * w[j]
    return answer

######

# IMPORTANT POLYTOPES:

# Returns G for the n-cube

def cube(n):
    G = matrix(2*n, 2*n)
    for f in range(2*n):
        for g in range(2*n):
            if g != f and g != f+n and f != g+n:
                G[f,g] = 1
    return G

# Returns G for the ideal regular octahedron

def octa():
    G = matrix(8,8)
    for f in range(8):
        for g in range(8):
            if (f <= 3 and g <= 3) or (f >= 4 and g >= 4):
                if (g-f) % 4 == 1 or (g-f) % 4 == 3:
                    G[f,g] = 1
            if abs (g-f) == 4: G[f,g] = 1
    return G

# Returns G for P3

def P3():
    G = Matrix(6, 6)
    for i in range(3):
        for j in range(3):
            if i!=j:
                G[i,j] = 1
                G[i+3,j+3] = 1
            if i==j:
                G[i,j+3] = 1
                G[i+3,j] = 1
    return G

# Returns G for P4

def P4():
    G = Matrix(10,10)
    for i in range(5):
        for j in range(5):
            if (i-j) % 5 in [1,4]:
                G[i,j] = 1
            if (i-j) % 5 in [2,3]:
                G[i+5, j+5] = 1
            if i!=j:
                G[i,j+5] = 1
                G[i+5,j] = 1
    return G

# Returns the 8 elements of the quaternion group

def Q8():
    units = [[1, 0, 0, 0], [-1, 0, 0, 0], [0, 1, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, -1, 0], [0, 0, 0, 1], [0, 0, 0, -1]]
    return units

# Returns the binary tetrahedral group, that is the facets of the ideal 24-cell

def binary_tetrahedral():
    elements = [[1, 0, 0, 0], [1/2, 1/2, 1/2, 1/2], [-1/2, 1/2, 1/2, 1/2]]
    units = Q8()
    answer = []
    for q in elements:
        for p in units:
            new_element = quat_mult(p,q)
            answer.append(new_element[:])
    return answer

# Returns the 24 vertices in R4 representing the facets of the 24-cell
# These are (pm 1, 0, 0, pm 1) and permutations

def C24_facets():
    facets = []
    for i in range(-1,2):
        for j in range(-1,2):
            for h in range(-1,2):
                for k in range(-1,2):
                    new_facet = [0,0,0,0]
                    new_facet[0] = i
                    new_facet[1] = j
                    new_facet[2] = h
                    new_facet[3] = k
                    if new_facet.count(0) == 2:
                        facets.append(new_facet[:])
    return facets

# Returns G for the ideal regular 24-cell

def C24():
    G = Matrix(24,24)
    facets = C24_facets()
    for f in range(24):
        fc = facets[f]
        for g in range(24):
            gc = facets[g]
            dsquared = (fc[0] - gc[0])^2 + (fc[1] - gc[1])^2 + (fc[2] - gc[2])^2 + (fc[3] - gc[3])^2
            if dsquared == 2:
                G[f,g] = 1
    return G

# Returns the binary icosahedral group, that is the facets of the 120-cell

def binary_icosahedral():
    tetra = binary_tetrahedral()
    phi = (1+sqrt(5))/2
    elements = [[1, 0, 0, 0], [0, 1/2, 1/(2*phi), phi/2], [0, 1/2, 1/(2*phi), -phi/2], [0, 1/2, -1/(2*phi), phi/2], [0, 1/2, -1/(2*phi), -phi/2]]
    answer = []
    valori_possibili = [0, 1, -1, 1/2, -1/2, phi/2, -phi/2, 1/(2*phi), -1/(2*phi)]
    for q in elements:
        for p in tetra:
            new_element = quat_mult(p,q)
            for i in range(4):
                for v in valori_possibili:
                    if abs(new_element[i] - v) < 0.1:
                        new_element[i] = v
            answer.append(new_element[:])
    return answer

# Returns G for the right-angled 120-cell

def C120():
    G = Matrix(120,120)
    ico = binary_icosahedral()
    phi = (1+sqrt(5))/2
    x = phi/2
    for i in range(120):
        v = ico[i]
        for j in range(i):
            w = ico[j]
            sp = sum([v[k] * w[k] for k in range(4)])
            if abs(sp - x) < 0.1:
                G[i,j] = 1
                G[j,i] = 1
    return G

# Returns G for P5

def P5():
    G = matrix(16,16)
    for i in range(8):
        for j in range(8):
            if (i-j) % 8 in [1,2,6,7]:
                G[i,j] = 1
            if (i-j) % 8 in [0,1,2,4,6,7]:
                G[i,j+8] = 1
                G[j+8,i] = 1
            if (i-j) % 8 in [2,3,5,6]:
                G[i+8,j+8] = 1
    return G

# Returns G for P6

def P6():
    G = matrix(27,27)
    basic_vertices = [[-1,0,0,0,0,0,0], [1,1,0,0,0,0,1], [1,0,1,1,1,1,2]]
    vertices = []
    for i in range(6):
        for j in range(3):
            new_vertex = [basic_vertices[j][(k+i) % 6] for k in range(6)]
            new_vertex.append(basic_vertices[j][6])
            vertices.append(new_vertex)
    more_vertices = [[1,0,1,0,0,0,1], [0,1,0,0,1,0,1], [0,0,0,1,0,1,1], [0,1,0,1,0,0,1], [0,0,1,0,0,1,1], [1,0,0,0,1,0,1], [0,0,1,0,1,0,1], [1,0,0,1,0,0,1], [0,1,0,0,0,1,1]]
    vertices.extend(more_vertices)
    for i in range(27):
        for j in range(27):
            somma = 0
            for k in range(6):
                somma = somma + vertices[i][k]*vertices[j][k]
            somma = somma - vertices[i][6]*vertices[j][6]
            if somma == 0:
                G[i,j] = 1
            else:
                G[i,j] = 0
    return G

def P6_old():
    G = matrix(27,27)
    for i in range(9):
        for j in range(9):
            if (i-j) % 9 in [1,2,7,8]:
                G[i,j] = 1
            if (i-j) % 9 in [2,3,4,5,6,7]:
                G[i,j+9] = 1
                G[j+9,i] = 1
            if (i-j) % 9 in [1,3,4,5,6,8]:
                G[i,j+18] = 1
                G[j+18,i] = 1
            if (i-j) % 9 in [1,4,5,8]:
                G[i+9, j+9] = 1
            if (i-j) % 9 in [1,2,3,6,7,8]:
                G[i+9, j+18] = 1
                G[j+18, i+9] = 1
            if (i-j) % 9 in [2,4,5,7]:
                G[i+18, j+18] = 1
    return G

# Returns the 240 vertices in R^8 representing the facets of P8

def P8_facets():
    answer = []
    for i in range(8):
        v = [0 for j in range(8)]
        v[i] = 1
        answer.append(v[:])
        v[i] = -1
        answer.append(v[:])
    lines = Fano_lines()
    binaries = list_of_binaries(4)
    for line in lines:
        v = [0 for j in range(8)]
        v_complement = [0 for j in range(8)]
        v[0] = 1/2
        for i in range(1,8):
            if i in line:
                v[i] = 1/2
            else:
                v_complement[i] = 1/2
        for binary in binaries:
            w = v[:]
            w_complement = v_complement[:]
            count_w = 0
            count_w_complement = 0
            for i in range(8):
                if w[i] != 0:
                    if binary[count_w] == 1:
                        w[i] = - w[i]
                    count_w = count_w + 1
                if w_complement[i] != 0:
                    if binary[count_w_complement] == 1:
                        w_complement[i] = - w_complement[i]
                    count_w_complement = count_w_complement + 1
            answer.append(w[:])
            answer.append(w_complement[:])
    return answer

# Returns G for P8

def P8():
    G = matrix(240,240)
    facets = P8_facets()
    for i in range(240):
        for j in range(i+1, 240):
            f = facets[i]
            g = facets[j]
            sp = sum([f[k] * g[k] for k in range(8)])
            if sp == 1/2:
                G[i,j] = 1
                G[j,i] = 1
    return G

# Returns G for P7

def P7():
    G8 = P8()
    G = get_facet(G8, 0)
    return G

#####

# ROUTINES THAT COLLECT GEOMETRIC INFORMATION FROM G:

# Returns the ordered list of facets adjacent to f
# If with_f == True it inserts f in the list

def adjacent_facets(G, f, with_f = False):
    n = G.dimensions()[0]
    answer = []
    for i in range(n):
        if G[f,i]:
            answer.append(i)
        if with_f and i == f:
            answer.append(i)
    return answer

# Returns the facet f as a new incident matrix

def get_facet(G, f):
    af = adjacent_facets(G, f)
    new_n = len(af)
    new_G = matrix(new_n, new_n)
    for i in range(new_n):
        for j in range(new_n):
            if G[af[i], af[j]] == 1:
                new_G[i,j] = 1
    return new_G

# Returns the (d-1)-cubic cusps of P. Here d is the dimension of P.
# A cusp of P is a list of 2(d-1) facets of P, those that are adjacent to the cusp.
# The list is produced so that i and i + (d-1) are opposite facets in the cube

def cusps(G, d, answer = [], level = 0, facets = [], n = 0):
    if n == 0:
        n = G.dimensions()[0]
        answer = []
        facets = []
    if facets == []:
        for i in range(2*(d-1)):
            facets.append(0)
    if level < 2*(d-1):
        start = 0
        if level > 0: start = max(facets[:level]) + 1
        for f in range(start, n):
            sum = 0
            for j in range(level):
                if G[f, facets[j]]: sum = sum + 1
            if sum >= level - 1:
                facets[level] = f
                cusps(G, d, answer, level + 1, facets, n)
        if level == 0:
            return answer
    else:
        is_cube = True
        for i in range(2*(d-1)):
            sum = 0
            for j in range(2*(d-1)):
                if G[facets[i], facets[j]]: sum = sum + 1
            if sum != 2*(d-2): is_cube = False
        if is_cube == True:

            # We reorder the faces

            ordered_facets = facets[:]
            for i in range(d-1):
                for j in range(i+1, 2*(d-1)):
                    if G[ordered_facets[i], ordered_facets[j]] == 0 and j != i + d-1:
                        ordered_facets[j], ordered_facets[i + d-1] = ordered_facets[i + d-1], ordered_facets[j]

            answer.append(ordered_facets[:])

# Returns a complete list of all the isometries of P

def get_isometries(G, level = 0, isometries = [], isometry = [], n = 0):
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
                get_isometries(G, level + 1, isometries, isometry, n)
    if level == 0:
        return isometries

# Returns the maximum number of pairwise disjoint facets

def max_disjoint_facets(G, answer = 0, facets = [], n = 0):
    a = -1
    if (n == 0):
        n = G.dimensions()[0]
    else:
        a = facets[-1]
    if len(facets) <= 3:
        print(facets)
    for i in range(a+1, n):
        ammissibile = True
        for f in facets:
            if G[f,i] == 1:
                ammissibile = False
                break
        if ammissibile == True:
            facets.append(i)
            if answer < len(facets):
                answer = len(facets)
            if answer == len(facets):
                print(facets)
            answer = max_disjoint_facets(G, answer, facets, n)
            facets.pop()
    return answer

#####

# COLOURS:

# A colouring for P with c colours is a list of n elements in {0, ..., c-1}, where n is the number of facets of P.
# A fixed colouring on P generates a manifold M by mirroring.

# Returns the induced colours on a facet f of P
# If renormalize == True, it fills the gaps transforming [0,2,2,3] into [0,1,1,2]

def get_colour_facet(G, f, colouring, renormalize = True):
    af = adjacent_facets(G, f)
    n = len(af)
    new_colouring = [colouring[af[i]] for i in range(n)]
    if renormalize == True:
        for i in range(max(new_colouring)):
            while i not in new_colouring and i < max(new_colouring):
                for j in range(n):
                    if new_colouring[j]>i: new_colouring[j] = new_colouring[j]-1
                # print (new_colouring, i)
    return new_colouring

# Returns the number of colours on all the cusps (= ideal vertices) of P

def get_colour_cusps(G, d, colouring):
    cuspidi = cusps(G, d)
    answer = []
    for cusp in cuspidi:
        colours = [colouring[cusp[i]] for i in range(2*(d-1))]
        num_colours = len(set(colours))
        answer.append(num_colours)
    return answer

# Returns the isometries that preserve a colouring (as a partition, not as individual colours)
# The routine is not very efficient, since it first constructs all the isometries of P.

def get_colouring_preserving_isometries(G, colouring):
    n = G.dimensions()[0]
    full_isometries = get_isometries(G)
    answer = []
    for isometry in full_isometries:
        preserves_colouring = True
        for f in range(n):
            for g in range(n):
                if colouring[f] == colouring[g] and colouring[isometry[f]] != colouring[isometry[g]]:
                    preserves_colouring = False
                    break
            if preserves_colouring == False:
                break
        if preserves_colouring == True:
            answer.append(isometry[:])
    return answer

# Returns the homology of the manifold M determined by G and its colouring

def get_homology(G, colouring):
    n = G.dimensions()[0]
    c = max(colouring) + 1
    binaries = list_of_binaries(c)
    betti = [1]

    graphs = []
    counter = 0
    found = 0
    for binary in binaries:
        columns = []
        for i in range(n):
            if binary[colouring[i]] == 1:
                columns.append(i)
        # print (columns)
        G_reduced = G[columns, columns]
        new_graph = Graph(G_reduced, format='adjacency_matrix')
        is_new = True
        for graph in graphs:
            if graph[0].is_isomorphic(new_graph):
                graph[1] = graph[1] + 1
                is_new = False
        if is_new == True:
            new_graph_with_mult = [new_graph, 0]
            new_graph_with_mult[1] = 1
            graphs.append(new_graph_with_mult)
            found = found + 1
        counter = counter + 1
        print (counter, found)
    print ("----")
    counter = 0
    for g in graphs:
        cc = g[0].clique_complex()
        for i in range(cc.dimension()+2):
            # print (i, betti)
            while len(betti) <= i: betti.append(0)
            betti[i] = betti[i] + cc.betti(i-1) * g[1]
        print (counter, betti)
        counter = counter + 1
    return betti

# Returns the 2-colouring of the ideal regular octahedron:

def octa_colouring():
    colouring = [0, 1, 0, 1, 1, 0, 1, 0]
    return colouring

# Returns the 3-colouring of P3:

def P3_colouring():
    colouring = [0, 1, 2, 1, 2, 0]
    return colouring

# Returns the 5-colouring of P4:

def P4_colouring():
    colouring = [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
    return colouring

# Returns the 3-colouring of the ideal 24-cell:

def C24_colouring():
    answer = []
    facets = C24_facets()
    for f in facets:
        new_colour = 0
        if f[0]*f[3] != 0 or f[1]*f[2] != 0: new_colour = 0
        if f[0]*f[2] != 0 or f[1]*f[3] != 0: new_colour = 1
        if f[0]*f[1] != 0 or f[2]*f[3] != 0: new_colour = 2
        answer.append(new_colour)
    return answer

# Returns the 5-colouring of the 120-cell:

def C120_colouring():
    answer = []
    for i in range(5):
        for j in range(24):
            answer.append(i)
    return answer

# Returns the 8-colouring of P5:

def P5_colouring():
    colouring = [0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7]
    return colouring

# Returns the 9-colouring of P6:

def P6_colouring():
    # old_colouring = [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8]
    colouring = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8]
    return colouring

# Returns the 15-colouring of P8:

def P8_colouring():
    answer = [0 for i in range(16)]
    facets = P8_facets()
    for j in range(7):
        for k in range(32):
            i = 16 + 32 * j + k
            colour = 2*j + 1
            parity = facets[i].count(-1/2) % 2
            if parity == 1:
                colour = colour + 1
            answer.append(colour)
    return answer

# Returns the 14-colouring of P7:

def P7_colouring():
    colouring = get_colour_facet(P8(), 0, P8_colouring())
    return colouring

#####

# STATES:

# A state for P is a list of n elements with entries in {0, 1}, where n is the number of facets of P.
# Each of the polytopes P3, ..., P8, octa, and 24-cell is equipped with a preferred "base state".

# Returns the induced state on a facet f of P:

def get_state_facet(G, f, state):
    af = adjacent_facets(G, f)
    n = len(af)
    new_state = [state[af[i]] for i in range(n)]
    return new_state

# Returns the base state of the ideal regular octahedron:

def octa_state():
    state = [1, 1, 1, 0, 1, 0, 0, 0]
    return state

# Returns the base state of P3:

def P3_state():
    state = [1, 1, 1, 0, 0, 0]
    return state

# Returns the base state of P4:

def P4_state():
    state = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
    return state

# Returns the base state of the right-angled ideal 24-cell:

def C24_state():
    state = [1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1]
    return state

# Returns the base state of the right-angled 120-cell:

def C120_state():
    base_state = [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]
    state = []
    for i in range(5):
        for j in range(24):
            state.append(base_state[j])
    return state

# Returns the base state of P5:

def P5_state():
    state = [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]
    return state

# Returns the base state of P6:

def P6_state():
    # state_old = [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    state = [1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]
    return state

# Returns the base state of P8:

def P8_state():
    facets = P8_facets()
    answer = [0 for i in range(240)]

    # Assigns the status 1 to \pm 1, \pm e_1, \pm e_2, \pm e_4:

    quaternions = [facets[i] for i in [0, 1, 2, 3, 4, 5, 8, 9]]

    for i in range(16):
        if facets[i] in quaternions:
            answer[i] = 1

    # In each hextet, assigns the status 1 to the 8 vertices obtained from a base element
    # by left-multiplying via \pm 1, \pm e_1, \pm e_2, \pm e_4

    for i in range(7):
        base = 16 + 32 * i
        for case in range(2):
            if case == 1:
                base = base + 2
            facet = facets[base]
            # print (facet) #ACHTUNG
            for j in range(8):
                factor = quaternions[j]
                v = octo_mult(factor, facet)
                pos = facets.index(v)
                answer[pos] = 1
    return answer

# Returns the base state of P7:

def P7_state():
    state = get_state_facet(P8(), 0, P8_state())
    return state

# Returns the orbit of a state, produced by the action of a given colouring.
# The output is a list of states.

def get_orbit(G, colouring, state):
    n = G.dimensions()[0]
    num_colours = max(colouring) + 1
    orbit = []
    binaries = list_of_binaries(num_colours)
    for binary in binaries:
        new_state = state[:]
        for i in range(n):
            if binary[colouring[i]] == 1:
                new_state[i] = 1 - new_state[i]
        orbit.append(new_state[:])
    return orbit

# Returns True if the given colouring and state induce a fibering on all the cusps

def all_cusps_fiber(G, d, colouring, state):
    cuspidi = cusps(G, d)
    answer = []
    for cusp in cuspidi:
        cusp_state = [state[cusp[i]] for i in range(2*(d-1))]
        cusp_colouring = [colouring[cusp[i]] for i in range(2*(d-1))]
        cusp_fibers = False
        for i in range(2*(d-1)):
            for j in range(i):
                if cusp_colouring[i] == cusp_colouring[j] and cusp_state[i] != cusp_state[j]:
                    cusp_fibers = True
        if cusp_fibers == False:
            return False
    return True

# Returns all the states of P, considered up to colour-preserving isometries of P
# and up to the action of the colouring
# Set the first k numbers to be zero (this can be often obtained by minimality, it speeds up the process)

def get_states(G, colouring, k=0):
    n = G.dimensions()[0]
    isometries = get_colouring_preserving_isometries(G, colouring)
    states = list_of_binaries(n-k)
    for state in states:
        for i in range(k):
            state.insert(0,0)
    print ("Isometries:", len(isometries))
    minimal_states = []

    terminal_output = open('/dev/stdout', 'w')
    for state in states:
        print (state, file=terminal_output)
        is_minimal = True
        for isometry in isometries:
            new_state = [state[isometry[i]] for i in range(n)]
            orbit = get_orbit(G, colouring, new_state)
            for st in orbit:
                if st < state:
                    print ("Not minimal", new_state, st, file=terminal_output)
                    is_minimal = False
                    break
            if is_minimal == False:
                break
        if is_minimal == True:
            minimal_states.append(state[:])
            print ("is minimal", file=terminal_output)
    return minimal_states


# ASCENDING AND DESCENDING LINKS
# Some routines to calculate and investigate the ascending and descending links

# Calculates the ascending and descending links of a state
# It returns a list of two graphs (the links are the clique complexes of these)

def adlinks(G, state):
    n = G.dimensions()[0]

    # case = 0 is ascending, case = 1 is descending

    links = []
    for case in range(2):
        columns = []
        for i in range(n):
            if state[i] == case:
                columns.append(i)
        G_reduced = G[columns, columns]
        links.append(Graph(G_reduced, format='adjacency_matrix'))
    return links

# Given a polyhedron P with a colouring and a state,
# it classifies all the ascending and descending links
# and returns them as graphs (the links are the clique complexes of these).

def get_links_graphs(G, colouring, state):
    n = G.dimensions()[0]
    orbit = get_orbit(G, colouring, state)
    print ("Number of states:", len(orbit))

    # The list 'graphs' will contain one representative in each isomorphism class
    # Each element of the list is a triple [graph, n0, n1] where 'graph' is the graph and
    # n0, n1 is the number of occurrences of it among the ascending and descending links, respectively

    graphs = []
    counter = 0
    found = 0
    for o in orbit:
        links = adlinks(G, o)

        # case = 0 is ascending, case = 1 is descending

        for case in range(2):
            new_graph = links[case]
            is_new = True
            for graph in graphs:
                if graph[0].is_isomorphic(new_graph):
                    graph[case + 1] = graph[case + 1] + 1
                    is_new = False
            if is_new == True:
                new_graph_with_mult = [new_graph, 0, 0]
                new_graph_with_mult[case + 1] = 1
                graphs.append(new_graph_with_mult)
                found = found + 1
            counter = counter + 1
            # print (counter, found)
    print ("Number of isomorphism types of graphs:", found)
    return graphs

# Given a polyhedron P with a colouring and a state,
# it classifies all the ascending and descending links
# and calculates their reduced Betti numbers and/or fundamental groups

def get_links_topology(G, colouring, state, homology = True, fundamental_group = True):
    graphs = get_links_graphs(G, colouring, state)
    print ("----")
    Euler_sum = 0
    for g in graphs:
        cc = g[0].clique_complex()
        Euler_contribution = (1 - cc.euler_characteristic()) * g[1]
        if homology:
            print (cc.homology())
        if fundamental_group:
            if cc.is_connected():
                print (cc.fundamental_group())
            else:
                print ("Not connected")
        print ("Occurrences:", g[1], g[2])
        print ("Contribution to Euler characteristic:", Euler_contribution)
        Euler_sum = Euler_sum + Euler_contribution
        print ("----")
    print ("Euler characteristic of the manifold:", Euler_sum)

# Given a simplicial complex g,
# it searches for a free face and remove it and all the faces that contain it.
# Returns the collapse of the simplicial complex or the same simplicial complex (if no free face is found)

def find_and_collapse(g):
    for k in range(0,g.dimension()):
        for f in g.faces()[k]:
            counter=0
            for mf in g.maximal_faces():
                if f.dimension() < mf.dimension() and set(f).issubset(set(mf)):
                    counter=counter+1
            if counter==1:
                g.remove_face(tuple(f))
                return g
    return g

# Given a simplicial complex g, it applies the function find_and_collapse() until it does not find any face to collapse

def collapse_until_possible(g):
    check = True
    while check == True:
        R = deepcopy(g)
        g = find_and_collapse(g)
        if R.is_isomorphic(g):
            check = False
    return g

# Given a polyhedron P with a colouring and a state,
# It checks whether all ascending and descending links collapse to a subcomplex of dimension <= d (default is d=1)
# If connected == True, also checks that links are connected
# If simple == True, also checks that links are connected and with first Betti number <= 1

def links_collapse(G, colouring, state, d = 1, connected = True, simple = False):
    graphs = get_links_graphs(G, colouring, state)
    collapsed = []
    for g in graphs:
        cc = deepcopy(g[0].clique_complex())
        cc = collapse_until_possible(cc)
        collapsed.append(cc)
    for cc in collapsed:
        if cc.dimension()>d:
            return False
        if simple == True or connected == True:
            euler = cc.euler_characteristic()
            b0 = cc.betti(0) + 1 # Non capisco perche is_connected() non funziona e perche' bisogna aggiungere 1
            if b0 != 1:
                return False
            if simple == True and euler < 0:
                return False
    return True


# Claudio's addition

# The following code goes through all cusps of G and checks for all pairs of opposite facets if they lie
# in different connected components of the subgraph corresponding to their colorings (inside the graph
# induced by G). If they lie in the same component it prints "True" and "False" if the don't. It also
# prints the cusp and its number of colours

# The code first creates the list of cusps cuspidi and a separate list cuspcol containing their colorings
# Then it conducts the following steps for every cusp cuspidi[i]
# It first prints the cups and its number of colours
# Step 1: We collect the (d-1) pairs of colours on opposite facets of the cusp.
# If the pair of opposite facets only has one color, then we create a 1-tuple, else a 2-tuple of colors.
# Step 2: We go through the list of colours of pairs of opposite facets and create the subgraph
# of the graph defined by G which contains the vertices corresponding to all facets with these one/two colors
# Step 3: We create the flag simplicial complex associated with this subgraph.
# This step is probably not necessary, but I didn't know how to check using sage if two vertices of a graph
# are in the same connected component
# Step 4: For every pair of opposite facets of cuspidi[i] we identify the two vertices of the flag
# simplicial complex that correspond to this pair (f1 and f2).
# Then we apply the function connected_component to the pair of 0-simplices they provide to check if they
# lie in the same connected component. The code now prints "False" if they do not lie in the same connected component
# and "True" if they do lie in the same connected component

def get_cusp_morphisms(G,d,colouring):
    cuspidi = cusps(G, d, answer = [], level = 0, facets = [], n = 0)
# Create the list of colours of the cusps
    cuspcol=[]
    for cusp in cuspidi:
        colours = [colouring[cusp[i]] for i in range(2*(d-1))]
        cuspcol.append(colours)
    opposite_colours=[]
    dim= G.dimensions()[0]
    n=len(cuspidi)
# The code will now go through all the cusps of G and run through various steps.
# We will now fix a cusp cuspidi[i] on which Steps 1 - 4 will be conducted
    for i in range(n):
# Start of Step 1
        cusp=cuspidi[i]
        print("The cusp is: ", cusp)
        print("The cusps colours are:", cuspcol[i])
        colours = [cuspcol[i][j] for j in range(2*(d-1))]
        print("Its number of colours is: ", len(set(colours)))
        print("Its i-th pair of opposite facets lies in the same connected component of the subcomplex induced by their colors:")
        opposite_colours=[]
        for j in range(d-1):
            if cuspcol[i][j]!=cuspcol[i][j+(d-1)]:
                new_opp_col=[cuspcol[i][j],cuspcol[i][j+(d-1)]]
                opposite_colours.append(new_opp_col)
            else:
                new_opp_col=[cuspcol[i][j]]
                opposite_colours.append(new_opp_col)
# Start of Step 2
        subcomplexes=[]
        for colpair in opposite_colours:
            new_subcomplex=[]
            for l in colpair:
                for j in range(dim):
                    if colouring[j] == l:
                        new_subcomplex.append(j)
            subcomplexes.append(new_subcomplex)
        graphs=[]
        for complex in subcomplexes:
            G_reduced = G[complex,complex]
            new_graph = Graph(G_reduced, format='adjacency_matrix')
            graphs.append(new_graph)
# Start of Step 3
        ccs=[]
        for g in graphs:
            new_cc = g.clique_complex()
            ccs.append(new_cc)
# Start of Step 4
        for j in range(d-1):
            complex=subcomplexes[j]
            current_cc=ccs[j]
            f1=0
            f2=0
            for l in range(len(complex)):
                if complex[l]==cusp[j]: f1=l
                if complex[l]==cusp[j+(d-1)]: f2=l
            print(j+1,current_cc.connected_component((f1,))==current_cc.connected_component((f2,)))
    return

