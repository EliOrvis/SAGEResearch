#### Create affine PPT graph for a given prime p

# Get solutions to a^2 + b^2 - c^2 mod p as vectors
def get_solutions(p):
    R.<a,b,c> = PolynomialRing(GF(p))
    Id = Ideal(a^2 + b^2 - c^2, a^p - a, b^p - b, c^p - c)
    L = Id.variety()
    L.remove({a:0,b:0,c:0})
    sols = list()
    for i in range(len(L)):
        sols.append(vector([L[i][a], L[i][b], L[i][c]]))
    return sols

 # Helper function
def immutable(A):             
    B=copy(A)
    B.set_immutable()
    return B

# Get edges given vertices - swap indicates whether we connect (a,b,c) with (b,a,c) or not, negate indicates whether we connect (a,b,c) with (-a,-b,-c)
def get_edges(p, verts, swap, negate):
    M1 = Matrix(GF(p), [[-1,2,2],[-2,1,2], [-2,2,3]])
    M2 = Matrix(GF(p), [[1,2,2],[2,1,2],[2,2,3]])
    M3 = Matrix(GF(p), [[1,-2,2],[2,-1,2],[2,-2,3]])
    M4 = Matrix(GF(p), [[0,1,0],[1,0,0],[0,0,1]])
    M5 = Matrix(GF(p), [[-1,0,0],[0,-1,0],[0,0,-1]])
    edges = list()
    for i in range(len(verts)):
        edges.append((immutable(verts[i]), immutable(M1*verts[i])))
        edges.append((immutable(verts[i]), immutable(M2*verts[i])))
        edges.append((immutable(verts[i]), immutable(M3*verts[i])))
        if swap == True:
            edges.append((immutable(verts[i]), immutable(M4*verts[i])))
        if negate == True:
            edges.append((immutable(verts[i]), immutable(M5*verts[i])))
    return edges
    
# Create graph given only p, whether to include swaps, whether to include negations, whether to allow multiple edges, whether to include loop vertices
def create_PPT_graph(p, swap, negate, multiedge):
    verts = get_solutions(p)
    edges = get_edges(p, verts, swap, negate)
    # This has to be done because Sage can't handle mutable vectors as vertices in a graph
    for i in range(len(verts)):
        verts[i].set_immutable()
    G = Graph([verts, edges], format='vertices_and_edges', multiedges = multiedge, loops = True)
    return G

#### Create projective version of PPT graph

#get projective solutions
def get_proj_sols(p):
    P.<a,b,c> = ProjectiveSpace(2, GF(p))
    V = P.subscheme([a^2 + b^2 - c^2])
    return V.rational_points()

# get projective edges
def get_proj_edges(p, verts):
    P.<a,b,c> = ProjectiveSpace(2, GF(p))
    M1 = P.hom([-a + 2*b + 2*c, -2*a + b + 2*c, -2*a + 2*b + 3*c], P)
    M2 = P.hom([a + 2*b + 2*c, 2*a + b + 2*c, 2*a + 2*b + 3*c], P)
    M3 = P.hom([a - 2*b + 2*c, 2*a - b + 2*c, 2*a - 2*b + 3*c], P)
    edges = list()
    for i in range(len(verts)):
        edges.append((verts[i], M1(verts[i])))
        edges.append((verts[i], M2(verts[i])))
        edges.append((verts[i], M3(verts[i])))
    return edges
    
# create projective PPT graph
def create_proj_PPT_graph(p, multiedge):
    verts = get_proj_sols(p)
    edges = get_proj_edges(p, verts)
    return Graph([verts, edges], multiedges = multiedge, loops = True)


 #### Helper function to get correct eigenvalues of adjacency matrices

 # Adjust adjacency matrix to count loops as degree 2
def get_adjusted_AM(G):
    AM = copy(G.adjacency_matrix())
    for i in range(AM.dimensions()[1]):
        if AM[i][i] > 0:
            AM[i,i] = 2*AM[i,i]
    return AM