### Creation of an IsogenyGraph class, so that we can package the relevant information
### i.e. prime, isogeny degree, vertex coloring information, etc. into one object and pass
### to other mehtods.

mpdb = ClassicalModularPolynomialDatabase()

# IsogenyGraph class inherits from sage.graphs.graph.Graph
class IsogenyGraph():
    def __init__(self, prime, isogeny_degree):
        self._graph = build_isogeny_graph_over_Fpbar(prime, isogeny_degree).to_undirected()

        # Relabel multiple edges to distinguish in path-length functions
        relabel_multiedges(self._graph)
        self._prime = prime
        self._isogeny_degree = isogeny_degree
    def __getattr__(self, attr):
        return getattr(self._graph, attr)
    def _repr_(self):
        return "Supersingular isogeny graph modulo %s with edges given by isogenies of degree %s"%(self._prime, self._isogeny_degree)
    def prime(self):
        return self._prime
    def isogeny_degree(self):
        return self._isogeny_degree

    # Method to return the maximal discriminants such that the elliptic curves with CM by the associated order
    # are supresingular modulo p. ubound and lbound are upper and lower bounds on discriminants returned, respectively.
    def embedded_fundamental_discriminants(self, ubound, lbound = 1):
        if lbound > ubound:
            raise ValueError("Lower bound is larger than upper bound")
        return [-d for d in [lbound..ubound] if is_fundamental_discriminant(-d) and (GF(self.prime())(-d).is_square() == False or d % self.prime() == 0)]
 
    # Method to return a random maximal discriminant such that the elliptic curves with CM by
    # the associated maximal order are supersingular mod p. ubound and lbound are upper and lower bounds on discriminant returned. 
    def random_fundamental_discriminant(self, ubound, lbound = 1):
        if lbound > ubound:
            raise ValueError("Lower bound is larger than upper bound")
        embedded_fundamental_discriminants = self.embedded_fundamental_discriminants(ubound, lbound)
        return embedded_fundamental_discriminants[randrange(0, len(embedded_fundamental_discriminants))] 

    # Method to return adjacency matrix of G, but with loops counted as two
    def adjusted_AM(self): 
        AM = copy(self.adjacency_matrix())
        for i in range(AM.dimensions()[1]):
            if AM[i][i] > 0:
                AM[i,i] = 2*AM[i,i]
        return AM




### Code from "Adventures in Supersingularland" https://arxiv.org/abs/1909.07779
### Authors: Sarah Arpin, Catalina Camacho-Navarro, Kristin Lauter, Joelle Lim, Kristina Nelson, Travis Scholl, Jana Sotáková.

### ELI: Minor modification to use Kohel database instead of hardcoded modular polynomials


def build_isogeny_graph_over_Fpbar(p, l, steps=oo):
    """
    Given a prime p, this function returns
    the l-isogeny graph of supersingular elliptic
    curves over bar(F_p).

    If given, this only gives the graph up to "steps" number
    of edges from an origin curve.

    The algorithm works by first finding any supersingular
    j-invariant via an algorithm of Broker, then walking the isogeny
    graph by a BFS.

    References:
     - Constructing Supersingular Elliptic Curves - Reinier Broker
    """
    # STEP 1: Find a single super-singular curve
    """
    Given a prime p >= 5, this step finds a supersingular elliptic
    curve over F_{p^2}.
    """
    q = next(q for q in Primes() if q%4 == 3 and kronecker_symbol(-q,p) == -1)
    K = QuadraticField(-q)
    H = K.hilbert_class_polynomial()
    j0 = H.change_ring(GF(p^2)).any_root()
    # STEP 2: Walk along the isogeny graph
    """
    Two elliptic curves E1,E2 are l-isogenous (over bar(F_p)) if and
    only if x=j(E1),y=j(E2) are a root of the l-modular polynomial
    Phi_l. Tables of Phi_l can be found online
    (see https://math.mit.edu/~drew/ClassicalModPolys.html) and in Sage
    via `ClassicalModularPolynomialDatabase()`.
    """
    phi = mpdb[l]
    def get_neighbors(j):
        """
        This function returns a list of all roots of Phi_l(j,X),
        repeated with appropiate multiplicity.
        """
        R.<x> = GF(p^2)[]
        return flatten([[j2]*k for j2,k in phi(j,x).roots()])
    G = DiGraph(multiedges=True,loops=True)
    visited = set()
    not_visited = set([j0])
    count = 0
    while not_visited:
        j1 = not_visited.pop()
        visited.add(j1)
        for j2 in get_neighbors(j1):
            if j1 < j2 or j1 == j2:
                G.add_edge([j1,j2])
            if j2 not in visited and j2 not in not_visited:
                not_visited.add(j2)
        count += 1
        if count == steps:
            break
    return G


### Function to relabel all multiple edges in a graph so that each edge has a unique (numeric) label
## This is useful when searching for paths in the SSl-I graph and disallowing backtracking, and is used in the IsogenyGraph constructor
## Inputs: G - graph object
## Outputs: None - the function directly modifies the input

def relabel_multiedges(G):

    # Iterator over all multiple edges in G
    i = 0

    # store new, labelled edges in list to append at the end
    new_edges = []

    # list of original, unlabelled multiedges, sorted so that there are no (u,v,None) and (v,u,None) 
    old_edges = G.multiple_edges()

    while i < len(old_edges):

        # iterator for sublist of edges between the same two vertices
        j = 0

        #reference vertices
        ref_vert_0 = old_edges[i][0]; ref_vert_1 = old_edges[i][1]

        # NOTE: This isn't a good way to do this, I'm expecting the list from G.multiple_edges() is nicely ordered
        #       A better solution would be to find the list of edges that have (u,v) or (v,u), iterate over all, deleting as I go
        #       and then go to the next reference vertex. I'll update this someday in the future.
        while (i < len(old_edges)) and (old_edges[i][0], old_edges[i][1]) == (ref_vert_0, ref_vert_1):
            G.delete_edge((old_edges[i][0], old_edges[i][1], None))
            new_edges.append((old_edges[i][0], old_edges[i][1], j))
            j += 1
            i += 1

    # Add in new (labelled) edges
    for e in new_edges:
        G.add_edge(e)


### ELI: Function to return vertices with endomorphisms by a given maximal order

## Inputs: G - isogeny graph, d - fundamental discriminant
def get_CM_vertices(G, d):
    p = G.prime()

    # Create quadratic field with discriminant d
    K = QuadraticField(d)

    # Create HCP for quadratic field above
    h = K.hilbert_class_polynomial()

    #Find roots of HCP in F_p^2
    hmodroots = h.change_ring(GF(p^2)).roots()

    #return vertices found above
    return [hmodroots[i][0] for i in (0..len(hmodroots)-1)]

### ELI: Function to color the vertices of an isogeny graph that have endomorphism rings by
###      the maximal order O_d1, O_d2 for fundamental discriminants d1, d2

## Inputs: G - isogeny graph, d1, d2 - fundamental discriminants, vertex_labels = whether to add labels or not in plot object
##         default input for d2 is None, if this is left, then only vertices from d1 are colored.
def color_isogeny_graph(G, d1, d2 = 0, vertex_labels = True):
    vertexset1 = get_CM_vertices(G, d1)
    if d2: vertexset2 = get_CM_vertices(G, d2)

    # Create dictionary for coloring
    if d2:
        d = {'#FF0000' : vertexset1, '#0000FF' : vertexset2}
    else:
        d = {'#FF0000' : vertexset1}

    # Return plot object with colors added
    return G.graphplot(vertex_colors = d, vertex_labels = vertex_labels)