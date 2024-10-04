### Creation of an IsogenyGraph class, so that we can package the relevant information
### i.e. prime, isogeny degree, vertex coloring information, etc. into one object and pass
### to other mehtods.

mpdb = ClassicalModularPolynomialDatabase()

# This loads James's isogeny graph code, which is used to construct isogeny graphs because it is much faster than the sage version.
# Note that there is a failsafe in the isogeny graph constructor, in case this fails for some reason (like the file moving or something).
try:
  load("~/Isogeny/ssl_pari.sage")
  load_success = True
except:
  print("Warning: Unable to load James' Pari code.")  
  load_success = False

# Main class for isogeny graphs
# IsogenyGraph class inherits from sage.graphs.graph.Graph
class IsogenyGraph():
    def __init__(self, prime, isogeny_degree, undirected = False):
        if undirected == True:
            self._graph = build_isogeny_graph_over_Fpbar(prime, isogeny_degree, True)
        else:
            self._graph = build_isogeny_graph_over_Fpbar(prime, isogeny_degree, False)

        # Relabel multiple edges to distinguish in path-length functions
        #frelabel_multiedges(self._graph)
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
    # Method to return the number of spine vertices. 
    # This can be faster and easier than finding all vertices on the spine directly for large graphs
    def n_spine_verts(self):
        return get_n_spine_verts_from_p(self.prime())

    # Method to return the subgraph induced by the F_p-vertices, i.e. the spine.
    def spine_subgraph(self):
        return self.subgraph(vertices = [v for v in self.vertices() if v^p == v])

    # Method to return the maximal discriminants such that the elliptic curves with CM by the associated order
    # are supresingular modulo p. ubound and lbound are upper and lower bounds on discriminants returned, respectively.
    def embedded_fundamental_discriminants(self, ubound, lbound = 1):
        if lbound > ubound:
            raise ValueError("Lower bound is larger than upper bound")
        return [-d for d in [lbound..ubound] if (-d).is_fundamental_discriminant() and (GF(self.prime())(-d).is_square() == False or d % self.prime() == 0)]
 
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

### ELI: Minor modification to use Kohel database instead of hardcoded modular polynomials, added undirected flag to give the option of an appropriate undirected
####     graph (somewhat messed up possibly at 0 and 1728)
####     Also, major modification to use James' code if possible, which saves substantial time on large primes.


def build_isogeny_graph_over_Fpbar(p, l, undirected, steps=oo):
    # First, create the graph with James' code if possible, otherwise, we use the slower SAGE implementation
    if load_success == True and undirected == False:
        return ssl_graph(p, l)
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
    if undirected == True:
        G = Graph(multiedges=True,loops=True)
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
        return G.to_undirected()
    else:
        G = DiGraph(multiedges=True,loops=True)
        visited = set()
        not_visited = set([j0])
        count = 0
        while not_visited:
            j1 = not_visited.pop()
            visited.add(j1)
            for j2 in get_neighbors(j1):
                G.add_edge([j1,j2])
                if j2 not in visited and j2 not in not_visited:
                    not_visited.add(j2)
            count += 1
            if count == steps:
                break
        return G

### This function builds the undirected isogeny graph for a given prime and isogeny degree
##  Inputs: p - a prime number; ell - a prime number other than <p>
##  Outputs: G - an undirected graph object with edges labelled by a representative for the isogeny class of that edge
#   NOTE: This is MUCH slower than the directed version, since in this version we compute all the isogenies. Use with caution.
def build_undirected_isogeny_graph_over_Fpbar(p, ell):
    # Note: This is attrocious, but it works. We compute the undirected graph, then get the isogenies, then figure out the equivalences
    #       I will fix this all someday, but right now it is left as is so that I have time to write a dissertation...
    und_graph = IsogenyGraph(p, ell)
    isogenies = get_isogenies(und_graph)
    print(len(isogenies))

    # The isogenies are already up to post-composition by an automorphism, so we only need to identify by pre-composition, then by duals
    # These functions are used for the identifications
    def id_by_precomp(phi1, phi2):
        auts = phi1.domain().automorphisms()

        # If phi1 is equal to phi2 after precomposing by an automorphism, return true, otherwise, return false
        if any([phi1.pre_compose(aut) == phi2 for aut in auts]):
            return True
        else:
            return False

    # THIS IS WHERE THE PROBLEM IS. IDENTIFYING UP TO DUALS IS NOT WORKING CORRECTLY SOMEHOW
    def id_by_dual(iso_class1, iso_class2):
        if any([phi.dual() in iso_class2 for phi in iso_class1]):
            return True
        else:
            return False

    # 
    edge_preaut_eq_classes = []
    for j1 in und_graph.vertices():
        for j2 in und_graph.vertices():
            j1j2_isos = [iso for iso in isogenies if iso.domain().j_invariant() == j1 and iso.codomain().j_invariant() == j2]
            # Put these into equivalence classes
            edge_preaut_eq_classes = edge_preaut_eq_classes + equal_classes(j1j2_isos, id_by_precomp)


    print(len(edge_preaut_eq_classes))
    # Now identify all the classes by duals - i.e. if two classes contain isogenies that are dual to one another, they are now considered the same
    dual_idd_edges = equal_classes(edge_preaut_eq_classes, id_by_dual)
    print(len(dual_idd_edges))

    # Now make the labelled edge set for the undirected graph
    edges = {}

    for j1 in und_graph.vertices():
        edges[j1] = {j2 : [edge_class[0][0] for edge_class in dual_idd_edges if edge_class[0][0].domain().j_invariant() == j1 and edge_class[0][0].codomain().j_invariant() == j2] for j2 in und_graph.vertices()}

    Gu = Graph(edges)
    return Gu

### Function to return the number of vertices in the SS isogeny graph over F_p-bar. Does the obvious thing from Silverman.
##  Inputs: p - prime number
##  Outputs: number of vertices in SS isogeny graph over F_p bar.
def n_ss_j_invars(p):

    base = floor(p/12)
    if p % 12 == 0:
        return base
    elif p % 12 in [5, 7]:
        return base + 1
    else:
        return base + 2

### Function to return the number of vertices  w/ j-invar in F_p in the SS isogeny graph over F_p-bar.
### Uses remark after Definition 2.4 in the Adventures in SSLand paper.
##  Inputs: p - prime number
##  Outputs: number of vertices in SS isogeny graph over F_p bar that are defined over F_p.
def n_Fp_j_invars(p):

    # Two cases for the discriminant from paper cited above.
    if p % 4 == 1:
        d = -4*p
        K = QuadraticField(d)
        return (1/2)*K.class_number()
    else:
        d = -p
        K = QuadraticField(d)
        if p % 8 == 7:
            return K.class_number()
        else:
            return 2*K.class_number()




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


### ELI: Function to return vertices with optimal embeddings of endomorphisms by a given order

## Inputs: G - isogeny graph, d - fundamental discriminant
def get_CM_vertices(G, d):
    p = G.prime()

    # Create HCP for quadratic discriminant above
    h = hilbert_class_polynomial(d)

    #Find roots of HCP in F_p^2
    hmodroots = h.change_ring(GF(p^2)).roots()

    #return vertices found above
    return [hmodroots[i][0] for i in (0..len(hmodroots)-1)]

###  Function to return the number of spine vertices, without computing the graph. See Adventures paper for justification

## Inputs: p - odd prime number;
## Outputs: integer number of spine vertices for G(p,ell). This number is independent of ell
def get_n_spine_verts_from_p(p):
    # minor validation step:
    if p.is_prime() == False or p == 2:
         raise TypeError("Input must be an odd prime.")
    
    # Compute quadratic field
    K = QuadraticField(-p)

    # Depending on p mod 8, we return either the class number, twice the class number, or half the class number
    if p % 4 == 1:
        return (1/2)*K.class_number()
    else:
        if p % 8 == 7:
            return K.class_number()
        if p % 8 == 3:
            return 2*K.class_number()

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


## ELI: Function to add color on the spine of the isogeny graph

## Inputs: G - isogeny graph; vertex_labels - whether to label vertices in output plot
## Outputs: graphplot object
def color_spine(G, vertex_labels = True):
    Fp_verts = [vert for vert in G.vertices() if vert^p == vert]

    return G.graphplot(vertex_colors = {'#868686' : Fp_verts}, vertex_labels = vertex_labels)

### Function to return the loops in the isogeny graph, as acutal isogenies
##  Inputs: G - Isogeny graph
##  Outputs: loops - list of isogenies that represent loops in the graph
def get_isogeny_loops(G):
    # For each vertex, find all isogenies, then return only the ones where the j-invariant of the image is the same as the domain
    # Note that this requires code from my EllipticCurveFunctions file

    verts = G.vertices()
    loops = []
    p = G.prime()
    for vert in verts:
        E = EllipticCurve_from_j(GF(p^2)(vert))
        isos = isogenies_of_degree_n(E, G.isogeny_degree())
        for iso in isos:
            if iso.codomain().j_invariant() == iso.domain().j_invariant():
                loops.append(iso)

    return loops

### Function to return a complete set of representatives for the isogenies in the SS ell-isogeny graph, as actual isogenies
##  Inputs: G - Isogeny graph
##  Outputs: isogenies - list of isogenies that represent all edges in the graph
##  NOTE: This requires the function <isogenies_of_degree_n> from my EllipticCurveFunctions.sage file
def get_isogenies(G):
    # For each vertex, find all isogenies
    # Note that this requires code from my EllipticCurveFunctions file

    verts = G.vertices()
    isogenies = []
    p = G.prime()
    # Get models for all curves in the isogeny graph
    models = [EllipticCurve_from_j(GF(p^2)(vert)).base_extend(GF(p^(12))) for vert in verts] 
    for E in models:
        isos = isogenies_of_degree_n(E, G.isogeny_degree(), model_set = models)
        isogenies.append(isos)

    return flatten(isogenies)    

### Function to return only the self-dual loops. Note that this currently has a small chance of crashing out because
### of a sage bug. We're working on it.
### NOTE: This requires more code from my EllipticCurveFunctions file.
##  Inputs: G - Isogeny graph
##  Otuputs: sd_loops - list of loops
def get_self_dual_isogeny_loops(G):

    # First get all loops
    loops = get_isogeny_loops(G)

    # Return the self-dual ones
    return [loop for loop in loops if is_self_dual_ss(loop)]

### Function to return the matrix of the action of duals on the isogenies in G
##  Inputs: G - Isogeny graph
##  Outputs: M - matrix of the actional of taking duals on the isogeny graph
def get_dual_matrix(G):
    # Get all isogenies in G
    isos = get_isogenies(G)

    # Make matrix
    M = Matrix(ZZ, len(isos), lambda i,j : 1 if any([isos[i].dual() == isos[j].post_compose(aut) for aut in isos[j].codomain().automorphisms()]) else 0)

    return M
