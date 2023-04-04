### Creation of an IsogenyGraph class, so that we can package the relevant information
### i.e. prime, isogeny degree, vertex coloring information, etc. into one object and pass
### to other mehtods.

# IsogenyGraph class inherits from sage.graphs.graph.Graph
class IsogenyGraph():
    def __init__(self, prime, isogeny_degree):
        self._graph = build_isogeny_graph_over_Fpbar(prime, isogeny_degree).to_undirected()
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
    def plot(self):
        return self.graph().plot()

    # Method to return the maximal discriminants such that the elliptic curves with CM by the associated order
    # are supresingular modulo p. ubound and lbound are upper and lower bounds on discriminants returned, respectively.
    def embedded_fundamental_discriminants(self, lbound, ubound):
        if lbound > ubound:
            raise ValueError("Lower bound is larger than upper bound")
        return [-d for d in [lbound..ubound] if is_fundamental_discriminant(-d) and GF(G.prime())(-d).is_square() == False]

    # Method to return a random maximal discriminant such that the elliptic curves with CM by
    # the associated maximal order are supersingular mod p. ubound and lbound are upper and lower bounds on discriminant returned. 
    def random_fundamental_discriminant(self, lbound, ubound):
        if lbound > ubound:
            raise ValueError("Lower bound is larger than upper bound")
        embedded_fundamental_discriminants = G.embedded_fundamental_discriminants(lbound, ubound)
        return embedded_fundamental_discriminants[randrange(0, len(embedded_fundamental_discriminants))] 




### Code from "Adventures in Supersingularland" https://arxiv.org/abs/1909.07779
### Authors: Sarah Arpin, Catalina Camacho-Navarro, Kristin Lauter, Joelle Lim, Kristina Nelson, Travis Scholl, Jana Sotáková.

### ELI: Minor modification to return an IsogenyGraph object as created above, rather than a sage.graphs.graph.Graph object -
###      I force the undirected graph, not sure why the original authors chose to create the graph directed.

### harcoded modular polynomials for ell = 2,3

R.<X,Y> = ZZ[]
ModularPoly = {}
ModularPoly[2] = X^3 - X^2*Y^2 + 1488*X^2*Y - 162000*X^2 + 1488*X*Y^2 + 40773375*X*Y + 8748000000*X + Y^3 - 162000*Y^2 + 8748000000*Y -157464000000000
ModularPoly[3] = X^4 - X^3*Y^3 + 2232*X^3*Y^2 - 1069956*X^3*Y + 36864000*X^3 + 2232*X^2*Y^3 + 2587918086*X^2*Y^2 + 8900222976000*X^2*Y +452984832000000*X^2 - 1069956*X*Y^3 + 8900222976000*X*Y^2 - 770845966336000000*X*Y + 1855425871872000000000*X + Y^4 + 36864000*Y^3 + 452984832000000*Y^2 + 1855425871872000000000*Y
ModularPoly[5] = X^6 - X^5*Y^5 + 3720*X^5*Y^4 - 4550940*X^5*Y^3 + 2028551200*X^5*Y^2 - 246683410950*X^5*Y + 1963211489280*X^5 +3720*X^4*Y^5 + 1665999364600*X^4*Y^4 + 107878928185336800*X^4*Y^3 + 383083609779811215375*X^4*Y^2 + 128541798906828816384000*X^4*Y + 1284733132841424456253440*X^4 - 4550940*X^3*Y^5 + 107878928185336800*X^3*Y^4 -     441206965512914835246100*X^3*Y^3 + 26898488858380731577417728000*X^3*Y^2 - 192457934618928299655108231168000*X^3*Y +     280244777828439527804321565297868800*X^3 + 2028551200*X^2*Y^5 + 383083609779811215375*X^2*Y^4 +     26898488858380731577417728000*X^2*Y^3 + 5110941777552418083110765199360000*X^2*Y^2 +     36554736583949629295706472332656640000*X^2*Y + 6692500042627997708487149415015068467200*X^2 - 246683410950*X*Y^5 +     128541798906828816384000*X*Y^4 - 192457934618928299655108231168000*X*Y^3 +     36554736583949629295706472332656640000*X*Y^2 - 264073457076620596259715790247978782949376 *X*Y +     53274330803424425450420160273356509151232000*X + Y^6 + 1963211489280*Y^5 + 1284733132841424456253440*Y^4 +     280244777828439527804321565297868800*Y^3 + 6692500042627997708487149415015068467200*Y^2 +     53274330803424425450420160273356509151232000*Y + 141359947154721358697753474691071362751004672000


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
    phi = ModularPoly[l]
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
def color_isogeny_graph(G, d1, d2, vertex_labels = True):
    vertexset1 = get_CM_vertices(G, d1); vertexset2 = get_CM_vertices(G, d2)

    # Create dictionary for coloring
    d = {'#FF0000' : vertexset1, '#0000FF' : vertexset2}

    # Return plot object with colors added
    return G.graph().graphplot(vertex_colors = d, vertex_labels = vertex_labels)