def kerpol_pushforward(h, phi):
    """Returns the kernel polynomial of the pushforward along phi of the kernel cut out by h.

    h is the kernel polynomial of a n-isogeny with kernel K, phi is an ell-isogeny with gcd(n, ell)=1. This function
    returns the kernel polynomial of the isogeny of E/kerphi with kernel phi(K).
    """
    phix = phi.x_rational_map()
    u = phix.numerator()
    A = u.parent()
    v = phix.denominator()
    F = h.base_ring()
    R.<T> = PolynomialRing(F)
    S.<X> = PolynomialRing(R)
    U = u(X)
    V = v(X)
    H = h(X)
    R = (H.resultant(T*V - U)) // (H.resultant(V))
    return A(R)

def make_orbit(E, h):
    """returns orbit of kernel polynomial h under automorphism group of E"""
    Fp2 = E.base_ring()
    if E.j_invariant() not in [Fp2(0), Fp2(1728)]:
        return [h]
    if E.j_invariant() == 0:
        Fp2 = E.base_ring()
        Fp2x.<x> = PolynomialRing(Fp2)
        d = h.degree()
        omega = (x^2+x+1).roots()[0]
        # remove duplicates of h if h is stable under Aut(E)
        orbit = [(omega^{-i * d}) * h(omega^i * x) for i in range(3)]
        orbit = list(dict.fromkeys(orbit))
        orbit.sort()
        return orbit
    if E.j_invariant() == 1728:
        Fp2 = E.base_ring()
        Fp2x.<x> = PolynomialRing(Fp2)
        d = h.degree()
        orbit = [((-1)^(i * d)) * h((-1)^i * x) for i in range(2)]
        # remove duplicates of h if h is stable under Aut(E)
        orbit = list(dict.fromkeys(orbit))
        orbit.sort()
        return orbit
def lift_to_borel(E, n):
    '''Returns all vertices in G(p,ell,n) above E, a vertex in G(p,ell).

    Vertices (E,C) are represented as (E, h) where h is a kernel.

    TODO: implement for non-prime / non-odd n. Probably does not work if p = 2 mod 3, n = 3 or p = 3 mod 4, n = 2.
    '''
    Fp2 = E.base_ring()
    p = Fp2.characteristic()
    assert p > 0, "characteristic must be positive"
    assert Fp2 == GF(p^2), "must work over finite field of p^2 elements"
    assert E.trace_of_frobenius() in [-2*p, 2*p], "must work in trace = +/- 2p isogeny class"
    isogenies = E.isogenies_prime_degree(n)
    kerpols = [phi.kernel_polynomial() for phi in isogenies]
    if E.j_invariant() not in [Fp2(0), Fp2(1728)]:
        return [(E, h) for h in kerpols]
    else:
        isomorphism_reps = []
        done = []
        for h in kerpols:
            if h in done:
                continue
            orbit = make_orbit(E, h)
            done = done + orbit
            isomorphism_reps.append((E, orbit[0]))
        return isomorphism_reps
        
def isogeny_graph_borel(p, ell, n):
    '''
    Returns the directed graph of ell-isogenies of supersingular elliptic curves with level n structure in characteristic p.

    Parameters:
        p: prime > 3, the characteristic
        ell: prime not equal to p, the degree of the isogenies
        n: the level, coprime to p*ell. Currently required to be prime.
    Returns:
        G: a directed graph
    '''
    Fp2 = GF(p^2)
    j = supersingular_j(Fp2)
    E = EllipticCurve(j=j)
    curves = {j:E}
    lifts = lift_to_borel(E, n)
    stack = lifts
    vertices = []
    graph_dict = {}
    while stack:
        (E, h) = stack.pop()
        isogenies = E.isogenies_prime_degree(ell)
        neighbors = []
        for phi in isogenies:
            Ephi = phi.codomain()
            j = Ephi.j_invariant()
            if j in curves:
                u = Ephi.isomorphism_to(curves[j])
                phi = u * phi
                Ephi = curves[j]
            else:
                curves[j] = Ephi
            hphi = kerpol_pushforward(h, phi)
            orbit_phi = make_orbit(Ephi, hphi)
            hphi = orbit_phi[0]
            if (Ephi, hphi) not in vertices:
                stack.append((Ephi, hphi))
            neighbors.append((Ephi.a4(), Ephi.a6(), hphi))
        vertices.append((E, h))
        graph_dict[(E.a4(), E.a6(), h)] = neighbors
    return DiGraph(graph_dict)