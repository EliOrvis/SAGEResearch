##### This file contains methods for investigating modular polynomials 

mpdb = ClassicalModularPolynomialDatabase()

### Function to find all lifts and associated primes of a j-invariant from a
### SSI graph to a j-invariant over C defined in the HCF associated to a maximal 
### order embedded in the endomorphism ring.

## Inputs: G - SSI graph, v - vertex, d - maximal order where the vertex has CM 
def find_CM_lifts(G, v, d):
	# Initiate everything
	Q.<sqrtd> = QuadraticField(d)
	HCP = Q.hilbert_class_polynomial()
	H.<j0> = Q.extension(HCP)

	HCF_primes = H.primes_above(G.prime())
	j_invars = j0.galois_conjugates(H)

	lifts = []

	# For each prime and j-invariant, we save (j, hom) if j reduces to v under hom.
	# An important caveat is that there is a choice of map from the Residue field to
	# F_p^2. We therefore save a list of lists, where each sublist contains all (j,hom) for a fixed
	# choice of map from the residue field at p to F_p^2. hom is the map from residue field to F_p^2.
	# Note that the prime of reduction can be recovered from hom via hom.domain().ideal().
	for p in HCF_primes:
		ResField = p.residue_field()
		Homs = Hom(ResField, GF(G.prime()^2))
		for j in j_invars:
			for hom in Homs:
				if hom(j) == v:
					lifts.append((j, hom))
	return lifts

### Find the number of primes dividing phi_ell^n(j_1, j_2), where j_2 ranges over all j-invariants in a given valley
## Inputs: G - isogeny graph object; d1, d2 - embedded fundamental discriminants; length - path-length to consider (i.e., ell^n = ell^length)
##			fixed_prime - whether to fix a (random) prime of reduction instead of fixing a j-invariant
## Outputs: list of triples (P, j_1, j_2) where P divides phi_ell^n(j_1,j_2)

# NOTE: Currently only implemented for when d1 == d2
def primes_dividing_modpoly_pseudonorm(G, d1, d2, length, fixed_prime = False):
	# Initiate everything
	Q.<sqrtd> = QuadraticField(d)
	HCP = Q.hilbert_class_polynomial()
	H.<j0> = Q.extension(HCP)

	HCF_primes = H.primes_above(G.prime())
	if fixed_prime == True:
		P = HCF_primes[randrange(0, len(HCF_primes))]

	j_invars = j0.galois_conjugates(H)

	phi = mpdb[G.isogeny_degree()^length]

	primes_dividing_pairs = []

	if fixed_prime == True:
		for j1 in j_invars:
			for j2 in j_invars:
				mod_poly_ideal = H.ideal(phi(j1, j2))
				if mod_poly_ideal.is_coprime(P) == False:
					primes_dividing_pairs.append((P, j1, j2))
	else:
		for j in j_invars:
			mod_poly_ideal = H.ideal(phi(j0,j))
			for prime in HCF_primes:
				if mod_poly_ideal.is_coprime(prime) == False:
					primes_dividing_pairs.append((prime, j0, j))

	return primes_dividing_pairs

