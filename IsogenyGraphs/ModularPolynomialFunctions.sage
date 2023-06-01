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
## Outputs: list of triples (P, j_1, j_2) where P divides phi_ell^n(j_1,j_2)

def primes_dividing_modpoly_pseudonorm(G, d1, d2, length):
	# Initiate everything
	Q1.<sqrtd1> = QuadraticField(d1); Q2.<sqrtd2> = QuadraticField(d2)
	HCP1 = Q1.hilbert_class_polynomial(); HCP2 = Q2.hilbert_class_polynomial()
	H1.<j01> = Q1.extension(HCP1)

	H1H2.<j02> = H1.extension(HCP2)

	HCF_primes = H1H2.primes_above(G.prime())

	j_invars_1 = j01.galois_conjugates(H1H2)
	j_invars_2 = j02.galois_conjugates(H1H2)

	phi = mpdb[G.isogeny_degree()^length]

	primes_dividing_pairs = []

	for j in j_invars_2:
		mod_poly = phi(j01,j)
		if mod_poly == 0:
			for prime in HCF_primes:
				primes_dividing_pairs.append((prime, j01, j))
		else:
			mod_poly_ideal = H1H2.fractional_ideal(mod_poly)
			for prime in HCF_primes:
				if mod_poly_ideal.is_coprime(prime) == False:
					primes_dividing_pairs.append((prime, j01, j))

	return primes_dividing_pairs

