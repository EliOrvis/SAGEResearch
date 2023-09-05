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
##           random_j1 - whether to randomly choose a j1 from H1 to be the fixed input to the modular polynomial
## Outputs: list of triples (P, j_1, j_2) where P divides phi_ell^n(j_1,j_2)
def primes_dividing_modpoly_pseudonorm(G, d1, d2, length, random_j1 = False):

	if d1 != d2:
		Q1.<sqrtd1> = QuadraticField(d1); Q2.<sqrtd2> = QuadraticField(d2)
		HCP1 = Q1.hilbert_class_polynomial()
		HCP2 = Q2.hilbert_class_polynomial()
		H1.<j1> = Q1.extension(HCP1)

		H1opt, from_H1opt, to_H1opt = H1.absolute_field(names = 'H1optgen').optimized_representation()

		# Factor HCP2 over H1opt so that we adjoin the root of an irreducible polynomial
		HCP2factor = HCP2.change_ring(H1opt).factor()[0][0]

		H1j2.<j2> = H1opt.extension(HCP2factor) 

		H1j2opt, from_H1j2opt, to_H1j2opt = H1j2.absolute_field(names = 'H1j2optgen').optimized_representation()

		# Factor x^2 - d2, in case sqrtd2 is already in H1j2opt
		R.<x> = ZZ[]
		f = x^2 - d2
		sqrtd2factor = f.change_ring(H1j2opt).factor()[0][0]

		H1H2.<sqrtd2> = H1j2opt.extension(sqrtd2factor)

		H1H2opt, from_H1H2opt, to_H1H2opt = H1H2.absolute_field(names = 'H1H2optgen').optimized_representation()

		# Find a j-invariant for each HCF in the new (optimized) composite field
		j01 = HCP1.change_ring(H1H2opt).any_root()
		j02 = HCP2.change_ring(H1H2opt).any_root()

	if d1 == d2:
		Q.<sqrtd> = QuadraticField(d1)
		HCP = Q.hilbert_class_polynomial()

		# When d1 == d2, we have only one Hilbert Class Filed, but we call it H1H2 to agree with the other name downstream
		H1H2.<j02> = Q.extension(HCP)

		H1H2opt, from_H1H2opt, to_H1H2opt = H1H2.absolute_field(names = 'H1H2optgen').optimized_representation()

		j01 = HCP.change_ring(H1H2opt).any_root()

	HCF_primes = H1H2opt.primes_above(G.prime())

	j_invars_2 = j02.galois_conjugates(H1H2opt)

	if random_j1 == True:
		j_invars_1 = j01.galois_conjugates(H1H2opt)
		j01 = j_invars_1[randrange(0, len(j_invars_1))] 

	phi = mpdb[G.isogeny_degree()^length]

	primes_dividing_pairs = []

	for j in j_invars_2:
		mod_poly = phi(j01,j)
		if mod_poly == 0:
			for prime in HCF_primes:
				primes_dividing_pairs.append((prime, j01, j))
		else:
			mod_poly_ideal = H1H2opt.fractional_ideal(mod_poly)
			for prime in HCF_primes:
				if mod_poly_ideal.is_coprime(prime) == False:
					primes_dividing_pairs.append((prime, j01, j))

	return primes_dividing_pairs


### Compute the correction factor between the number of pairs and the count of primes dividing the modular polynomial pseudonorm
## Inputs: d1, d2 - embedded fundamental discriminants;
## Outputs: integer correction factor

# Note: d1 and d2 must be input in the same order that they were input to the pair coutning function and the prime counting function for (silly) technical reasons.
def prime_pair_correction_factor(d1, d2):
	# when d1 == d2, no correction is needed
	if d1 == d2:
		return 1

	# Otherwise, the correction factor is given by 2^{r-1}/h(d_2), where 2^r is the degree of the intersection H_1 \cap H_2
	else:
		K2 = QuadraticField(d2)
		hd2 = K2.class_number()

		H1.<j1> = QuadraticField(d1).hilbert_class_field()

		H2.<j2> = K2.hilbert_class_field()

		# Find degree of intersection 
		H2poly = H2.absolute_polynomial()

		composite_extension_degree = H2poly.change_ring(H1).factor()[0][0].degree()

		# Degree of intersection is given by degree of H2 over degree of a factor of HCP2 over H1
		intersection_degree = (2*hd2 / composite_extension_degree)

		# Correction factor is given by intersection degree over twice class number of K2
		return intersection_degree / (2*hd2)


