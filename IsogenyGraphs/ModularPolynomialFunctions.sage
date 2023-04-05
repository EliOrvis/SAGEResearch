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

	v_min_poly = minimal_polynomial(v)

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