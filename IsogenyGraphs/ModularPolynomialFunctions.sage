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

	# For each prime and j-invariant, we fix a map from the residue field at p to the finite field Fp^2
	# We save (j, p) if j reduces to v under this particular map.
	# Note: Since there are two maps from ResField to GF(p^2), choosing instead to save all pairs simply gives
	##      twice as many pairs (j,P), where we have now included both (j,P) and (sigma(j),sigma(P))
	for p in HCF_primes:
		ResField = p.residue_field()
		Homs = Hom(ResField, GF(G.prime()^2))
		Fixed_hom = Homs[0]
		for j in j_invars:
			if Fixed_hom(j) == v:
				lifts.append((j, p))

	return lifts