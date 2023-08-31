#### This file contains implementations of formulae from Gross--Zagier's On Singular Moduli, which are helpful for studying isogeny graphs

## This helper function implements the delta function in Theorem 4.7 of On Singular Moduli
# Inputs: x - integer input to delta function; d - discriminant
# Outputs: delta(x)

def GZ_delta(x, d):  
	if x % -d == 0:
		return 2
	else:
		return 1

## This helper function implements the r_I function in Theorem 4.7 of On Singular Moduli
# Inputs: n - positive integer input to r_I; I - fractional ideal
# Outpus: r_I(n) - the number of representations of n by an ideal in the class of I
def repnum(n, I):

	# If n = 0, return 1 over the number of roots of unity in I.number_field()
	if n == 0:
		return 1/len(I.number_field().roots_of_unity())

	# Get the number field that I lives in 
	F = I.number_field()

	# Find all ideals of norm n in that field
	# Note that this is HIGHLY non-optimized, as it computes ALL ideals of
	#      norm up to level and then only considers the ideals of norm equal to level
	#      if performance every becomes an issue and/or before publication this should be fixed
	ideals_of_norm_n = F.ideals_of_bdd_norm(n)[n]

	# Get ideal class of I
	G = F.class_group()
	I_class = G(I)

	# Set a counter
	counter = 0

	# For each ideal of norm n, add one to counter if that ideal is in the class of I
	for ideal in ideals_of_norm_n: 
		if I_class == G(ideal):
			counter = counter + 1

	return counter


## This helper function computes the Artin isomorphism on a fractional ideal of a quadratic field
# Inputs: I - fractional ideal of K; hom - embedding of K into H
# Output: element of the Galois group of hom.codomain()
def artin_iso(I, hom):
	# Some validation:
	if hom.domain() != I.number_field():
		raise ValueError("Domain of hom must be number field in which I lives.")

	if hom.codomain().is_isomorphic(I.number_field().hilbert_class_field(names = 'b')) == False:
		raise ValueError("Hom must have codomain isomorphic to the HCF of I.")

	K = hom.domain()
	H = hom.codomain()

	clG = K.class_group()

	rep_prime = clG(I).representative_prime()

	H_prime = H.primes_above(hom(rep_prime))[0]

	GalG = H.galois_group()

	return GalG.artin_symbol(H_prime)

## This function does the computation of Theorem 4.7 from GZ on input B, A, level, p, d
def GZ_computation(B, A, level, p, d):
	# new var names are GZ notation, input variables are my notation 
	mp = level*(-d)

	mp_mod_l = ZZ(mp % p)

	max_n_multiple = floor(mp/p)

	total_val = 0

	for i in [0..max_n_multiple]:
		n_val = 0

		n = mp_mod_l + i*p

		if (mp - n) != 0:
			k_max = (mp - n).ord(p)

			for k in [1..k_max]:
				frac = (mp - n)*(1/(p^k))

				n_val = n_val + GZ_delta(n,d)*repnum(n, B.inverse())*repnum(frac, B*A^2)

			total_val = total_val + n_val

		else:
			total_val = total_val + GZ_delta(n,d)*repnum(n, B.inverse())*repnum(0, B*A^2)

	return total_val



## This function computes the P-order of phi_m(j_1, j_2), where j_1 and j_2 are assumed to be supersingular modulo p, and have CM by O_sqrt(-d) where d is a prime
# Inputs: G - isogeny graph; d - discriminant; 
#         B - fractional ideal of K taking real j-invariant w/ CM by d to j2 via CM/Galois action; 
#         level - level of modular curve; P - prime of H
# Outputs: ord_P(phi_m(j1,j2))  
# def GZ_thm4_7(G, d, B, level, P): 

# 	### Validate that the assumptions of the theorem are true.

# 	## Discriminant must be the negative of a prime
# 	if (-d).is_prime() == False:
# 		raise ValueError("Discriminant must be -p for a prime p.")

# 	## p must not split in the field Q(sqrt(d))	
# 	if legendre_symbol(d, G.prime()) == 1:
# 		raise ValueError("CM curves with discriminant d must be supersingular modulo p.")
	
# 	## There must not be an ideal of norm level in the ideal class of B
# 	## Note that this is HIGHLY non-optimized, as it computes ALL ideals of
# 	##      norm up to level and then only considers the ideals of norm equal to level
# 	##      if performance every becomes an issue and/or before publication this should be fixed

# 	# Get the number field that B lives in (this should be the HCF of d)
# 	# Note that I should be given as a *relative* ideal in H, so that we can recover the base field later
# 	H = B.number_field()
# 	if H.is_isomorphic(QuadraticField(d).hilbert_class_field()) == False:
# 		raise ValueError("B must be an ideal in the HCF of Q(sqrt(d)).")

# 	# Get all ideals of norm equal to level in H
# 	ideals_of_norm_level = H.ideals_of_bdd_norm(level)[level]

# 	# Get the class group of H
# 	G = H.class_group(); B_class = G(B)

# 	# For each ideal of norm level, check if it is equivalent to B
# 	for ideal in ideals_of_norm_level: 
# 		if B_class == G(ideal):
# 			raise ValueError("There cannot be any ideals of norm equal to level in the ideal class of B.")

# 	### Implement formula if no errors were raised

# 	## First, we need to find the ideal \mathfrak{a}^2. This can be done by fixing any
	