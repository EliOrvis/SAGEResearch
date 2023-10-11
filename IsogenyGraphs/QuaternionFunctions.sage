#### This file contains functions that are useful for explicitly working with the quaternion algebra B_{p,infinity}.


### This function is from the paper Deuring For the People, by Jonathan, Jana, Lorenz, and someone else. Not my code.
##  Input - a quaternionic order in Bpinf for some prime p
##  Output - an integral left O-ideal of this order.
def random_ideal(O):
		r"""
		Sample an integral left O-ideal whose class is approximately uniform.
		"""
		Q = O.quaternion_algebra()
		i,j,k = Q.gens()
		p = ZZ(Q.discriminant())
		assert p.is_pseudoprime()

		I = O.unit_ideal()
		for it in reversed(range(55+2*ceil(log(p,2)))):   #TODO figure out good bound
				O = I.right_order()
				while True:
						beta = sum(randint(1,10**4)*a for a in O.basis())
						if gcd(beta.reduced_norm(), 4) == 2:
								break
				J = O*2 + O*beta
				assert J.norm() == 2
				I *= J
				if not it or I.norm() >= p**(ZZ(2)/3):
						# find an ideal of smaller norm in the class
						_,mn = I.quadratic_form().__pari__().qfminim(None,None,1)
						el = sum(ZZ(c)*g for c,g in zip(mn, I.basis()))
						I *= el.conjugate() / I.norm()
#        print(f'{it:4}', I)
		return I

#### This function returns the type number of Bpinf
##   Inputs - Quaternion Algebra with prime discriminant
##   Outpus - type number
def type_number(B):
  # Get the discriminant of B
  p = B.discriminant()

  # Validate that the discriminant is a prime
  if p.is_prime() == False:
    raise ValueError("Quaternion algebra must have prime discriminant.")

  # Use standard formulas to determine the type number of B
  base = floor(p/12)
  if p % 12 == 0:
    return base
  elif p % 12 in [5, 7]:
    return base + 1
  else:
    return base + 2

#### This function returns the type number of Bpinf
##   Inputs - ele - an element in a quaternion algebra;
##            basis - basis to write ele in terms of (with coefficients in Q)
##   Outpus - vector of coefficients of ele with respect to basis.
def quaternion_basis_rep(ele, basis):
  # Get parent quaternion algebra - this is why it is important
  # that ele is actually an element of the quaternion algebra
  # (e.g. pass QA(1) rather than 1 to write 1 in terms of a given basis).
  QA = ele.parent()

  # If the user provides a basis, use that. Otherwise, use O.basis().
  # This gives a list of vectors whose entries are the coefficients of 
  # each v in terms of the standard basis 1, i, j, k.
  B = [vector(QA(v)) for v in basis]

  # Check that the basis is actually linearly independent over Q
  if matrix(QQ, B).determinant() == 0:
    raise ValueError("%s is not a linearly independent set"%(basis))

  # Construct Q-vector space of rank 4 on basis B
  V = QQ**4
  W = V.submodule_with_basis(B)

  # Return coefficient tuple of ele in W
  return W.coordinate_vector(ele.coefficient_tuple())


#### This function returns representatives for the type set of Bpinf
#### WARNING: This function just picks random ideals until we get sufficiently many to have one order for the whole graph.
####          There must be a better way to do this.
##   Inputs - B, a quaternion algebra with prime discriminant
##   Outputs - list of maximal orders in B, none of which are contained in any others, of length equal to the size of the type set of B.
# def type_set_representatives(B):
#   # Validate discriminant condition on B
#   p = B.discriminant()
#   if p.is_prime() == False:
#     raise ValueError("Quaternion algebra must have prime discriminant.")

#   # Compute the type number of B, this will be used to know when to stop our loop 
#   type_num = type_number(B)

#   # Start with the maximal order that SAGE gives by default
#   O = B.maximal_order()
#   type_reps = [O]

#   # Until we have enough distinct maximal orders, keep finding new ones
#   while len(type_reps) != type_num:
#     # Take a random integral O ideal
#     I = random_ideal(O)

#     # Take the right order of I
#     OI = I.right_order()
    
#     # Validate that the right order of I has discriminant p (i.e. is maximal).
#     # This should always be true, but I can't seem to find a reference for this right now.
#     assert OI.discriminant() == p

#     # If we don't already have OI, we add it to type_reps
#     # CHECK THIS BIT WITH JAMES. DO WE KNOW THAT IN OUR CASE WE WILL ACTUALLY HAVE ISOMETRY RATHER THAN TWISTED SIMILARITY?
#     if not any([O.ternary_quadratic_form().is_rationally_isometric(OI.ternary_quadratic_form()) for O in type_reps]):
#       type_reps.append(OI)
#       O = OI
#     print(len(type_reps))

#   return type_reps


