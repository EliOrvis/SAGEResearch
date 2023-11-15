#### This file contains functions that are useful for explicitly working with the quaternion algebra B_{p,infinity}.
import itertools
import pdb

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

#### This function returns the expanded type number of Bpinf (i.e. the number of SS curves)
##   Inputs - Quaternion Algebra with prime discriminant
##   Outpus - expanded type number
def expanded_type_number(B):
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

#### This function returns the norm form on the sublattice of trace 0 elements for a quaternion order
##   Inputs - O, a quaternion order
##   Outpus - list: first element is norm form as a sage quadratic form object;second is the basis for the trace 0 elements 
def trace_0_norm_form(O):
  # Get the standard 1,i,j,k basis for the QA and re-order to i,j,k,1
  B = O.quaternion_algebra()
  std_basis = B.basis()

  reordered_basis = (std_basis[1],std_basis[2],std_basis[3],std_basis[0])

  # The next goal is to find the HNF of the basis elements with respect to the reordered basis
  # The first three rows & columns of this will give Z-generators for O^0

  # Make the matrix whose columns are the generators for O with respect to the reordered basis
  row1 = quaternion_basis_rep(O.basis()[0], reordered_basis)
  row2 = quaternion_basis_rep(O.basis()[1], reordered_basis) 
  row3 = quaternion_basis_rep(O.basis()[2], reordered_basis)
  row4 = quaternion_basis_rep(O.basis()[3], reordered_basis)

  # We have to have this be over Z, so we clear denominators:
  denom = lcm([x.denominator() for x in flatten([list(row1), list(row2), list(row3), list(row4)])])
  M = matrix(ZZ, [denom*row1, denom*row2, denom*row3, denom*row4])

  # Get HNF from pari, rescaling to fix the extra factor of denom
  HNF = (1/denom)*pari(M.transpose()).mathnf().sage()

  # Create the quaternion elements that correspond to these columns
  reordered_reps = HNF[:3,:3].columns()
  O0_basis = tuple([sum([a*b for (a,b) in zip(reordered_basis[:3],column)]) for column in reordered_reps])

  # Find coefficients for the norm form on this lattice
  # Note that Sage uses lexicographic ordering, so coeff1 is on x0^2, 2 is on x0*x1 etc.
  coeff1 = O0_basis[0].reduced_norm()
  coeff2 = (O0_basis[0]*(O0_basis[1].conjugate())).reduced_trace()
  coeff3 = (O0_basis[0]*(O0_basis[2].conjugate())).reduced_trace()
  coeff4 = O0_basis[1].reduced_norm()
  coeff5 = (O0_basis[1]*(O0_basis[2].conjugate())).reduced_trace()
  coeff6 = O0_basis[2].reduced_norm()

  # Return the quadratic form with these coefficients
  return [QuadraticForm(ZZ, 3, [coeff1, coeff2, coeff3, coeff4, coeff5, coeff6]), O0_basis]




#### This function returns the representation of a quaternion element in a given base
##   Inputs - ele - an element in a quaternion algebra;
##            basis - basis to write ele in terms of (with coefficients in Q)
##   Outputs - vector of coefficients of ele with respect to basis.
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

#### This function returns the subspace of elements (as vectors) wrt the standard base
####   conjugating ele1 to ele2
##   Inputs - ele1, ele2 - elements in a quaternion algebra;
##   Outpus - vector space of elements (with respect to the standard basis) conjugating ele1 to ele2.
#     Note: conjugating here means that x ele1 = ele2 x, not ele1 x = x ele2.
def conjugation_subspace(ele1, ele2):
  # Validate that ele1 and ele2 are elements of the same quaternion algebra
  assert ele1.parent() == ele2.parent()

  # Get coefficients of ele1 and ele2 wrt to standard basis as variables
  c0, c1, c2, c3 = ele1.coefficient_tuple()
  d0, d1, d2, d3 = ele2.coefficient_tuple()

  # Get invariants of parent QA
  s1, s2 = ele1.parent().invariants()

  # Create matrix whose solution space is the coefficients of elements
  # x = a + bi + cj + dij such that x ele1 = ele2 x. This was found by hand,
  # and the calculation is easy, if tedious.
  row1 = [c0 - d0, c1*s1 - s1*d1, c2*s2 - s2*d2, d3*s1*s2 - c3*s1*s2]
  row2 = [c1 - d1, c0 - d0, -d3*s2 - c3*s2, c2*s2 + d2*s2]
  row3 = [c2 - d2, c3*s1 + d3*s1, c0 - d0, -c1*s1 - d1*s1]
  row4 = [c3 - d3, c2 + d2, -d1 - c1, c0 - d0]
  # The first entry is the base field. For some reason this doesn't seem easily 
  # accessible directly from the QA, and we have to take the underlying VS first.
  M = matrix(ele1.parent().vector_space().base_field(), [row1, row2, row3, row4])

  return M.kernel()

#### This function returns the norm form of a quaternion algebra
##   Inputs - B, a quaternion algebra;
##   Outputs - Quadratic Form that is the norm form for the quaternion algebra
def quaternion_norm_form(B):
  # Get invariants of B
  a, b = B.invariants()

  # Return the norm form
  return DiagonalQuadraticForm(B.vector_space().base_field(), [1, -a, -b, a*b])

#### This function returns all elements of given norm (and trace, optional) from an order 
#### in a definite quaternion algebra over Q
##   Inputs - O, an order in a definite quaternion algebra; norm, integer to be the reduced norm
##            trace, integer to be the reduced trace (optional)
##   Outputs - list of elements of B with given norm (and trace, optional)
def quaternion_elements_by_minpoly(O, norm, trace = None):
  # Validate that B is definite over Q
  B = O.quaternion_algebra()
  assert B.vector_space().base_field() == Rationals()
  assert len(B.ramified_primes()) % 2 == 1
  
  # Get quaternion norm form
  norm_form = O.quadratic_form()

  # Get all elements of norm equal to norm
  # This returns the vector coefficients for the standard basis
  ele_vectors = norm_form.representation_vector_list(norm + 1)[norm]
  basis = O.basis()
  # This line constructs the actual elements
  eles = [sum([a*b for (a,b) in zip(tup,basis)]) for tup in ele_vectors]

  if trace:
    # Return only those elements with specified reduced trace
    return [ele for ele in eles if ele.reduced_trace() == trace]
  else:
    return eles

#### This function determines whether two quaternion orders are isomorphic, and (optionally)
#### returns the linear subspace (wrt the standard basis) of elements that conjugate one to the other
##   Inputs - O1 & O2, orders in a definite quaternion algebra; explicit, optional flag, default is false
##            if true, returns the vector space described above, rather than a boolean
##   Outputs - boolean, or vector space of elements (wrt the standard basis) conjugating O1 to O2
def quaternion_order_isomorphisms(O1, O2, explicit = False):
  # Validate that O1 and O2 are from the same QA and that algebra is definite over Q
  B = O1.quaternion_algebra()
  assert B.vector_space().base_field() == Rationals()
  assert len(B.ramified_primes()) % 2 == 1
  assert B == O2.quaternion_algebra()
  
  # Find a short basis for the trace 0 elements in O1
  O1_tr0_nf, O1_nf_basis = trace_0_norm_form(O1)
  # This gives the coefficients of short vectors with respect to the basis for tr0 elements
  O1_short_basis_reps = O1_tr0_nf.basis_of_short_vectors()
  # We now convert those into quaternion elements

  O1_short_basis = tuple([sum([a*b for (a,b) in zip(rep, O1_nf_basis)]) for rep in O1_short_basis_reps])

  # Find the norm form for the trace 0 elements in O2
  O2_tr0_nf, O2_nf_basis = trace_0_norm_form(O2)

  # Check whether O2 has the same number of elements as O1 of the two shortest norms in O1
  # Along the way, we compute the actual elements so that if these counts do agree,
  # we don't have to re-compute them.
  # We only do this for two elements of O1 because if we get two L.I. elements that both conjugate into O2,
  # then the intersection of the conjugating subspaces will be dimension 0 or 1.
  O1_norms = [b.reduced_norm() for b in O1_short_basis[:2]]
  O1_short_vectors = [O1_tr0_nf.representation_vector_list(nrd + 1)[nrd] for nrd in O1_norms]
  O2_short_vectors = [O2_tr0_nf.representation_vector_list(nrd + 1)[nrd] for nrd in O1_norms]

  if any([len(O1_short_vectors[i]) != len(O2_short_vectors[i]) for i in [0..len(O1_norms)-1]]):
    # The condition above means there are a different number of norm N vectors in O2 than O1 for some N, so they are non-isomorphic
    return False

  # For the next step, we will need the actual elements for the short vectors in O2
  O2_short_elements = [[sum(a*b for (a,b) in zip(vector, O2_nf_basis)) for vector in vector_list] for vector_list in O2_short_vectors]

  # If the above is true, we need to know whether we can simultaneously conjugate the first two elements of the short
  # basis for O1 to two elements of O2.
  # For each conjugate in O2, get the vector space of elements that conjugate the related basis element to that element.
  conjugate_spaces_v1 = [conjugation_subspace(O1_short_basis[0], ele) for ele in O2_short_elements[0]]
  conjugate_spaces_v2 = [conjugation_subspace(O1_short_basis[1], ele) for ele in O2_short_elements[1]]
  conjugate_spaces = [conjugate_spaces_v1, conjugate_spaces_v2]

  # Iterate over one choice of each subspace, and if any such choice gives a nonzero intersection over all the subspaces, store the subspace
  valid_spaces = []
  for spaces in itertools.product(*conjugate_spaces):
    base_space = spaces[0]
    for it in range(len(spaces)): 
      base_space = base_space.intersection(spaces[it])
    if base_space.dimension() > 0:
      valid_spaces.append(base_space)
  
  # If there are no choices where the intersection of all spaces has dimension at least 1, return False
  if len(valid_spaces) == 0:
    return False

  # Otherwise, check whether a non-zero element of any of the valid subspaces conjugates all of O1 into O2
  elif explicit == False:
    # Get one representative quaternion element from each valid subspace
    valid_conjugate_reps = [space.basis()[0] for space in valid_spaces]
    valid_conjugate_eles = [sum([a*b for (a,b) in zip(vector_rep, B.basis())]) for vector_rep in valid_conjugate_reps]
    for ele in valid_conjugate_eles:
      # If, for any of these elements, all the conjugates of a basis for all of O1 is in O2, return True
      if all([ele*basis_ele in O2*ele for basis_ele in O1.basis()]):
        return True
    # If none of them work, then nothing will work, so we return False
    return False

  # Finally, if explicit is true, find the vector space generated by all valid subspaces
  # THIS STILL NEEDS TO BE IMPLEMENTED


### This function returns representatives for the type set of Bpinf
#   Inputs - B, a quaternion algebra with prime discriminant
#   Outputs - list of maximal orders in B, none of which are contained in any others, of length equal to the size of the type set of B.
def type_set_representatives(B):
  # Validate discriminant condition on B
  p = B.discriminant()
  if p.is_prime() == False:
    raise ValueError("Quaternion algebra must have prime discriminant.")

  # Compute the expanded type number of B, this will be used to know when to stop our loop 
  type_num = expanded_type_number(B)

  # Start with the maximal order that SAGE gives by default
  O = B.maximal_order()
  type_reps = [O]

  Frob_eles_in_O = quaternion_elements_by_minpoly(O, p, trace = 0)

  # Until we have enough distinct maximal orders, keep finding new ones
  # We start with either counting O as 2 (corresponds to Fp2 vertices) or 1 (corresponds to an Fp vertex)
  if len(Frob_eles_in_O) > 0:
    counter = 1
  else:
    counter = 2

  while counter != type_num:
    # Take a random integral O ideal
    I = random_ideal(O)
    print(counter)

    # Take the right order of I
    OI = I.right_order()
    
    # Validate that the right order of I has discriminant p (i.e. is maximal).
    # This should always be true, but I can't seem to find a reference for this right now.
    assert OI.discriminant() == p

    # If OI is not isomorphic to any of the orders we already have, append it
    if all([quaternion_order_isomorphisms(OI, order) == False for order in type_reps]):
      type_reps.append(OI)
      Frob_eles_in_OI = quaternion_elements_by_minpoly(OI, p, trace = 0)
      if len(Frob_eles_in_OI) > 0:
        counter += 1
      else:
        counter += 2

  return type_reps

