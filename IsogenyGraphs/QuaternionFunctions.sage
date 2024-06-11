#### This file contains functions that are useful for explicitly working with the quaternion algebra B_{p,infinity}.
import itertools
import pdb
from sage.matrix.matrix_integer_dense_hnf import hnf

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

  # This gives a list of vectors whose entries are the coefficients of 
  # each v in terms of the standard basis 1, i, j, k.
  base = [vector(QA(v)) for v in basis]

  # Check that the basis is actually linearly independent over Q
  if matrix(QQ, base).determinant() == 0:
    raise ValueError("%s is not a linearly independent set"%(basis))

  # Construct Q-vector space of rank 4 on basis B
  V = QQ**4
  W = V.submodule_with_basis(base)

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

  # Get trace 0 norm form - this will speed things up *dramatically* when trace is 0
  t0_norm_form, trace_0_basis = trace_0_norm_form(O)

  # Get all elements of norm equal to norm
  # This returns the vector coefficients for the standard basis

  # If trace is 0, then we use <trace_0_norm_form>
  if trace == 0:
    ele_vectors = t0_norm_form.representation_vector_list(norm + 1)[norm]
    eles = [sum([a*b for (a,b) in zip(tup,trace_0_basis)]) for tup in ele_vectors]
    return eles

  else:

    ele_vectors = norm_form.representation_vector_list(norm + 1)[norm]
    basis = O.basis()
    # This line constructs the actual elements
    eles = [sum([a*b for (a,b) in zip(tup,basis)]) for tup in ele_vectors]

    if trace != None:
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


#### This function returns a list of discriminants that embed (not necessarily optimally) into a quaternion order O
#### NOTE: This is a very simplistic approach. There is possibly a much better approach.
####       This also assumes (and requires) that the quaternion algebra is B p infinity for some prime p
##   Inputs - O, a quaternion order; ubound, an upper bound on the discriminant search; lbound (optional), lower bound
##   Outputs - List of discriminants
def quaternion_order_embedded_discs(O, ubound, lbound = 3):
  # Validation
  assert(O.quaternion_algebra().discriminant().is_prime())

  embedded_discs = []

  # Loop over all numbers from lbound to ubound, add to list if there are non-integer elements with the given discriminant 
  for d in [lbound..ubound]:
    # Check if -d is an imaginary quadratic discriminant
    if -d % 4 in [0,1]:
      # Find quaternion elements in O with min poly of discriminant -d
      eles = quaternion_elements_by_minpoly(O, (ZZ(d % 2) + d)/4, -ZZ(d % 2))
      # Add -d to the list if any of these elements are not in Z
      if any([x not in ZZ for x in eles]):
        embedded_discs.append(-d)

  # Return the norm form
  return embedded_discs

#### This function is a very naive approach to finding the number of optimal embeddings of a discriminant with prime-power 
####  conductor in a quaternion order of Bpinfinity
##   Inputs - O, a quaternion order; d - a fundamental discriminant; ell - a prime; n - non-negative integer;
##            optimal - optional parameter set by default to <True>, will count optimal embeddings if true
##   Outputs - number of optimal embeddings of <(ell^2)^i d> in <O> for each i up to n 

def naive_embedding_count(O, d, ell = 1, n = 0, optimal = True):
  # Make sure we are in Bpinfinit
  assert(O.quaternion_algebra().discriminant().is_prime())

  # This will store the number of optimal embeddings of <ell^n*d> for each n
  optimal_n_embeds = []
  # Loop up to n, and count the number of embeddings of discriminant <ell^n d>
  for i in [0..n]:
    disc = ((ell^2)^i)*d
    eles = quaternion_elements_by_minpoly(O, (ZZ(disc % 2) - disc)/4, -ZZ(disc % 2))
    # If we are interested in all embeddings, then we just need to divide the number of elements we found by 2
    if optimal == False:
      optimal_n_embeds.append(len(eles) / 2)
    else:
      # Number of optimal embeddings with conductor ell^(2n) is the number of conjugate pairs
      #  of elements with that discriminant, minus the number of all previous such pairs, since
      #  each previous pair gives one non-optimal pair.
      optimal_n_embeds.append(len(eles)/2 - sum(optimal_n_embeds))

  # If n = 0, then we are just looking for the number of embeddings of a particular order, so return that, rather than a list
  if n == 0:
    return optimal_n_embeds[0]

  # else, return the whole list
  return optimal_n_embeds

#### This function returns the number of square roots of <-ell> in maximal orders of B_{p, \infty},
#### up to automorphisms. This is done via an embedding count.
##   Inputs - p - a prime other than ell; ell - a prime other than p
##   Outputs - number of square roots of <-ell>
#    NOTE: This requires functions from my arithmetic_functions.sage file
def sqrt_ell_count(p, ell):
  # Use results from Chapter 30 in Voight to count the number of embeddings of sqrt(-ell), and
  # (1 + sqrt(-ell))/2 if necessary.
  if ell % 4 == 3:
    return imaginary_quadratic_order_class_number(-4*ell)*(1 - legendre_symbol(-ell,p)) + imaginary_quadratic_order_class_number(-ell)*(1 - legendre_symbol(-ell,p))
  else:
    return imaginary_quadratic_order_class_number(-4*ell)*(1 - legendre_symbol(-ell, p))

#### This function returns a basis for a quaternion order given a spanning set for that order.
##   Inputs - B: a quaternion algebra, span_set: a list of integral elements in <B> that span some order
##   Outputs - basis: a basis of the order spanned by the elements of <span_set>
def quaternion_order_basis_from_spanning_set(B, span_set):
  # Algorthim is pretty simple - just make the matrix whose columns are the spanning set with
  # respect to the default basis, then take the hermite normal form, and use the columns of that as
  # the coefficients of the basis.

  std_basis = B.basis()
  coeff_list = [quaternion_basis_rep(ele, std_basis) for ele in span_set]

  # Scale matrix to make columns into integers - this is necessary for HNF to work
  M = Matrix(QQ, coeff_list).transpose()
  M2 = MatrixSpace(ZZ, 4, len(span_set))(M.denominator()*M)

  # Take the zero entry in the hnf output b/c that is the actual matrix, and scale back, then take appropriate columns
  H = (1/M.denominator())*pari(M2).mathnf().sage()

  # Use the columns of the HNF to get a basis for the order spanned by <span_set>
  order_basis = [sum([a*b for (a,b) in zip(col, std_basis)]) for col in H.columns()]

  return order_basis

#### This function finds the intersection number of two optimal embeddings in a quaternion order
##   Inputs - O: a quaternion order in which <opt1> and <opt2> live;
##            opt1, opt2: elements of O which generate the quadratic orders we are intersecting (as rings);
##   Outputs - intersection_num: the intersection number of the embeddings <opt1>, <opt2>
def intersection_number_embeddings(O, opt1, opt2):

  # Validation
  assert(opt1 in O and opt2 in O)
  assert(O.quaternion_algebra().discriminant().is_prime())

  # Store p, B for future use
  B = O.quaternion_algebra()
  p = B.discriminant()

  # Initiate counter at 1, and set flag for when to stop
  intersection_num = 1
  still_intersecting = True
  while still_intersecting:
    # Make quaternion orders corresponding to phi(d1) + p^{n-1}O and same for phi(d2)
    opt1_O_span_set = flatten([[(p^intersection_num)*ele for ele in O.basis()],[B(1), opt1]]) 
    opt1_O = B.quaternion_order(quaternion_order_basis_from_spanning_set(B, opt1_O_span_set))

    opt2_O_span_set = flatten([[(p^intersection_num)*ele for ele in O.basis()],[B(1), opt2]])
    opt2_O = B.quaternion_order(quaternion_order_basis_from_spanning_set(B, opt2_O_span_set))

    # If these are still the same, we add one to the intersection count, which also increments the power of p
    if opt1_O == opt2_O:
      intersection_num += 1

    # Otherwise, we end the loop 
    else:
      still_intersecting = False

  return intersection_num

#### This function checks whether two pairs of quaternion elements are "simultaneously conjugate"
##   Inputs - O: quaternion order; pair1, pair2 - pairs of elements in O with the same min polys
##   Outputs - result: Boolean
def simul_conj(pair1, pair2):
  # Some validation
  assert(pair1[0].reduced_norm() == pair2[0].reduced_norm() and pair1[0].reduced_trace() == pair2[0].reduced_trace())
  assert(pair1[1].reduced_norm() == pair2[1].reduced_norm() and pair1[1].reduced_trace() == pair2[1].reduced_trace())

  # get units of O
  units = quaternion_elements_by_minpoly(O,1)

  # check for simultaneous conjugacy
  # NOTE: This is a little extraneous since we check both u and -u, but I'm too lazy to fix that right now
  for u in units:
    if pair1[0]*(u*pair2[0]*(1/u)) == (u*pair2[0]*(1/u))*pair1[0] and pair1[1]*(u*pair2[1]*(1/u)) == (u*pair2[1]*(1/u))*pair1[1]:
      return True

  # If we got through all units and did not return True, then we return False
  return False

#### This function intersection numbers for two discriminants in a quaternion order
##   Inputs - O: quaternion order; d1, d2 - embedded discriminants
##   Outputs - intersection_nums: list of intersection numbers for all pairs of embeddings
##             with the given discriminants, up to simultaneous conjugation.
def intersection_number_discriminants(O, d1, d2):
  # Get generators for the quadratic orders with discriminants d1, d2 in O
  pos_d1 = quaternion_elements_by_minpoly(O, (ZZ(d1 % 2) - d1)/4, -ZZ(d1 % 2))
  pos_d2 = quaternion_elements_by_minpoly(O, (ZZ(d2 % 2) - d2)/4, -ZZ(d2 % 2))

  # Make all pairs of generators
  gen_pairs = [[ele1, ele2] for ele1 in pos_d1 for ele2 in pos_d2]

  # Create equivalence relation function for simultaneous conjugation
  eq_classes = equal_classes(gen_pairs, simul_conj)

  # Get equivalence reps
  eq_reps = [cl[0] for cl in eq_classes]

  intersection_nums = []
  for pair in eq_reps:
    intersection_nums.append(intersection_number_embeddings(O, pair[0], pair[1]))

  return intersection_nums 


#### This function returns all suborders of index n in a given quaternion order
##   Inputs - O: quaternion order; n - positive integer
##   Outputs - suborders: list of suborders of index <n>
#    NOTE: This requires a function from my <arithmetic_functions.sage> file
def suborders_by_index(O, n):
  # Get the parent of O, and a basis
  QA = O.quaternion_algebra()
  O_basis = O.basis()

  # Every suborder is, in particular, a sublattice. So we find all of those with the right index
  sublattices = sublattices_by_index(n, len(O_basis))

  # For each sublattice, we check whether the columns give a suborder of <O>. If so, we keep it, otherwise, we move on
  suborders = []
  for M in sublattices:
    try: 
      quat_eles = [sum([coeff*basis_ele for (coeff, basis_ele) in zip(col, O_basis)]) for col in M.columns()]
      suborders.append(QA.quaternion_order(quat_eles))
    except:
      continue


  return suborders


