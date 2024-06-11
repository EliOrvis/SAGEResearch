#### This file contains miscellaneous arithmetic functions that are helpful for working with any code base.
import itertools


# Function to give a random prime in a particular range with a particular modulo condition. Code from David Lowry-Duda

def random_prime_in_a_mod_b(lowerbound, upperbound, a, b, proof=False):
  """
  Returns a random prime in [lowerbound, upperbound]
  that is congruent to a mod b.
  """
  if not gcd(a, b) == 1:
      raise ValueError("a and b are not coprime")
  if not lowerbound < upperbound:
      raise ValueError("lowerbound is not smaller than upperbound")
  if proof:
      prime_test = is_prime
  else:
      prime_test = is_pseudoprime
  while True:
      n = randint(lowerbound // b, upperbound // b)
      p = n * b + a
      if p < lowerbound or p > upperbound:
          continue
      if prime_test(p):
          return p

# Function to get the coefficients x and y in ax - by == gcd(a,b) in Euclidean algorithm 
# Copied from https://www.rookieslab.com/posts/extended-euclid-algorithm-to-find-gcd-bezouts-coefficients-python-cpp-code

def extended_euclid_gcd(a, b):
    """
    Returns a list `result` of size 3 where:
    Referring to the equation ax + by = gcd(a, b)
        result[0] is gcd(a, b)
        result[1] is x
        result[2] is y 
    """
    s = 0; old_s = 1
    t = 1; old_t = 0
    r = b; old_r = a

    while r != 0:
        quotient = old_r//r 
        old_r, r = r, old_r - quotient*r
        old_s, s = s, old_s - quotient*s
        old_t, t = t, old_t - quotient*t

    # ADDED: negate everything if the gcd returned here is negative
    # This way we always have gcd(a,b) = <old_r> = <old_s*a + old_t*b>, and gcd(a,b) > 0
    if old_r < 0:
      old_r, old_s, old_t = -old_r, -old_s, -old_t
    return [old_r, old_s, old_t]

# Function to return the genus number of an imaginary quadratic discriminant.
# This is the number of two torsion elements in the class group. For proof, see Cox Proposition 3.11.
def genus_number_of_d(d):
  # Validation.
  if d >= 0 or d % 4 in [2,3]:
    raise ValueError("%s is not an imaginary quadratic discriminant."%(d))

  primes_dividing_d = ZZ(abs(d)).prime_divisors()

  odd_primes_dividing_d = [div for div in primes_dividing_d if div % 2 == 1]

  # r is the number of odd primes dividing d
  r = len(odd_primes_dividing_d)

  # mu is defined in cases
  if d % 4 == 1:
    mu = r
  else:
    n = -d/4
    if n % 8 == 0:
      mu = r + 2
    elif n % 8 == 3 or n % 8 == 7:
      mu = r
    else:
      mu = r + 1

  return 2^(mu - 1)

# This function returns the class number of an imaginary quadratic order from the discriminant
# Note that the order need not bex` maximal, which is why this is more general than QuadraticField(d).class_number()

def imaginary_quadratic_order_class_number(d):
  # Create the order of discriminant d
  pd = ZZ(d % 2)
  O.<z> = EquationOrder(x^2 + pd*x + (pd - d)/4)

  return O.class_number()


### This function returns the genus field of an imaginary quadratic field
##  Inputs: K - an imaginary quadratic field number field
##  Outputs: M - the genus field of K, as a relative number field extending K

def genus_field(K):

  # Create a polynomial ring so that we can use the variable x later
  R.<x> = ZZ[]

  # Get the discriminant of K and prime factors
  d = K.discriminant()
  prime_list = d.prime_divisors()

  # Adjust all the odd prime factors to be 1 mod 4:
  p_stars = [(-1)^((p - 1)/2)*p for p in prime_list if p != 2]

  # Set M = K so that we can extend K prime-by-prime
  M = K

  # Loop over the pstars and add them to M if necessary
  for pstar in p_stars:
    
    # If <x^2 + pstar> is already reducible over M, then we don't need to extend M,
    # so we only do anything if <x^2 + pstar> is irreudiclbe
    if (x^2 - pstar).change_ring(M).is_irreducible():
      
      # The try except here is to catch if the user already has a generator with one of these names
      try:
        M = M.extension(x^2 - pstar, 'a' + str(p_stars.index(pstar)))
      except:
        M = M.extension(x^2 - pstar, 'b' + str(p_stars.index(pstar)))

  # Minor error checking
  assert M.absolute_degree() == 2*genus_number_of_d(d)

  return M


### This function implements Dirichlet composition on positive definite quadratic forms.
### All of the theory can be found on pgs 39 - 40 of Cox's book.
##  Inputs: qf1, qf2 - two binary quadratic forms
##  Outputs: qf3 - the reduced form equivalent to the Dirichlet composition of these two forms

def Dirichlet_compose(qf1, qf2):

  # Some error checking:
  if qf1.discriminant() != qf2.discriminant():
    raise ValueError("Both forms must have the same discriminant.")

  d = qf1.discriminant()

  # To make sure that Dirichlet composition is defined on the level of forms, we have to change <qf2> to
  # a properly equivalent form whose leading coefficient is coprime to the leading coefficient of <qf1>

  # First, get the coefficients of the first form
  coeffs1 = qf1.polynomial().coefficients()
  a = coeffs1[0]
  if len(coeffs1) == 2:
    b = 0
    c = coeffs1[1]
  else:
    b = coeffs1[1]
    c = coeffs1[2]

  # Find a prime represented by qf2 that does not divide <a>
  # Start the search with the second smallest value represented by qf2
  #  in order to avoid excessivley long searches
  p = qf2(0,1).next_prime()
  found = False
  while found == False:
    if qf2.solve_integer(p) != None and p.divides(a) == False:
      found = True
      # These are stored as they will be used later
      r,s = qf2.solve_integer(p)
    else:
      p = p.next_prime()

  # Transform <qf2> to an equivalent form whose leading coefficient is p
  euclid_coeffs = extended_euclid_gcd(r,s)
  t, u = euclid_coeffs[1], -euclid_coeffs[2]

  # Coefficients of the second form are used for the transformation
  coeffs2 = qf2.polynomial().coefficients()
  ap_old = coeffs2[0]
  if len(coeffs2) == 2:
    bp_old = 0
    cp_old = coeffs2[1]
  else:
    bp_old = coeffs2[1]
    cp_old = coeffs2[2]

  # coefficients of new form can be obtained directly
  ap = qf2(r,s)
  bp = (2*ap_old*r*u + bp_old*r*t + bp_old*u*s + 2*cp_old*s*t)
  cp = qf2(u,t)

  # Compute B
  B = [x for x in [0..2*a*ap] if (x % (2*a) == b % (2*a) and x % (2*ap) == bp % (2*ap) and x^2 % (4*a*ap) == d % (4*a*ap))][0]

  # Return form
  return BinaryQF([a*ap, B, (B^2 - d)/(4*a*ap)]).reduced_form()

### This function returns the genera of binary quadratic forms with discriminant d.
##  Inputs: d - negative discriminant
##  Outputs: genera - list of lists, with each list one genus

def binary_qf_genera(d):

  # First, get all reduced binary quadratic forms of discriminant d
  qfs = BinaryQF_reduced_representatives(d)

  # Find the principal genus, which is the set of all squares
  principal_genus = []
  for qf in qfs:
    qfsq = Dirichlet_compose(qf,qf)
    if qfsq not in principal_genus:
      principal_genus.append(qfsq)

  # Now return list of cosets:
  genera = []
  # Sets are used to quickly check inclusion of lists without accounting for order
  genera_sets = []
  for qf in qfs:
    qf_coset = [Dirichlet_compose(qf, pqf) for pqf in principal_genus]
    qf_coset_set = set(qf_coset)
    if qf_coset_set not in genera_sets:
      genera.append(qf_coset)
      genera_sets.append(qf_coset_set)

  return genera


#### This function creates equivalence classes from a list of elements and a relation
##   Inputs: lst - a list of elements; R - a binary function that returns true or false on
##            two elements of <lst>
##`  Outputs: equal_classes - a list of lists that represent all equivalence classes under `R`
##   NOTE: There is little validation that `R` is an equivalence relation. Be careful.
def equal_classes(lst, R):
  # Storage for our equivalence classes 
  equal_class_list = []

  # Go over each item, and put it in the appropriate class
  for item in lst:
    # Get the existing equivalence class containing item, if it exists
    existing_class = [sub_lst for sub_lst in equal_class_list if any([R(item, b) for b in sub_lst])]
    # If there is such a class, validate that there is only one, and add item
    if len(existing_class) > 0:
      assert(len(existing_class) == 1)
      existing_class[0].append(item)
    # Otherwise, make a new class with item in it
    else:
      equal_class_list.append([item])

  return equal_class_list


#### This function finds matrix representative of the sublattices of index n for a generic lattice of dimension d
##   Inputs - d, n: positive integers
##   Outputs - sublattices: matrix representatives of all sublattice of index n in a generic lattice of dimension d
#    NOTES: <n> is currently restricted to a prime power
#           This is based on the theory described here: https://journals.iucr.org/a/issues/1997/06/00/au0106/au0106.pdf
def sublattices_by_index(n, d):
  # Validate that n is a prime power
  assert(n.is_prime_power())

  ### Enumerate all possible diagonals
  diagonals = []

  # Get valuation and prime of n
  prime, valuation = n.factor()[0]


  # Get partitions of the valuation with length at most d
  val_parts = Partitions(valuation, max_length = d).list()

  # For each partition, if the length is less than d, pad with zeros and add all permutations to diagonals
  for part in val_parts:
    if len(part) < d:
      part = part + [0]*(d - len(part))
    diagonals = diagonals + list(itertools.permutations(part))

  # Get only unique values
  diagonals = list(set(diagonals))

  # Alter diagonals to have the correct power of <prime>
  diagonals = [[prime^i for i in ele] for ele in diagonals]

  ### Construct sets of all possible rows
  sublattices = []
  for diag in diagonals:
    pos_rows = []
    for i in [0..d-1]:
      # Create a list of possible values for each entry in the row
      entry_limits = [[0] for j in [0..i-1]] + [[diag[i]]] + [[0..entry - 1] for entry in diag[i+1:]]
      # The possibilities for row i for this diagonal are the product of all the entry limits
      pos_rows.append(list(itertools.product(*entry_limits)))

    # Add all matrices made from one choice of each possible row to the sublattices
    sublattices = sublattices + list(itertools.product(*pos_rows))

  return [Matrix(sl) for sl in sublattices]

