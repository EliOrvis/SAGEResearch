#### This file contains miscellaneous arithmetic functions that are helpful for working with any code base.

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
# Note that the order need not be maximal, which is why this is more general than QuadraticField(d).class_number()

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
