#### This file contains implementations of formulae from Gross--Zagier's On Singular Moduli, which are helpful for studying isogeny graphs

## This helper function implements the delta function in Theorem 4.7 of On Singular Moduli
# Inputs: x - integer input to delta function; d - discriminant
# Outputs: delta(x)

def GZ_delta(x, d):  
  if x % -d == 0:
    return 2
  else:
    return 1

## This helper function implements the epsilon function of On Singular Moduli
# Inputs: d1, d2 - discriminants; p - prime
# Outputs: epsilon(x)

def GZ_epsilon(d1,d2,p):  
  assert(p.is_prime())
  if gcd(d1, p) == 1:
    return legendre_symbol(d1, p)
  elif gcd(d2, p) == 1:
    return legendre_symbol(d2, p)
  else:
    raise ValueError("%s divides both %s and %s"%(p, d1, d2))

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

  # Find domain and codomain number fields - assumption is that I is an ideal in K and hom is a map to the HCF
  K = hom.domain()
  H = hom.codomain()

  # Construct the class group of G
  clG = K.class_group()

  # Pick a representative prime in order to compute the Artin map
  rep_prime = clG(I).representative_prime()

  # Pick a random prime of H above the representative prime in I
  H_prime = H.primes_above(hom(rep_prime))[0]

  # Smallest integer in H_prime is the prime that H_prime lies over
  p = H_prime.smallest_integer()

  # Validation to check whether we randomly picked one of the finitely many ramified primes.
  # Not sure yet how to make a second choice of representative prime if we run into this issue
  if p.divides(K.discriminant()):
	  raise ValueError("Rep prime in artin_iso is ramified, I have to fix this, why oh why?")

  # Assumes K is Galois so that we can pick whichever prime over p, this should always be true for me
  norm = K.primes_above(p)[0].norm()

  # Get Galois group of H - note that H has to be Galois for this to work, right now I don't have any validation guaranteeing
  GalG = H.galois_group()

  gens = H.ring_of_integers().ring_generators()	

  t = []

  # This is taken from source code, and just replaced "p" by "norm" so that we can get the Hilbert Symbol over other ground fields
  for s in GalG.decomposition_group(H_prime):
    w = [(s(g) - g**(norm)).valuation(H_prime) for g in gens]
    if min(w) >= 1:
      t.append(s)
  if len(t) > 1:
    raise ValueError("%s is ramified" % H_prime)
  return t[0]

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

## This function returns a table with the factorization data within a given valley
## Note that no validation of the necessary assumptions in OSM is implemented yet, this will come later.
def factorization_data_within_valley(level, p, d):
  # Construct relevant fields, polynomials, etc., making sure everything is Galois 
  K = QuadraticField(d)
  H.<gen> = K.hilbert_class_field('j').galois_closure()

  HCP = K.hilbert_class_polynomial()

  R.<t> = ZZ[]
  f = t^2 - d
  f_root  = f.change_ring(H).any_root()

  js = [root[0] for root in HCP.change_ring(H).roots()]



  # Get quadratic subfield and embedding in H - note that this assumes there is only ONE quadratic subfield of H
  # This assumption should be guaranteed in our case since the class number is odd, and therefore there is 
  # 	only one subgroup of Gal(H/K) with index 2. In a more general situation we would have to figure out which quadratic
  #   subfield is the correct one.
  QuadSub, QuadSub_embedding = [sub for sub in H.subfields() if sub[0].absolute_degree() == 2][0][:2]

  CLG = QuadSub.class_group()

  primes_in_H = H.primes_above(p)

  GalG = H.galois_group()

  counter = 0

  for j in js:
		
    for cl in CLG:
      B = cl.representative_prime()

      for prime in primes_in_H:

        # Find A for the particular prime
        # Again, everything is assuming there is only one of each, which is ok in our case as above
        tau = [g for g in GalG if g(j) == j and g.order() == 2][0]
        P = [prim for prim in primes_in_H if tau(prim) == prim][0]

        actors = [g for g in GalG if g(prime) == P]

        sigmaA = [g for g in actors if g(f_root) == f_root][0]

        A = [cl for cl in CLG if artin_iso(cl.representative_prime(), QuadSub_embedding) == sigmaA][0].representative_prime()


        if GZ_computation(B,A,level,p,d) > 0:
          counter += 1

  return counter


## This function is a MUCH faster version, which relies on a simplification when counting
##  over all primes P and conjugats j' of one fixed j. See RJ 9/5 or 9/6 for the details.
##  The output here should be the exact number of paths within the valley (of course with all)
##  of the usual assumptions in place, which still are not validated.

def count_edges_within_valley(d, p, l):

  # Initialize total counter
  count = 0

  # useful to reduce amount of writing
  ld = (-l)*d

  # Initialize class group
  clG = QuadraticField(d).class_group()

  # This could be sped up by iterating only over n in the right residue class
  for n in [1..ld]:
    if n % p != ld % p:
      count = count + 0
    else:
      # This case is necessary for the .ord call, could probably improve this with error handling
      if (ld - n) != 0:
        # This is the highest k value where the fraction is at least 1
        k_max = (ld - n).ord(p)

        for k in [1..k_max]:
          # Initialize counters
          rep_n_count = 0
          rep_frac_count = 0

          # Loop over all elements of the class group to count how many represent n, (ld-n)/p^k
          for cl in clG:
            # Pick an element of the ideal class to pass to rep_num
            rep_prime = cl.representative_prime()

            if repnum(n, rep_prime) > 0:
              rep_n_count = rep_n_count + 1
            if repnum((ld - n)*(1/p^k), rep_prime) > 0:
              rep_frac_count = rep_frac_count + 1

        # Add the product of the counts to our total
        count = count + rep_n_count*rep_frac_count

      else:
        # In this case, loop over just the class group
        rep_n_count = 0
        rep_0_count = 0
        for cl in CLG:
          rep_prime = cl.representative_prime()
          if repnum(n, rep_prime) > 0:
            rep_n_count = rep_n_count + 1
          if repnum(0, rep_prime) > 0:
            rep_0_count = rep_0_count + 1
        count = count + rep_n_count*rep_0_count 

  return count


#### This function the creates the order R(fraca, lambda) from Section 6 of Lauter Viray
##   Inputs - B, a quaternion algebra; fraca, a fractional ideal of the maximal order of K = Q(i),
##                                      where i is the first generator of <B>; 
##            lamb, an element of <K> satisfying the conditions of section 6.2 of LV;
##            embedding (optional), an element of <B> that can be the image of sqrt(d) in K. 
##                  If not given, then <i> is used as the target;
##            fracq (optional), fractional ideal of norm q containing <lamb>, 
##                  If not given, the program finds one (this could take a long time)
##   Outputs - O, a maximal order in <B>
##   NOTES: This is currently only implemented for the maximal order of <K>, and for the case where
##          the prime p (discriminant of <B>) is inert in <K>.
##          This also uses a function from my <QuaternionFunctions.sage> file
def Lauter_Viray_maximal_order(B, fraca, lamb, embedding = None, fracq = None):
  ### First, validation

  # <B> must have prime discriminant
  p = B.discriminant()
  assert(p.is_prime())

  # Discriminant of <K> must be the same as the field containing <fraca>
  K = fraca.number_field()
  tau = K.gens()[0]
  i,j,k = B.gens()
  d = K.discriminant()
  diff = K.different()
  assert(tau^2 == d)
  assert(i^2 == d)
  q = ZZ(-j^2/p)


  # Create <fracabar> and <fracainv> to validate conditions on <lamb>
  fracabar = K.ideal([gen.conjugate() for gen in fraca.gens()])
  fracainv = fraca.inverse()
  abarainv = fracabar*fracainv

  if fracq == None:

    # NOTE: right now, this is the expensive step. This can be made much faster by
    #       fiding the factorization of q, then finding the prime ideals of the right norms
    #       then trying the different combinations of them until we find the correct fracq
    qs = K.ideals_of_bdd_norm(q)[q]
    for qi in qs:
      if lamb in qi*(abarainv.inverse()):
        fracq = qi
        break

  else:

    # Validate that <fracq> and <fraca> are from the same field
    assert(fracq.number_field() == K)

    # Validate that <lamb> is in the right ideal
    assert(lamb in fracq*(abarainv.inverse()))

  # Validate that <lamb> is in q^{-1}(abar^{-1} a) for some q:
  assert(fracq)

  # Validate the condition on the norm of <lamb>:
  assert(lamb.norm() % d == -p*q % d)
  assert(lamb.is_integral())

  # Create homomorphism from <K> to <B> to transfer structure
  if embedding == None:
    phi = K.hom([i])
  else:
    assert(embedding^2 == d)
    phi = K.hom([embedding])

  # Create basis obtained as in my research journal entry on 12/4/2024
  beta1, beta2 = [phi(base_ele) for base_ele in (fracq.inverse()*diff.inverse()*abarainv).basis()] 
  basis = [1, (d + phi(tau))/2, phi(lamb)*beta1 + beta1*j, phi(lamb)*beta2 + beta2*j]
 
  O = B.quaternion_order(basis)

  return O

#### This function the creates the hom set Hom(E,E') given in Corollary 4.2.3 of Scullard's thesis
##   Inputs - B, a quaternion algebra; fracb, a fractional ideal of the maximal order of K = Q(i),
##                                      where i is the first generator of <B>; 
##            lamb (optional), an element of <K> satisfying the conditions of section 4.2.2 of Scullard's thesis;
##            embedding (optional), an element of <B> that can be the image of sqrt(d) in K. 
##                  If not given, then <i> is used as the target;
##            fracq (optional), fractional ideal of norm q containing <lamb>, 
##                  If not given, the program finds one (this could take a long time)
##   Outputs - an ideal in <B> representing the End(E)-module Hom(E, E')
##   NOTES: This is currently only implemented for the maximal order of <K>, and for the case where
##          the prime p (discriminant of <B>) is inert in <K>.
##          This also uses a function from my <QuaternionFunctions.sage> file
##          We also assume that the choice of fraca is (1), since this covers all cases anyway
def Scullard_Homset(B, fracb, lamb, fraca = None, embedding = None, fracq = None):
  ### First, validation

  # <B> must have prime discriminant
  p = B.discriminant()
  assert(p.is_prime())

  # Discriminant of <K> must be the same as the field containing <fraca>
  K = fracb.number_field()
  tau = K.gens()[0]
  i,j,k = B.gens()
  d = K.discriminant()
  diff = K.different()
  assert(tau^2 == d)
  assert(i^2 == d)
  q = ZZ(-j^2/p)


  # Create <fracabar> and <fracainv> to validate conditions on <lamb>
  fracbbar = K.ideal([gen.conjugate() for gen in fracb.gens()])
  fracbinv = fracb.inverse()
  bbarbinv = fracbbar*fracbinv

  # If the user does not supply fraca, make it trivial
  if fraca == None:
    fraca = K.maximal_order().fractional_ideal(1)

  fracabar = K.ideal([gen.conjugate() for gen in fraca.gens()])
  fracainv = fraca.inverse()
  abarainv = fracabar*fracainv

  if fracq == None:

    # NOTE: right now, this is the expensive step. This can be made much faster by
    #       fiding the factorization of q, then finding the prime ideals of the right norms
    #       then trying the different combinations of them until we find the correct fracq
    qs = K.ideals_of_bdd_norm(q)[q]
    for qi in qs:
      if lamb in qi*(bbarbinv.inverse())*(abarainv.inverse()):
        fracq = qi
        break

  else:

    # Validate that <fracq> and <fraca> are from the same field
    assert(fracq.number_field() == K)

    # Validate that <lamb> is in the right ideal
    assert(lamb in fracq*(bbarbinv.inverse())*(abarainv.inverse()))

  # Validate that <lamb> is in q^{-1}(abar^{-1} a) for some q:
  assert(fracq)

  # Validate the condition on the norm of <lamb>:
  assert(lamb.norm() % d == -p*q % d)
  assert(lamb.is_integral())

  # Create homomorphism from <K> to <B> to transfer structure
  if embedding == None:
    phi = K.hom([i])
  else:
    assert(embedding^2 == d)
    phi = K.hom([embedding])

  # Create basis obtained as in the case of maximal orders
  alpha1, alpha2 = [phi(base_ele) for base_ele in (diff.inverse()*fracb).basis()]
  beta1, beta2 = [phi(base_ele) for base_ele in (fracq.inverse()*diff.inverse()*abarainv*fracbbar).basis()] 
  basis = [alpha1, alpha2, phi(lamb)*beta1 + beta1*j, phi(lamb)*beta2 + beta2*j]
 
  I = B.ideal(basis)

  return I