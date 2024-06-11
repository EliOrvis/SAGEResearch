from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

### Function to obtain the isogeny phi' from a separable isogeny phi, known to be m times phi'. 
### This is based on Proposition 2.6 in McMurdy's, Explicit Representation of the Endomorphism Rings of Supersingular Elliptic Curves (see McMurdy for definitions of terms in comments)
## Inputs: phi - an isogeny; m - an integer to divide by
## Outputs: phi' - an isogeny  
def isogeny_division_by_m(phi, m):
  # Initialize domain and codomain curves
  E1 = phi.domain()
  E2 = phi.codomain()

  # Get F(x) and G(x) 
  F = phi.rational_maps()[0]
  G = phi.rational_maps()[1]

  # Get c_F, P(x), Q(x)
  ## First we find cF, by getting leading coefficents of numerator and denominator of F(x)
  top_lead = F.numerator().coefficents()[0]
  bot_lead = F.denominator().coefficents()[0]
  cF = top_lead / bot_lead

  # To find P(x), just multiply the numerator by 1/leading coefficient
  P = F.numerator() / top_lead

  # To find Q(x), divide out either the kernel polynomial, or the kernel polynomial squared
  ker_poly = E1.division_polynomial(m)
  if m == 2:
  	Q = (F.denominator()/bot_lead).quo_rem(ker_poly)[0]
  else:
  	Q = (F.denominator()/bot_lead).quo_rem(ker_poly^2)[0]

  # get x-function of multiplication-by-m on both sides
  m1 = E1.scalar_multiplication(m).rational_maps()[0]
  m2 = E2.scalar_multiplication(m).rational_maps()[1]

  # construct p(x) and q(x)
  return(P,Q)


### Function to obtain all cyclic endomorphisms of degree n (up to post-composition by an automorphism) of an elliptic curve over a finite field.
##  Inputs: E - Elliptic Curve; n - an integer 
##  Outputs: isos - a list of all endomorphisms of E of degree n 
def endomorphisms_of_degree_n(E, n):

  # Get the base field of the elliptic curve
  F = E.base_field()

  ## Step one is to extend the field of definition of E to a field where all of the n-torsion is defined
  #  This is done by taking the splitting field of the division polynomial
  div_poly = E.division_polynomial(n)
  K.<a> = div_poly.splitting_field()

  # Extend E to be defined over a quadratic extension of K
  # NOTE: This is where we use that the field is finite, since in this case all quadratic extensions are the same.
  #       The quadratic extension is necessary in the first place to obtain the y-values.
  E2 = E.base_extend(K.extension(2))
    
  ## Step two is to compute all of the cyclic subgroups of order n in E
  #  This is done by finding a basis for the n-torsion, and then using standard methods
  
  # Find the points of order <n>
  n_points = [p for p in E2(0).division_points(n) if p.order() == n]

  # Find the number of subgroups of order <n>
  n_subgroups = len([(a,b) for a in [0..n-1] for b in [0..n-1] if gcd(n,gcd(a,b)) == 1])/euler_phi(n)

  # Find generators for the subgroups of order <n> by adding points to the list until we have enough
  subgroup_gens = []
  total_subgroup_elements = []
  i = 0
  while len(subgroup_gens) < n_subgroups:
    if n_points[i] not in flatten(total_subgroup_elements):
      subgroup_gens.append(n_points[i])
      total_subgroup_elements.append([k*n_points[i] for k in [0..n-1]])
    i += 1

  ## Step three is to compute all of the isogenies for these subgroups
  #  We try to compute the isogenies to E, if this fails, we simply move on
  isos = []
  for p in subgroup_gens:
    try:
      iso = E2.isogeny(p, codomain = E2)
      if iso not in isos:
        isos.append(iso)
    except:
      continue

  # return the result
  return isos


### Function to obtain all cyclic isogenies of degree n (up to post-composition by an automorphism) out of an elliptic curve
##  Inputs: E - Elliptic Curve; n - an integer; model_set - a set of models in which there exists a model of all the codomains 
##  Outputs: isos - a list of all cyclic isogenies out of E of degree n 
def isogenies_of_degree_n(E, n, model_set = None):

  # Get the base field of the elliptic curve
  F = E.base_field()

  ## Step one is to extend the field of definition of E to a field where all of the n-torsion is defined
  #  This is done by taking the splitting field of the division polynomial
  K = E.division_field(n)

  # Extend E to be defined over K
  E2 = E.base_extend(K)

  # If we have explicit models, we need to change the field of definition for the models
  if model_set != None:
    model_set = [E.base_extend(K) for E in model_set]
      
  ## Step two is to compute all of the cyclic subgroups of order n in E
  #  This is done by finding a basis for the n-torsion, and then using standard methods
  
  # FIX THIS LATER
  n_points = [p for p in E2(0).division_points(n) if p.order() == n]

  ## Step three is to compute all of the isogenies for these subgroups
  #  We try to compute the isogenies to E, if this fails, we simply move on
  isos = []
  for point in n_points:
    # if the user supplies a set of models to consider, create isogenies into those models
    if model_set:
      found = False
      for model in model_set:
        try:
          iso = E2.isogeny(point, codomain = model)
          found = True
        except:
          if found == True:
            continue
          else:
            iso = None
      # After the loop, assert that at least one model worked. Otherwise, we got a bad model set
      assert(iso != None)
    # If no model set is supplied, just make the isogeny to wherever
    else:
      iso = E2.isogeny(point)
    if iso not in isos:
      isos.append(iso)

  # return the result
  return isos


### Function to check whether an endomorphism of a supersingular elliptic curve is self-dual
##  Inputs: phi - an endomorphism of a supersingular elliptic curve
##  Outputs: is_self_dual - boolean
def is_self_dual_ss(phi):
  # First, some input verification
  assert(phi.codomain().j_invariant() == phi.domain().j_invariant()) # Check that phi is really an endomorphism
  assert(phi.domain().is_supersingular()) # Check that the EC is ss (this is necessary to guarantee later steps work)

  # Extend to a field where all the endomorphisms of the domain/codomain are defined.
  # This guarantees that there is a Weierstrass isomorphism between the domain and codomain 
  codom = phi.codomain()
  dom = phi.domain()
  F = phi.base_ring()
  ker_poly = phi.kernel_polynomial()
  if dom.j_invariant() == F(0):
    Fext = F.extension(3)
  elif dom.j_invariant() == F(1728):
    Fext = F.extension(2)
  else:
    Fext = F

  domext = dom.base_extend(Fext)
  codomext = codom.base_extend(Fext)
  phiext = domext.isogeny(ker_poly, codomain = codomext)

  # Get a Weierstrass isomorphism from the codomain to the domain
  iso = WeierstrassIsomorphism(codomext, None, domext)
  phinew = phiext.post_compose(iso)

  # Check if phi is self-dual (up to post-composition by an automorphism)
  auts = domext.automorphisms()
  if any([phinew == phinew.dual().post_compose(aut) for aut in auts]):
    is_self_dual = True
  else:
    is_self_dual = False

  return is_self_dual