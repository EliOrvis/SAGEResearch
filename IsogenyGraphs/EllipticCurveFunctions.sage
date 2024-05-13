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


### Function to obtain all cyclic endomorphisms of degree n (up to post-composition by an automorphism) of an elliptic curve
##  Inputs: E - Elliptic Curve; n - an integer 
##  Outputs: isos - a list of all endomorphisms of E of degree n 
def endomorphisms_of_degree_n(E, n):

  # Get the base field of the elliptic curve
  F = E.base_field()

  ## Step one is to extend the field of definition of E to a field where all of the n-torsion is defined
  #  This is done by taking the splitting field of the division polynomial
  div_poly = E.division_polynomial(n)
  K.<a> = div_poly.splitting_field()

  # Extend E to be defined over K
  E2 = E.base_extend(K)
      
  ## Step two is to compute all of the cyclic subgroups of order n in E
  #  This is done by finding a basis for the n-torsion, and then using standard methods
  
  # FIX THIS LATER
  n_points = [p for p in E2(0).division_points(n) if p.order() == n]

  ## Step three is to compute all of the isogenies for these subgroups
  #  We try to compute the isogenies to E, if this fails, we simply move on
  isos = []
  for p in n_points:
    try:
      iso = E2.isogeny(p, codomain = E2)
      if iso not in isos:
        isos.append(iso)
    except:
      continue

  # return the result
  return isos


### Function to obtain all cyclic isogenies of degree n (up to post-composition by an automorphism) out of an elliptic curve
##  Inputs: E - Elliptic Curve; n - an integer 
##  Outputs: isos - a list of all cyclic isogenies out of E of degree n 
def isogenies_of_degree_n(E, n):

  # Get the base field of the elliptic curve
  F = E.base_field()

  ## Step one is to extend the field of definition of E to a field where all of the n-torsion is defined
  #  This is done by taking the splitting field of the division polynomial
  K = E.division_field(n)

  # Extend E to be defined over K
  E2 = E.base_extend(K)
      
  ## Step two is to compute all of the cyclic subgroups of order n in E
  #  This is done by finding a basis for the n-torsion, and then using standard methods
  
  # FIX THIS LATER
  n_points = [p for p in E2(0).division_points(n) if p.order() == n]

  ## Step three is to compute all of the isogenies for these subgroups
  #  We try to compute the isogenies to E, if this fails, we simply move on
  isos = []
  for p in n_points:
    iso = E2.isogeny(p)
    if iso not in isos:
      isos.append(E2.isogeny(p))

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