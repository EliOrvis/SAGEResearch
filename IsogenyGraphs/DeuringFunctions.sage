##### This file contains methods for computing things related specifically
##### to the interplay between quaternion algebras and SS curves, i.e. the Deuring Correspondence

### Obtains the j-invariant of an Fp curve with endomorphism ring O
### Returns both j-invariants if the curve is defined over Fp^2
### NOTE: This is the most naive way I could think to do this, there is almost surely something better
## Inputs: O - order in a quaternion algebra ramified at a single prime p
## Output: complete list of supersingular j-invariants with endomorphism ring isomorphic to O
def j_from_order(O):
  p = O.quaternion_algebra().discriminant()
  assert p.is_prime()

  # Simply loop over discriminants, checking whether the discriminant embeds in O,
  # and then storing the possible values of j. Repeat until we have either 1 or 2 options.

  pos_js = []
  done = False
  d = -1
  while done == False:
    d = d - 1
    # Only consider fundamental discriminants with supersingular reductions modulo p
    if d.is_fundamental_discriminant() and legendre_symbol(d, p) != 1:
      pd = ZZ(d % 2)
      # If there are elements of O with discriminant d, then j is a root of H_D modulo p
      if len(quaternion_elements_by_minpoly(O, (pd - d)/4, trace = -pd)) > 0:
        H = hilbert_class_polynomial(d)
        roots = [tup[0] for tup in H.change_ring(GF(p^2)).roots()]
        if len(pos_js) == 0:
          pos_js = roots
        else:
          pos_js = [j for j in pos_js if j in roots]

      # Check if pos_js has only one or two elements, and stop if so
      if len(pos_js) == 1:
        # Verification that if we have only one root, it is in Fp
        assert pos_js[0]^p == pos_js[0]
        done = True


      elif len(pos_js) == 2:
        # If either is in Fp then we should continue, otherwise verify that we found conjugates and move on
        if pos_js[0]^p != pos_js[0] and pos_js[1]^p != pos_js[1]:        
          # Verification that if we have two roots, they are Frobenius conjugates
          assert pos_js[0]^p == pos_js[1]
          done = True

  return pos_js