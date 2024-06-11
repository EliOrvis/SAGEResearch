###### Functions for working with the Ihara zeta functions of supersingular isogeny graphs.
###### Note: there are many dependencies to other files of mine.

### Function to compute the Ihara zeta function of a supersingular isogeny graph
### For explanations of why the formula implemented below works, see the upcoming paper by myself, TM, GS, JBL, and LZ
##  Inputs: G - supersingular isogeny graph object
##  Outputs: zeta - rational polynomial in one variable that is the Ihara zeta function
def ihara_zeta(G):
  # Construct rational function space
  R.<u> = ZZ[]
  S = R.fraction_field()

  # Get adjacency matrix and isogeny degree from G
  A = G.adjacency_matrix()
  ell = G.isogeny_degree()
  p = G.prime()

  # Get the number of edges of G
  n_edges = len(G.vertices())*(ell + 1)

  # Get the number of self-dual loops of G (this is done using the formula, not my code for finding all self-dual loops)
  if ell % 4 == 3:
    n_embeddings = (1/2)*(imaginary_quadratic_order_class_number(-4*ell) + imaginary_quadratic_order_class_number(-ell))*(1 - legendre_symbol(-ell, p))
  else:
    n_embeddings = (1/2)*imaginary_quadratic_order_class_number(-4*ell)*(1 - legendre_symbol(-ell, p))

  if ell == 2 and p % 4 == 3:
    n_sd = n_embeddings + 1
  else:
    n_sd = n_embeddings 

  # Implement formula 
  if p % 12 == 1:
    eps = 0
  elif p % 12 == 5:
    eps = 2*(ell + 1 - len([loop for loop in G.loops() if loop[0] == 0]))/3
  elif p % 12 == 7:
    eps = (ell + 1 - len([loop for loop in G.loops() if loop[0] == 1728]))/2
  elif p % 12 == 11:
    if 0 in G.neighbors(1728 % p):
      raise NotImplementedError("Ihara zeta only implemented when 0 and 1728 are not neighbors in G.")
    else:
      eps = 2*(ell + 1 - len([loop for loop in G.loops() if loop[0] == 0]))/3 + (ell + 1 - len([loop for loop in G.loops() if loop[0] == 1728]))/2

  num_exp = ((2 - ell - 1)*len(G.vertices()) + n_sd + eps)/2

  return (1 - u^2)^num_exp / (((1 - A*u + ell*u^2).determinant())*(1 + u)^n_sd)

### Function to directly find the non-backtracking cycles of a given length in the SSI graph
##  Inputs:  G - SSI graph object; n - length of cycles to count; isos (Optional) - list of isogenies in G
##  Outputs: cycles - non-backtracking cycles of length n
#   NOTES:   (1) Currently only implemented for n = 2. (2) Computing the isogenies and duals is the bottleneck, so passing isos
#            is SUBSTANTIALLY faster.
def get_non_backtracking_cycles(G, n = 2, isos = None):
  # First, get the isogenies if not already given
  if isos == None:
    isos = get_isogenies(G)

  # Verify that the number of isogenies equals the number of edges.
  # This is just a simple check in case someone passes the wrong isogeny set for G
  assert(len(G.edges()) == len(isos))

  # Temporary until I find time to implement for general n
  assert(n == 2)

  # For each vertex count the number of paths of length two starting and ending at that vertex
  cycles = []
  for v in G.vertices():
    for iso in isos:
      if iso.domain().j_invariant() == v:
        # Get all isogenies back from the codomain to the domain
        for isoback in isos:
          if isoback.codomain() == iso.domain() and isoback.domain() == iso.codomain():
            # Add one only if isoback is not the dual of iso and iso is not the dual of isoback
            # The first step is checking for backtracking, the second is checking for tails
            if all([isoback.post_compose(aut) != iso.dual() for aut in iso.domain().automorphisms()]):
              if all([iso.post_compose(aut) != isoback.dual() for aut in iso.codomain().automorphisms()]):
                cycles.append((iso, isoback))

  return cycles


### Function to return the matrix of the non-backtracking operator on edges of an isogeny graph
##  Inputs: G - Isogeny graph
##  Outputs: M - matrix of the action of the non-backtracking operator
def get_non_backtracking_operator(G):
    # Get all isogenies in G
    isos = get_isogenies(G)

    # Make matrix
    M = Matrix(ZZ, len(isos), lambda i,j : 1 if (isos[j].domain() == isos[i].codomain() and not any([isos[i].dual() == isos[j].post_compose(aut) for aut in isos[j].codomain().automorphisms()])) else 0)

    return M