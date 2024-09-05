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
      eps = 2*(ell + 1 - len([loop for loop in G.loops() if loop[0] == 0]))/3 + (ell + 1 - len([loop for loop in G.loops() if loop[0] == 1728]))/2 
    else:
      eps = 2*(ell + 1 - len([loop for loop in G.loops() if loop[0] == 0]))/3 + (ell + 1 - len([loop for loop in G.loops() if loop[0] == 1728]))/2

  num_exp = ((2 - ell - 1)*len(G.vertices()) + n_sd + eps)/2

  return (1 - u^2)^num_exp / (((1 - A*u + ell*u^2).determinant())*(1 + u)^n_sd)

### Function to directly find the non-backtracking cycles of a given length in the SSI graph
##  Inputs:  G - SSI graph object; n - length of cycles to count; isos (Optional) - list of isogenies in G
##  Outputs: tailless_cycles - non-backtracking, tailless cycles of length n
##  NOTES:   This is exponential in <n>. Tread carefully.
def get_non_backtracking_cycles(G, n, isos = None):
  # First, get the isogenies if not already given
  if isos == None:
    isos = get_isogenies(G)

  # Verify that the number of isogenies equals the number of edges.
  # This is just a simple check in case someone passes the wrong isogeny set for G
  assert(len(G.edges()) == len(isos))
  assert(n > 1)

  # For each vertex, find non-backtracking paths of length <n> beginning and ending at that vertex
  nb_paths_of_length_n = []
  for v in G.vertices():
    # Counter for how many steps we need to take.
    steps_left = n - 1
    # List of current paths
    current_paths = [[iso] for iso in isos if iso.domain().j_invariant() == v] 
    # Holder for new_paths
    new_paths = []

    while steps_left > 0:
      # For each path in current_paths, add on all non-backtracking options
      for path in current_paths:
        # Find the last iso in path
        last_iso = path[-1]
        # Over all isomorphisms, find the non-backtracking ones

        for iso in isos:
          # Only look at isos out of the end of previous isogeny
          if iso.domain() == last_iso.codomain():
            # Add iso only if it is not backtracking:
            # This case is when you don't return to the same vertex
            if iso.codomain() != last_iso.domain():
              new_paths.append(flatten([path, iso]))
            # This case is when you return to the same vertex, but you are not backtracking
            elif all([iso.post_compose(aut) != last_iso.dual() for aut in iso.codomain().automorphisms()]):
              new_paths.append(flatten([path, iso]))

      # Reset <current_paths>
      current_paths = new_paths
      new_paths = []

      # Subtract one from steps_left
      steps_left = steps_left - 1

    # When while loop is done, add all the paths of given length from v to nb_paths
    nb_paths_of_length_n = nb_paths_of_length_n + current_paths

  #Find the cycles among the paths
  cycles = [path for path in nb_paths_of_length_n if path[0].domain() == path[-1].codomain()]

  #Finally, check for tails
  tailless_cycles = [cycle for cycle in cycles if all([cycle[0].post_compose(aut) != cycle[-1].dual() for aut in cycle[0].codomain().automorphisms()])]

  return tailless_cycles


### Function to return the matrix of the non-backtracking operator on edges of an isogeny graph
##  Inputs: G - Isogeny graph
##  Outputs: M - matrix of the action of the non-backtracking operator
def get_non_backtracking_operator(G):
    # Get all isogenies in G
    isos = get_isogenies(G)

    # Make matrix
    M = Matrix(ZZ, len(isos), lambda i,j : 1 if (isos[j].domain() == isos[i].codomain() and not any([isos[i].dual() == isos[j].post_compose(aut) for aut in isos[j].codomain().automorphisms()])) else 0)

    return M


### Function to return the generating function u (d/du)log(zeta) for the number of non-backtracking, tailess paths of length n
### directly from G
##  Inputs: G - Isogeny graph
##  Outputs: gen_fun
def get_cycle_gen_fun(G):

  R.<t> = ZZ[]

  ## Compute powers for the zeta function, exactly as we do in the code to compute the zeta function
  # Get adjacency matrix and isogeny degree from G
  A = G.adjacency_matrix()
  ell = G.isogeny_degree()
  p = G.prime()

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

  # Use formula, derivation is currently on paper, but hopefully will be in the overleaf soon (if this works)
  gen_fun = ((-2*num_exp + n_sd)*(t^2/(1 - t^2)) - n_sd*(t/(1 - t^2)))- t*(((1 - A*t + ell*t^2).inverse())*(-A + 2*ell*t)).trace()
  var('u')
  gen_fun = gen_fun(t = u)

  return gen_fun


### Function to return the number of isogeny cycles using the formulas in Orientations and Cycles in Supersingular Isogeny Graphs
##  Inputs: G - Isogeny graph, r - integer giving the length of the cycle
##  Outputs: n - number of isogeny cycles
def count_cycles_by_orientations_formula(G, r):

  p = G.prime()
  ell = G.isogeny_degree()

  # This is just because I only implemented the simplest case of their formulas
  assert(ell^r < p)

  # Get the discriminants and class numbers giving cycles of length <r>
  # Note that this uses code from my <Distribution_of_cycles.sage> file
  class_ns = [tup[1] for tup in get_discriminants_by_ell_order_fast(ell, r, p = p)]

  # Implement formula
  return (2/r)*sum(class_ns)