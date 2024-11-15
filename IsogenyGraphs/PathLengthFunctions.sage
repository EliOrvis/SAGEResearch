##### This file contains methods for extracting path lengths between different
##### sets of vertices in isogenygraphs

import pandas as pd

pd.set_option('display.max_columns', None)

### Obtains the set of shortest path lengths between valleys with given fundamental discriminants
## Inputs: G - IsogenyGraph object, d1 & d2 - fundamental discriminants, table_format - if true, return pd Dataframe with nice indices for viewing, otherwise dictionary 
## Output: dictionary with keys (j1, j2) where j1 is a j-invariant of the SMALLER (in abs value) discriminant 
##         and j2 is a j-invariant of the LARGER (in abs value) discriminant, and values the shortest path length between them
def path_lengths_between_valleys(G, d1, d2, table_format = False):
	if abs(d1) <= abs(d2):
		d_small = d1
		d_large = d2
	else:
		d_small = d2
		d_large = d1

	vertexset1 = get_CM_vertices(G, d_small); vertexset2 = get_CM_vertices(G, d_large)

	lengths = {}
	
	#Iterate over vertices in vertexset1 and vertexset2 to get path lengths
	for w in vertexset2:
		# Put all lengths from valley 1 into valley 2 into a list
		w_lengths = []

		for v in vertexset1:
			w_lengths.append(G.shortest_path_length(v,w))

		# Add lengths to w under the key w
		lengths[w] = w_lengths

	if table_format == True:
		# Convert into pandas dataframe for easier viewing
		df = pd.DataFrame(data = lengths)	
		# Set indices to be vertices from valley 1
		s = pd.Series(vertexset1)

		return df.set_index(s)

	else:
		return lengths


### Counts the number of paths of given length between two given vertices; 
###    the only paths counted are those that do not backtrack along edges, i.e. if we take 
###    an edge e from u to v, the subsequence edge cannot be the same edge back from v to u, 
###    but it can be a different edge from v to u. This corresponds to counting cyclic isogenies,
###    except maybe at the vertices with extra automorphisms, where there could be additional cyclic isogenies.
## Inputs: G - a graph object; length - length of path to search for
## Outputs: Matrix of path lengths, with indices given by the vertices of G in sorted order

# Note: This function does not actually find the paths (at all). Instead, it uses combinatorics to quickly count such paths.
# Note 2: On 1/25/2024, long after writing this, I noticed that the matrix sometimes has negative entries. I don't know what these mean! 
#         They could even be messing up the count of connections between valleys that is done later on, but I thought I had tested this?


def paths_of_length_n_matrix(G, length):

	# Get adjusted (i.e. loops counted as 2) adjacency matrix for G 
	A = G.adjacency_matrix()
	# Make diagonal matrix with vertex degree on diagonal
	# Note that sort = True ensures this is in the same order as in making the AM
	verts = G.vertices(sort=True)
	D = diagonal_matrix([G.degree(v) for v in verts])
	# Make identity matrix
	I = identity_matrix(len(verts))

	#create variables Atemp,Dtemp,z
	Atemp,Dtemp,z = var('Atemp,Dtemp,z')
	# This is the generating function for the # NBW matrix when we set Itemp = I, Atemp = A, and Dtemp = D 
	# See RJ 5/5/2023 - this is from Kempton 2016
	genfun = -(z^2 - 1)/((Dtemp - 1)*z^2 - Atemp*z + 1)

	breakpoint()
	# Get coefficient on x^length in series expansion
	# Evaluating this will give the matrix of NBW counts
	coeff = genfun.series(z, length + 1).coefficient(z,length)	

	# Create function field in 4 variables, this is necessary because sage won't evaluate symbolic expressions on matrices
	R.<Atemp,Dtemp,z> = PolynomialRing(QQ,3)
	F = R.fraction_field()
	coeff_func = F(coeff)

	# Get matrix of counts
	M = coeff_func(Atemp = A, Dtemp = D)

	#return requested entry
	return M


### Returns the pairs of vertices between two valleys separated by a path of a given length
## Inputs: G - an isogeny graph object; d1, d2 - embedded fundamental discriminants (these can be the same); length - length to search for
## Outputs: list of pairs of vertices with optimally embedded copy of maximal order in d1, d2 respectively, separated by a path of requested length

def CM_pairs_separated_by_path_length(G, d1, d2, length):
	vertexset1 = get_CM_vertices(G, d1); vertexset2 = get_CM_vertices(G, d2)

	pairs = []

	# Get sorted vertices so that indices match the matrix indices
	verts = G.vertices(sort = True)

	M = paths_of_length_n_matrix(G,length)

	for v in vertexset1:
		for w in vertexset2:
			if M[verts.index(v)][verts.index(w)] > 0:
				pairs.append((v,w))

	return pairs

### This is a helper function to iterate over all edges out of a vertex in a particular (fixed) order
##	Inputs: G - an isogeny graph; v - a vertex of G
##	Outputs: None - yields (w, i) for i between 0 and number of edges from v to w minus 1 (no output w/ w if there are no such edges)
def iterate_edges(G, v):
	# Turn G into a dictionary
	G_dict = G.to_dictionary(multiple_edges = True, edge_labels = True)

	for key in G_dict:
		for sub_key in G_dict[key]:
			for i in [0..len(G_dict[key][sub_key]) - 1]:
				yield (key, i)


### Returns the number of (non-backtracking) paths between two vertices of a given length n, counted directly from the graph
## Inputs: G - an isogeny graph object; start, target - vertices of G; n - path lengths to search for, isos - set of isogenies for <G>
##         indexed_paths - Boolean for whether to return explicit paths as list of indices in <isos> or actual isogenies
## Outputs: n - number of paths of specified length
def non_backtracking_paths(G, start, target, n, isos = None, indexed_paths = False):
  # First, get the isogenies if not already given
  if isos == None:
    isos = get_isogenies(G)

  # Verify that the number of isogenies equals the number of edges.
  # This is just a simple check in case someone passes the wrong isogeny set for G
  assert(len(G.edges()) == len(isos))

  # Find non-backtracking paths of length <n> beginning at <start>
  nb_paths_of_length_n = []

  # Counter for how many steps we need to take.
  steps_left = n - 1
  # List of current paths
  current_paths = [[iso] for iso in isos if iso.domain().j_invariant() == start] 
  # Holder for new_paths
  new_paths = []

  if steps_left == 0:
    nb_paths_of_length_n = current_paths

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

  #Find the paths from <start> to <target>
  # If <target> is None, return all non-backtracking paths from <start> of length <n>
  if target == None:
    paths_temp = nb_paths_of_length_n
  else:
    paths_temp = [path for path in nb_paths_of_length_n if path[-1].codomain().j_invariant() == target]

  if indexed_paths:
    #Relabel isogenies in paths with index in isos to be able to compare across paths
    paths = [[isos.index(iso) for iso in path] for path in paths_temp]
    return paths

  else:
    return paths_temp

### Returns a dictionary of the number of paths of a given length between the vertices in two isogeny valleys
## Inputs: G - an isogeny graph object; d1, d2 - embedded discriminants for <G>; n - path lengths to search for, isos - set of isogenies for <G>
##         explicit_paths - boolean indicating whether to return the actual paths, or just the count
## Outputs: out_dict - dictionary of path counts
def n_paths_between_valleys(G, d1, d2, n, isos = None, explicit_paths = False):
  # First, get the isogenies if not already given
  if isos == None:
    isos = get_isogenies(G)

  # Get relevant vertices
  vertexset1 = get_CM_vertices(G, d1); vertexset2 = get_CM_vertices(G, d2)

  out_dict = {}
  for v in vertexset1:
    for w in vertexset2:
      if explicit_paths:
        out_dict[(v,w)] = non_backtracking_paths(G, v, w, n, isos = isos, indexed_paths = True)
      else:
        out_dict[(v,w)] = len(non_backtracking_paths(G,v,w,n,isos= isos))

  return out_dict