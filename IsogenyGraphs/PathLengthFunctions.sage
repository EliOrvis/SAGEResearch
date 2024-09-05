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
## Inputs: G - an isogeny graph object; start, target - vertices of G; n - path lengths to search for
## Outputs: n - number of paths of specified length
#  NOTE: This doesn't work correctly right now. In particular, to get the non-backtracking, we first need to understand the action of duals
#        on the directed graph.
def n_paths_of_length_n(G, start, target, n, seen = None):
	if seen is None:
		seen = {}

	if n == 0:
		if start == target:
			return 1
		else:
			return 0
	else:
		answer = 0
		for v, i in iterate_edges(G, start):
			if (v, start, i) not in seen:
				if (start, v, i) in seen:
					seen.add((start, v, i))
					answer += n_paths_of_length_n(G, v, target, n-1, seen)
					seen.pop((start, v, i))
				else:
					answer += n_paths_of_length_n(G, v, target, n-1, seen)
		return answer