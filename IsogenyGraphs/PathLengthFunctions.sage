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
		df = pd.DataFrame(data = path_lengths)	
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
## Inputs: G - a graph object; v, w - vertices of G; length - length of path to search for
## Outputs: integer count of paths

# Note: This function does not actually find the paths (at all). Instead, it uses combinatorics to quickly count such paths.

def n_paths_of_length(G, v, w, length):

	# Get adjusted (i.e. loops counted as 2) adjacency matrix for G 
	A = G.adjusted_AM()
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

	# Get coefficient on x^length in series expanssion
	# Evaluating this will give the matrix of NBW counts
	coeff = genfun.series(z, length + 1).coefficient(z,length)	

	# Create function field in 4 variables, this is necessary because sage won't evaluate symbolic expressions on matrices
	R.<Atemp,Dtemp,z> = PolynomialRing(QQ,3)
	F = R.fraction_field()
	coeff_func = F(coeff)

	# Get matrix of counts
	M = coeff_func(Atemp = A, Dtemp = D)

	#return requested entry
	return M[verts.index(v)][verts.index(w)]


### Returns the pairs of vertices between two valleys separated by a path of a given length
## Inputs: G - an isogeny graph object; d1, d2 - embedded fundamental discriminants (these can be the same); length - length to search for
## Outputs: list of pairs of vertices with optimally embedded copy of maximal order in d1, d2 respectively, separated by a path of requested length

def CM_pairs_separated_by_path_length(G, d1, d2, length):
	vertexset1 = get_CM_vertices(G, d1); vertexset2 = get_CM_vertices(G, d2)

	pairs = []

	for v in vertexset1:
		for w in vertexset2:
			if n_paths_of_length(G,v,w,length) > 0:
				pairs.append((v,w))

	return pairs