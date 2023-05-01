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


### Determines whether there exists a path of given length between two given vertices
## Inputs: G - a graph object; v, w - vertices of G; length - length of path to search for
## Outputs: True/False boolean
# NOTE: currently allows for cycles, so, for example, every vertex is separated from itself by a "path" of length 2

def exists_path_of_length(G, v, w, length, visited = []):

	if length == 0:
		if v == w:
			return True
		else:
			return False


	else:
		neighbors = list(set(G.neighbors(v)) - set(visited))
		while len(neighbors) > 0:
			visited.append(v)
			neigh = neighbors.pop()
			if exists_path_of_length(G,neigh,w, length - 1, visited):
				return True
			visited = []
		return False


### Returns the pairs of vertices between two valleys separated by a path of a given length
## Inputs: G - an isogeny graph object; d1, d2 - embedded fundamental discriminants (these can be the same); length - length to search for
## Outputs: list of pairs of vertices with optimally embedded copy of maximal order in d1, d2 respectively, separated by a path of requested length

def CM_pairs_separated_by_path_length(G, d1, d2, length):
	vertexset1 = get_CM_vertices(G, d1); vertexset2 = get_CM_vertices(G, d2)

	pairs = []

	for v in vertexset1:
		for w in vertexset2:
			if exists_path_of_length(G,v,w,length):
				if (v,w) not in pairs and (w,v) not in pairs:
					pairs.append((v,w))

	return pairs