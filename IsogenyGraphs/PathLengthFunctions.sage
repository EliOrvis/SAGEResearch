##### This file contains methods for extracting path lengths between different
##### sets of vertices in isogenygraphs

import pandas as pd

### Obtains the set of path lengths between valleys with given fundamental discriminants
## Inputs: G - IsogenyGraph object, d1 & d2 - fundamental discriminants 
## Output: dictionary with keys (j1, j2) where j1 is a j-invariant of the SMALLER (in abs value) discriminant 
##         and j2 is a j-invariant of the LARGER (in abs value) discriminant, and values the shortest path length between them
def path_lengths_between_valleys(G, d1, d2):
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
			w_lengths.append(G.graph().shortest_path_length(v,w))

		# Add lengths to w under the key w
		lengths[w] = w_lengths

	# Convert into pandas dataframe for easier viewing
	df = pd.DataFrame(data = lengths)	
	# Set indices to be vertices from valley 1
	s = pd.Series(vertexset1)

	return df.set_index(s)