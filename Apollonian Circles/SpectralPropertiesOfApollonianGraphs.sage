#### Create Apollonian graph modulo p, starting from a given quadruple

def create_Apollonian_graph(p, quadruple, projective=False):
    Fp = GF(p)
    if projective == True:
    	P.<a,b,c,d> = ProjectiveSpace(3, Fp)
    	quadrup = P(quadruple)
    else:
    	P.<a,b,c,d> = AffineSpace(4, Fp)
    	quadrup = P(quadruple)
    S1 = P.hom([a,b,c, 2*a + 2*b + 2*c - d], P)
    S2 = P.hom([-a + 2*b + 2*c + 2*d, b, c, d], P)
    S3 = P.hom([a, 2*a - b + 2*c + 2*d,c,d], P) 
    S4 = P.hom([a, b, 2*a + 2*b - c + 2*d, d], P)
    vertices = [quadrup]
    edges = []
    visited = []
    unvisited = copy(vertices)
    i = 0
    while unvisited != []:
        for v in unvisited:
            vertices.append(S1(v))
            vertices.append(S2(v))
            vertices.append(S3(v))
            edges.append((v,S1(v)))
            edges.append((v,S2(v)))
            edges.append((v,S3(v)))
            visited.append(v)
        unvisited = [x for x in vertices if x not in visited]
    return Graph([vertices, edges], format='vertices_and_edges', loops = True)

 # Adjust adjacency matrix to count loops as degree 2
def get_adjusted_AM(G):
    AM = copy(G.adjacency_matrix())
    for i in range(AM.dimensions()[1]):
        if AM[i][i] > 0:
            AM[i,i] = 2*AM[i,i]
    return AM