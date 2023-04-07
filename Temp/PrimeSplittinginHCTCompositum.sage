### Temp file to write some functions that will be helpful for searching
### for the possibility of inertia between H_1H_2 and H_1.

def construct_compositum(d1, d2):
	K1 = QuadraticField(d1); K2 = QuadraticField(d2)

	HCP1 = K1.hilbert_class_polynomial(); HCP2 = K2.hilbert_class_polynomial()

	H1.<j1> = K1.extension(HCP1)

	H2.<j2> = H1.extension(HCP2)

	return H2