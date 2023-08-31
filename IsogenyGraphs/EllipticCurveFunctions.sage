### Function to obtain the isogeny phi' from a separable isogeny phi, known to be m times phi'. 
### This is based on Proposition 2.6 in McMurdy's, Explicit Representation of the Endomorphism Rings of Supersingular Elliptic Curves (see McMurdy for definitions of terms in comments)
## Inputs: phi - an isogeny; m - an integer to divide by
## Outputs: phi' - an isogeny  
def isogeny_division_by_m(phi, m):
	# Initialize domain and codomain curves
	E1 = phi.domain()
	E2 = phi.codomain()

	# Get F(x) and G(x) 
	F = phi.rational_maps()[0]
	G = phi.rational_maps()[1]

	# Get c_F, P(x), Q(x)
	## First we find cF, by getting leading coefficents of numerator and denominator of F(x)
	top_lead = F.numerator().coefficents()[0]
	bot_lead = F.denominator().coefficents()[0]
	cF = top_lead / bot_lead

	# To find P(x), just multiply the numerator by 1/leading coefficient
	P = F.numerator() / top_lead

	# To find Q(x), divide out either the kernel polynomial, or the kernel polynomial squared
	ker_poly = E1.division_polynomial(m)
	if m == 2:
		Q = (F.denominator()/bot_lead).quo_rem(ker_poly)[0]
	else:
		Q = (F.denominator()/bot_lead).quo_rem(ker_poly^2)[0]

	# get x-function of multiplication-by-m on both sides
	m1 = E1.scalar_multiplication(m).rational_maps()[0]
	m2 = E2.scalar_multiplication(m).rational_maps()[1]

	# construct p(x) and q(x)
	return(P,Q)




