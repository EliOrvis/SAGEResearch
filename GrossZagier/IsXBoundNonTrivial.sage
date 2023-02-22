### Find explicit examples where the bound |x| < m \sqrt{D} gives more information about the possibilities for D mod p than just that D must be a square

## Inputs are p, m, where p is the prime and m is the prime power walk of interest

def Find_possible_D(p, m, small_bound = True, compare_BonehLove = False):
	# Only look at discriminants at most (p-1)^2/12m^2:
		# This bound is chosen because 3 is the smallest fundamental imaginary quadratic discriminant.
		# If d_1 is greater than this bound then d_1 d_2 will be guaranteed to be larger than (p-1)^2 / 4m^2,
		# and will therefore be thrown out in the final step.
	if small_bound == False:
		discriminant_bound = floor((p-1)^2)*((12*m^2)^(-1)) # This step appears to make a HUGE difference in the runtime, using an arbitrary smaller bound gives much faster performance.
	#This default below is not optimal but is better for testing:
	else:
		discriminant_bound = floor(p/3)

	#Make Z/pZ to check whether p is inert, i.e. whether O_d is in B_p,infinity
	R = IntegerModRing(p)

	#Make list of fundamental discriminants
	fun_discriminants = [n for n in (1..discriminant_bound) if is_fundamental_discriminant(-n) == True and R(-n).is_square() == False]

	#Products of discriminants, including d_1, d_2 for reference
	products = [[d_1, d_2, d_1*d_2] for d_1 in fun_discriminants for d_2 in fun_discriminants if d_1 <= d_2 and gcd(d_1,d_2) == 1]

	#Make list of possible D values within bound, including d_1, d_2 for reference
	#This bound guarantees that the x_bound in the following function will be small enough that not every square mod p will appear as m^{-1}x^2 for some x within the x_bound.
	#See Research Journal February 12, 2023 for more information.
	#The second condition ensures that p can divide a number of the form (mD - x^2) / 4, i.e. p needs to be smaller than the largest such number, which is mD - x^2.
	#  Unclear if we need this second condition to guarantee only true negatives, since if none of the squares line up with D mod p, then p can't divide any such number (regardless of the size of D)
	#  This bound might be important for positives though, if we have an interest in those in the future.
	possible_D = [l for l in products if l[2] < ((p-1)^2)*(4*m^2)^(-1) and p < (m*l[2]/4)]

	#If comparing against Boneh-Love paper, remove any pairs of discriminants where Boneh-Love already tells us there is no path of length m between d_1 and d_2
	if compare_BonehLove == True:
		possible_D = [l for l in possible_D if ceil(sqrt(p)/(2*determine_minimum_norm(-l[0],-l[1]))) <= l[2]]

	return possible_D


### Given d_1, d_2, p, m, check whether there can be a path of length m from a curve with endomorphisms by d_1 to a curve with endomorphisms by d_2
### BE VERY CAREFUL: True means that there might be a path of the length in question, not that there is a path
## Inputs are D_tuple = [d_1, d_2, d_1*d_2], p, m
def exists_path_possible(D_tuple, p, m):
	# Find bound for |x| < m\sqrt{D}
	x_bound = floor(m*sqrt(D_tuple[2]))

	#Find squares for x within range
	squares_in_range = [x^2 for x in (1..x_bound)]

	#Make Z/pZ to check congruence conditions for D
	R = IntegerModRing(p)

	#Adjust squares by m^{-2} mod p
	squares_in_range_adj = [R(1/(m^2))*R(t) for t in squares_in_range]

	#Return true if D is in adjusted squares and false otherwise
	return R(D_tuple[2]) in squares_in_range_adj


### Helper function that determines the smallest upper bound on the minimum norms of the non-integral elements in the two quadratic rings
def determine_minimum_norm(d_1,d_2):
	if d_1 % 4 == 0:
		d1_bound = d_1/(-4)
	else:
		d1_bound = (1 - d_1)/4
	if d_2 % 4 == 0:
		d2_bound = d_2/(-4)
	else:
		d2_bound = (1 - d_2)/4
	return max(d1_bound, d2_bound)
