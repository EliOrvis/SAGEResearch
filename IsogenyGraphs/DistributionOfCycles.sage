#### This file contains functions used to collect data on the distribution of cycles in the ell-isogeny graph
## WARNING: Several of my internal files are required for this function to work. It is assumed that the running
##          copy of sage has these loaded or attached.


# We will use pandas to construct dataframes and work with data, and use pickle to save data
import pandas as pandas
import pickle as pkl
import numpy as np
import itertools

## This function creates data on the odd cycles of a given length for a range of primes and ell values
#  Inputs:  prime_range - range of primes (p in isogeny notation) to create data for, given as list of ints;
#           isogeny_range - range of isogeny degrees (ell in iso notation), given as list of ints;
#           cycle_length - length of cycles to consider (as number of edges), given as int;
#           outfile_name - string path to where output file should be saved;
#           one_mod_twelve - bool flag, if true, we only consider primes that are 1 mod 12;
#  Outpus:  None - data is stored at outfile_name

def create_cycle_data(prime_range, isogeny_range, cycle_length, outfile_name, one_mod_twelve = False):  

    # Create column headers for cycles. These have to be dynamic to account for different cycle lengths
    cycle_name = 'n_' + str(cycle_length) + '_cycles'
    spine_cycle_name = 'n_' + str(cycle_length) + '_cycles_spine'
    ratio_name = 'ratio_length_' + str(cycle_length)

  # Create a dictionary whose keys will be column headers in the output
    cycle_data = {'p': [], 'ell': [], 'n_vertices' : [], 'n_spine_vertices' : [], cycle_name : [], spine_cycle_name : [],
                  ratio_name : []}
    
    # If one_mod_twelve == True, restrict to primes that are 1 mod 12
    if one_mod_twelve:
        prime_range = [p for p in prime_range if p % 12 == 1]

    # Loop over primes in prime_range, and add data to dictionary
    for ell in isogeny_range:
        for p in prime_range:
        
            # Create isogeny graph
            # NOTE: When Isogeny Graph constructor is fixed to have the option of making the actual directed graph, change this line to do that
            G = IsogenyGraph(prime = p, isogeny_degree = ell).to_directed()
        
            # Compute vertices in spine
            verts = G.vertices()
            spine_verts = [v for v in verts if v^p == v]

            # Find all cycles of length equal to cycle length (sage counts cycles by number of vertices, so we add one to number of edges)
            # NOTE: This is non-optimized, Sage has this function to return all simple cycles up to some length
            #       but not exactly some length.
            #       Also, this counts direction of traversal as well. This doesn't affect the ratio but is good to keep in mind.
            cycles = [c for c in G.all_simple_cycles(max_length = cycle_length + 1) if len(c) == cycle_length + 1]

            # Find cycles intersecting the spine
            spine_cycles = [c for c in cycles if any([vert in spine_verts for vert in c])]

            # Add data to dictionary
            cycle_data['p'].append(p)
            cycle_data['ell'].append(ell)
            cycle_data['n_vertices'].append(len(verts))
            cycle_data['n_spine_vertices'].append(len(spine_verts))
            cycle_data[cycle_name].append(len(cycles))
            cycle_data[spine_cycle_name].append(len(spine_cycles))
            if len(cycles) > 0:
              cycle_data[ratio_name].append(len(spine_cycles)/len(cycles))
            else:
              cycle_data[ratio_name].append(np.nan)
             
    # Save results as a pandas dataframe
    df = pd.DataFrame(data = cycle_data).set_index('p')
    
    # Write results as .pkl file for future use.
    df.to_pickle(outfile_name + ".pkl")
    
    return

  ## This function returns all imaginary quadratic discriminants where <ell> has order <order>
  #  Inputs: ell - prime; order - positive integer
  #  Outputs: discs - list of tuples of the form (d, N), where d is an imaginary quadratic discriminant
  #                   (not necessarily fundamental) such that 
  #                   (1) <ell> does not divide the conductor of O_d
  #                   (2) The order of <ell> in the class group of O_d is <order>
  #                   (3) <ell> splits in O_d,
  #                   and N is the class number of O_d

def get_discriminants_by_ell_order(ell, order):

  power = ell^order

  # Compute discriminant bound - this uses that the smallest number other than 1
  #   represented by x^2 + p_d xy + (p_d - d)/4 y^2 is (p_d - d)/4.
  # This is the bound such that any D representing <power> has -D <= <bound>
  disc_bound = 4*(power) - 1

  # We will store relevant discriminants here
  discs = []

  # Loop through discriminants:
  for d in [1..disc_bound]:

    # skip if -d is not 0 or 1 mod 4
    if -d % 4 in [2,3]:
      continue

    # Set binary quadartic form
    pd = ZZ(d % 2) 
    Q = QuadraticForm(ZZ, 2, [1, pd, (pd + d)/4])

    # Skip if order wtih discriminant d is not ell fundamental
    conductor = (-d / fundamental_discriminant(-d))^(1/2)

    if ell.divides(conductor):
      continue

    # Since the order is ell-fundamental, we can use Dedekind-Kummer to check if ell splits.
    #  Otherwise skip.
    R.<x> = ZZ[]
    poly = R(Q.polynomial()(x, 1))

    ell_factors = poly.change_ring(GF(ell)).factor()
    
    if len(ell_factors) == 1:
      continue

    # Skip if binary quadratic form represents a power of <ell> less than <ell^order>.
    # These are the cases where the order of ell is less than order.

    # These are vectors of all representations for values at most ell^i
    representations = Q.representation_vector_list(power + 1)

    # These are the representations of powers of ell less than ell^order
    ell_reps = flatten([representations[ell^i] for i in [1..order - 1]])

    # check if any of the representations are proper (no zeros and coprime entries), and skip if so:
    if any([(0 not in vect) and gcd(vect) == 1 for vect in ell_reps]):
      continue
    

    # Return (-d, N) if all else has skipped, and Q represents ell^order properly:
    power_reps = representations[power]
    
    if any([(0 not in vect) and gcd(vect) == 1 for vect in power_reps]):
      # This is the case we want, so we compute the class number of the order of discriminant -d
      # <poly> has discriminant -d, so we generate the order defined by poly
      O.<a> = EquationOrder(poly)
      
      # Validation to make sure I didn't mess up any of the cases:
      if O.discriminant() != -d:
        raise RuntimeError("There is an issue with the implementation: the order for the class number does not have the right discriminant!")

      discs.append((-d, O.class_number()))

  return discs 

## This function takes in a an <ell> value, and a list of cycle lengths, and returns moduli conditions on p that will guarantee no cycles of those lengths
## WARNING: THIS TAKES A REALLY LONG TIME. IF YOU JUST WANT TO CHECK WHETHER THERE EXIST SOLUTIONS, SET <check_n> TO BE THE NUMBER OF SOLUTIONS YOU WANT
#  Inputs: ell - positive prime representing isogeny degree; cycle_list - list of positive integers >= 3 representing cycle lengths
#          check_n - Boolean flag, if True then the program will just return the modulus and solution as soon as it finds one solution, rather than looking for all of them
#  Outputs: res_dict - dictionary whose key is a modulus with values the values that guarantee if p is equivalent to one of these values mod the modulus, 
#                       then there are no cycles of length any of the lenghts in <cycle_list>.
def p_conditions_to_prohibit_cycles(ell, cycle_list, check_n = 20):

  # Get discriminant list.
  discs = []
  for length in cycle_list:
    disc_info = get_discriminants_by_ell_order(ell, length)
    just_discs = [d[0] for d in disc_info]
    discs = discs + just_discs

  # Pass this to p_split_residues_from_discs
  res_dict = p_split_residues_from_discs(discs, check_n = check_n)

  return res_dict


## This function returns the list of residues where all the discriminants in <disc_list> will have p split 
## Note that this uses the much longer (but conceptually simple) function p_split_residues (below)
#  Inputs: disc_lst - list of imaginary quadratic discriminants, just_check - Boolean, if True the function will
#                     just return the modulus and one solution once a solution is found.
#  Outputs: res_dict - dictionary whose keys are moduli and values are residues in that modulus that will guarantee
#                      that p splits in all the discriminants
def p_split_residues_from_discs(disc_list, check_n):

  # Create a dictionary whose keys are the modulus condition for each discriminant, and whose values are the moduli making p split
  disc_dict = {}
  for d in disc_list:
    d_dict = p_split_residues(d)
    for key in d_dict:
      # The line below is in case the same modulus shows up for two different orders
      # This can happen, but I think we'll always end up with the same conditions anyway, so
      #   this is probably unnecessary.
      if key in disc_dict.keys():
        disc_dict[key] += d_dict[key]
        # This line is just removing duplicates
        disc_dict[key] = list(set(disc_dict[key]))
      else:
        disc_dict[key] = d_dict[key]

  # New modulus is the lcm of the moduli for the inputs
  key_list = [key for key in disc_dict]
  modulus = lcm(key_list)

  res_dict = {modulus : []}

  # For each combination of one value from each modulus, we need to check if that combo has a CRT solution, and save it if so
  values_list = [disc_dict[key] for key in disc_dict]

  for value_set in itertools.product(*values_list):
    try:
      value = CRT_list(list(value_set), key_list)
      res_dict[modulus].append(value)
      if len(res_dict[modulus]) >= check_n:
        return res_dict
    except:
      continue

  return res_dict




## This function returns the list of residues mod the squarefree part of d or 4 or 8 times this where d is a square mod p
#  Inputs: d - imaginary quadratic discriminant
#  Outputs: res_dict - dictionary <{modulus: [list of residues]}>


def p_split_residues(d):
  
  # Get squarefree part of -d 
  sqrf_d = (-d).squarefree_part()

  # All squarefree factors odd factors. We remove 2 if it occurs since we handle this separately
  #  in the quadratic reciprocity arguments
  sqrf_factors = [d[0] for d in sqrf_d.factor() if d != 2]

  # These will store lists of values mod factors where legendre symbol is 1 (or 0), -1, respectively
  # NOTE: we exclude legendre symbol 0 in the first case since we don't want p to be equal to one of the factors
  #       (this would mean p is ramified in O_d, so O_d appears in the SS-ellI graph) 

  fact_quad_res = {fact: [x for x in range(fact) if kronecker(x, fact) == 1] for fact in sqrf_factors}
  fact_quad_nonres = {fact: [x for x in range(fact) if kronecker(x, fact) == -1] for fact in sqrf_factors}

  # To determine congruence conditions, we will check how many primes dividing <sqrf_d> are 3 mod 4
  three_mod_four_factors = [fact for fact in sqrf_factors if fact % 4 == 3]

  # If this is not divisible by 2, we get congruence conditions modulo either <sqrf_d> or <4*sqrf_d>
  if 2.divides(sqrf_d) == False:

    # If there are an odd number of factors that are 3 mod 4, then the congruence conditions
    # determining whether d is a square mod p are just mod <sqrf_d>, otherwise mod <4*sqrf_d>
    #   See RJ 09/24/2023 (or purple notebook at home if I don't get around to writing this up)  
    if len(three_mod_four_factors) % 2 == 1:

      # This dictionary will store outputs
      res_dict = {sqrf_d : []}

      # Loop over all subsets of even size of factors, and choose any quadratic non-residues from
      # those factors, and quadratic residues from all other factors. Then find crt congruences 
      # from these lists.

      # This iterator will let us loop over 0, 2, 4, etc.
      for i in range(floor(len(sqrf_factors)/2) + 1):

        # This iterator loops over all choices of a particular set of size 2*i of factors
        for comb in itertools.combinations(sqrf_factors, 2*i):
          
          # Choose the quadratic non-residues for the factors in comb and residues for the others
          residues = [fact_quad_res[factor] if factor not in comb else fact_quad_nonres[factor] for factor in sqrf_factors]
          
          # For every combination of residues, add the modulus given by the chinese remainder theorem
          for ress in itertools.product(*residues):
            res_dict[sqrf_d].append(CRT_list(list(ress), sqrf_factors))

    # This case has a similar structure, but now the congruences are modulo 4*sqrf_d, because the cases are different depending
    # on p mod 4.
    else:

      res_dict = {4*sqrf_d : []}

      # If p is 1 mod 4, we do the same process as the last case, but include the condition that p is 1 mod 4 in the final CRT call
      for i in range(floor(len(sqrf_factors)/2) + 1):
        for comb in itertools.combinations(sqrf_factors, 2*i):
          residues = [fact_quad_res[factor] if factor not in comb else fact_quad_nonres[factor] for factor in sqrf_factors]
          for ress in itertools.product(*residues):
            # Include 1 in the values and 4 in the moduli in this call to make sure we are taking p % 4 == 1
            res_list = list(ress)
            res_list.append(1)
            sqrf_list = [fact for fact in sqrf_factors]
            sqrf_list.append(4)
            res_dict[4*sqrf_d].append(CRT_list(res_list, sqrf_list))

      # The other case is p % 4 == 3, in this case we want to choose an odd number of -1s instead of even
        for comb in itertools.combinations(sqrf_factors, 2*i + 1):
          residues = [fact_quad_res[factor] if factor not in comb else fact_quad_nonres[factor] for factor in sqrf_factors]
          for ress in itertools.product(*residues):
            # Include 3 in the values and 4 in the moduli in this call to make sure we are taking p % 4 == 3
            res_list = list(ress)
            res_list.append(3)
            sqrf_list = [fact for fact in sqrf_factors]
            sqrf_list.append(4)
            res_dict[4*sqrf_d].append(CRT_list(res_list, sqrf_list))

  # If <sqrf_d> is divisible by 2, we get congruence conditions modulo <8*sqrf_d>
  if 2.divides(sqrf_d) == True:
    
    res_dict = {8*(sqrf_d/2) : []}

    # If there are an odd number of square free factors that are 3 mod 4, then we want a +1 in the product when p is 1, -1 mod 8 and -1 otherwise
    #   See RJ 09/24/2023 (or purple notebook at home if I don't get around to writing this up)  
    if len(three_mod_four_factors) % 2 == 1:

      # These are the 1, -1 mod 8 cases, we want an even number of -1s here
      for i in range(floor(len(sqrf_factors)/2) + 1):
        for comb in itertools.combinations(sqrf_factors, 2*i):
          residues = [fact_quad_res[factor] if factor not in comb else fact_quad_nonres[factor] for factor in sqrf_factors]
          for ress in itertools.product(*residues):
            # Include 1, 7 in the values and 8 in the moduli in this call to make sure we are taking p % 8 == pm 1
            res_list_1 = list(ress)
            res_list_7 = list(ress)
            res_list_1.append(1)
            res_list_7.append(7)
            sqrf_list = [fact for fact in sqrf_factors]
            sqrf_list.append(8)
            res_dict[8*(sqrf_d/2)].append(CRT_list(res_list_1, sqrf_list))
            res_dict[8*(sqrf_d/2)].append(CRT_list(res_list_7, sqrf_list))

      # These are the 3, -3 mod 8 cases, we want an odd number of -1s here
        for comb in itertools.combinations(sqrf_factors, 2*i + 1):
          residues = [fact_quad_res[factor] if factor not in comb else fact_quad_nonres[factor] for factor in sqrf_factors]
          for ress in itertools.product(*residues):
            # Include 3, 5 in the values and 8 in the moduli in this call to make sure we are taking p % 8 == pm 1
            res_list_3 = list(ress)
            res_list_5 = list(ress)
            res_list_3.append(3)
            res_list_5.append(5)
            sqrf_list = [fact for fact in sqrf_factors]
            sqrf_list.append(8)
            res_dict[8*(sqrf_d/2)].append(CRT_list(res_list_3, sqrf_list))
            res_dict[8*(sqrf_d/2)].append(CRT_list(res_list_5, sqrf_list))

    # If there are an even number of square free factors that are 3 mod 4, then we want a +1 in the product when p is 1, 3 mod 8 and -1 otherwise
    #   See RJ 09/24/2023 (or purple notebook at home if I don't get around to writing this up)  
    else:

       # These are the 1, 3 mod 8 cases, we want an even number of -1s here
      for i in range(floor(len(sqrf_factors)/2) + 1):
        for comb in itertools.combinations(sqrf_factors, 2*i):
          residues = [fact_quad_res[factor] if factor not in comb else fact_quad_nonres[factor] for factor in sqrf_factors]
          for ress in itertools.product(*residues):
            # Include 1, 3 in the values and 8 in the moduli in this call to make sure we are taking p % 8 == pm 1
            res_list_1 = list(ress)
            res_list_3 = list(ress)
            res_list_1.append(1)
            res_list_3.append(3)
            sqrf_list = [fact for fact in sqrf_factors]
            sqrf_list.append(8)
            res_dict[8*(sqrf_d/2)].append(CRT_list(res_list_1, sqrf_list))
            res_dict[8*(sqrf_d/2)].append(CRT_list(res_list_3, sqrf_list))

      # These are the 5, 7 mod 8 cases, we want an odd number of -1s here
        for comb in itertools.combinations(sqrf_factors, 2*i + 1):
          residues = [fact_quad_res[factor] if factor not in comb else fact_quad_nonres[factor] for factor in sqrf_factors]
          for ress in itertools.product(*residues):
            # Include 5, 7 in the values and 8 in the moduli in this call to make sure we are taking p % 8 == pm 1
            res_list_5 = list(ress)
            res_list_7 = list(ress)
            res_list_5.append(5)
            res_list_7.append(7)
            sqrf_list = [fact for fact in sqrf_factors]
            sqrf_list.append(8)
            res_dict[8*(sqrf_d/2)].append(CRT_list(res_list_5, sqrf_list))
            res_dict[8*(sqrf_d/2)].append(CRT_list(res_list_7, sqrf_list))

  return res_dict

## This function takes in a discriminant, prime, and returns the number of cycles along the spine, as long as the prime is larger than the discriminant 
#  Inputs: d - an imaginary quadratic discriminant; p - a prime number larger than the absolute value of d
#  Outputs: n - the number of cycles that the discriminant produces along the spine.
#           Note that the result is the number of *undirected* cycles, and therefore half of the number returned by G.all_simple_cycles()
def n_spine_cycles(d, p):

  if abs(d) >= p:
    raise ValueError("For accuracy, absolute value of d must be less than or equal to p.")

  # If p is not inert in the order of discriminant d, then d gives no cycles (along the spine or otherwise).
  if legendre_symbol(d,p) == 1:
    return 0

  # Get prime factors of d
  prime_factors = d.prime_factors()

  # d gives no cycles if the coditions of Chen, Xue are not satisfied:
  for div in prime_factors:
    if div == 2 and not (div % 8 == 7 or (-p + (d/4)) % 8 in [0,1,4] or (-p + d) % 8 == 1):
      return 0
    if div != 2 and legendre_symbol(-p,div) != 1:
      return 0

  # If none of these conditions hold, then return the genus number of d
  return genus_number_of_d(d)

## This function takes in a discriminant, class number, prime, ell, and cycle length, and returns the number of cycles produced by that discriminant 
#  Inputs: d - an imaginary quadratic discriminant; cln - class number of d; p - a prime number larger than the absolute value of d; length - cycle length
#  Outputs: n - the number of cycles that the discriminant produces in G_p.
#           Note that the result is the number of *undirected* cycles, and therefore half of the number returned by G.all_simple_cycles()
def n_cycles_from_d(d, cln, p, ell, length = False):

  if abs(d) >= p:
    raise ValueError("For accuracy, absolute value of d must be less than or equal to p.")
  if not length.divides(cln):
    raise ValueError("This function should only be run for discriminants where the length of the cycle divides the class number.")

  # If p is not inert in the order of discriminant d, then d gives no cycles (along the spine or otherwise).
  if legendre_symbol(d,p) == 1:
    return 0

  # If length != False, the user gave a specific length, so we validate that ell has that order in the class group of discriminant d
  K = QuadraticField(d)
  frakells = K.primes_above(ell)
  # validate the ell splits in K
  if len(frakells) != 2:
    raise ValueError("ell must split in the quadratic field with discriminant d.")

  # If everything passes, then we return cln / length.
  return cln / length

## This function takes as input a range of primes, an ell value, and an r value, and returns a dataset of counts of r cycles in the ell isogeny graph, with number along the spine
## Note that this differs from create_cycle_data, in that here we use the results of my cycles paper, so we can consider MUCH larger primes.
#  Inputs: p_range - <[lower, upper]>, where <lower> is lower bound on primes to consider, and <upper> is upper bound; ell - isogeny degree; r - cycle length
#  Outputs: df - pandas dataframe with columns p, ell, r, nt (total # of cycles), ns (total # along the spine), total_verts (# of vertices in graph), spine_verts (# of spine vertices),
#           and freq_ratio ((ns/spine_verts)/(nt/total_verts)).
def create_cycle_data_from_class_numbers(p_range, ell, r):

  # First, verify that I have the requested set of discriminants:
  # NOTE: This assumes that the discriminant data is stored in the right place in my file system
  try:
    with open("/mnt/c/Users/welio/Documents/Boulder/Research/SAGE/Data" + "/IsogenyGraphs/cycle_discriminants_" + str(ell) + "_" + str(r) + ".pkl", "rb") as f:
      # Load discriminants - we assume discriminants are stored as a list of tuples (disc, class number)
      discriminant_data = pkl.load(f)
  except:
    raise ValueError("No discriminant data for ell = %s and r = %s"%(ell, r))

  # Validate that user gave inputs that are large enough for my theorems to hold
  just_discriminants = [d for (d, cln) in discriminant_data]
  assert p_range[0] > max(just_discriminants)*(sorted(just_discriminants, reverse = True)[1])/4

  # Create dictionary to store data in
  dict_data = {"p" : [], "ell" : [], "r" : [], "nt" : [], "ns" : [], "total_verts" : [], "spine_verts" : [], "freq_ratio" : []}

  # Just to keep me sane.
  counter = 0
  print("Goal: %s"%(p_range[1]))
  for p in prime_range(p_range[0], p_range[1]):
    counter += 1
    if counter == 50:
      print(p)
      counter = 0

    # Total number of cycles is computed by adding up counts for each discriminant 
    nt = sum([n_cycles_from_d(d,cln,p,ell, length = r) for (d, cln) in discriminant_data])
    # Number of spine cycles is computed by adding up counts for each discriminant
    ns = sum([n_spine_cycles(d,p) for (d,cln) in discriminant_data])

    # Total vertices and spine vertices are computed directly
    total_verts = n_ss_j_invars(p)
    spine_verts = n_Fp_j_invars(p)

    # Ratio is computed as a float:
    freq_ratio = float((ns/spine_verts)/(nt/total_verts))

    # Data is stored for this prime:
    dict_data["p"].append(p)
    dict_data["ell"].append(ell)
    dict_data["r"].append(r)
    dict_data["nt"].append(nt)
    dict_data["ns"].append(ns)
    dict_data["total_verts"].append(total_verts)
    dict_data["spine_verts"].append(spine_verts)
    dict_data["freq_ratio"].append(freq_ratio)

  return pd.DataFrame(data = dict_data)
