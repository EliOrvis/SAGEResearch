#### This file contains functions used to collect data on the distribution of cycles in the ell-isogeny graph
## WARNING: Several of my internal files are required for this function to work. It is assumed that the running
##          copy of sage has these loaded or attached.


# We will use pandas to construct dataframes and work with data, and use pickle to save data
import pandas as pandas
import pickle as pkl
import numpy as np

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
    poly = Q.polynomial()(x, 1)

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
        print(d)
        raise RuntimeError("There is an issue with the implementation: the order for the class number does not have the right discriminant!")

      discs.append((-d, O.class_number()))

  return discs 