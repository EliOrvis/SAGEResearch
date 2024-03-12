#### Functions for investigating Lehmer's conjecture for elliptic curves

##   Import relevant modules 
import itertools
import pandas as pd
import pickle as pkl

### Function to compute heights of points on an EC, for all x in the maximal order with coefficients up to <xbound>
##  Input: E - elliptic curve; xbound - bound on x (in abs value); outfile_name - path to attach to <data_path> for saving data 
#   Note: This currently searches for only integer values of x. This is possibly highly non-optimal
def get_height_data(E, xbound, outfile_name):
  # Get defining polynomial for E. This is given in homogeneous, so we dehomogenize
  poly = E.defining_polynomial()
  poly = poly(z = 1)

  K = E.base_field()
  Ok = K.maximal_order()

  # Output dictionaries
  data = {'height' : [], 'x_coefficients' : []}
  curve_data = {'curve' : E, 'Ok_basis' : Ok.gens()}

  # Iterate over tuples
  for tup in itertools.product(*[[-xbound..xbound] for i in [1..Ok.rank()]]):
    # Create x-value
    x_val = sum(a*b for (a,b) in zip(tup, Ok.gens()))
    f = poly(x = x_val)
    
    g = f.univariate_polynomial().change_ring(K)
    # If g is irreducible, we need to move to a quadratic extension of Q to find our point
    if g.is_irreducible():
      # Quadratic extension where g has a point
      M.<temp> = K.extension(g)
      L.<extgen> = M.absolute_field()
      from_L, to_L = L.structure()
      invars = E.ainvs()
      E_ext = EllipticCurve([L(invar) for invar in invars])

      # Compute height function & x point
      h = E_ext.height_function()
      x_point = to_L(M(x_val))

      # Compute height at each corresponding y-value
      gal1, gal2 = M.automorphisms()
      y_points = [to_L(gal1(temp)), to_L(gal2(temp))]
      # Compute both corresponding y-values
      for y in y_points:

        # Save data
        P = E_ext(x_point, y)
        # For torsion points, just add 0
        if P.order() < oo:
          data['height'].append(0)
          data['x_coefficients'].append(tup)
        else:
          # Get the degree of Q(x,y) to multiply by that
          f1 = x_point.minpoly(); f2 = y.minpoly()
          if f1.degree() > 1:
            D.<alpha> = NumberField(f1)
            # f2 is always quadratic, so we just multiply by 2 if irred
            if f2.change_ring(D).is_irreducible():
              multiplier = 2*D.degree()
            else:
              multiplier = D.degree()
          else:
            if f2.degree() == 1:
              multiplier = 1
            else:
              D.<alpha> = NumberField(f2)
              multiplier = D.degree()

          # Save height data
          data['height'].append(multiplier*h(E_ext(x_point, y)))
          data['x_coefficients'].append(tup)


    # In this case, we don't need to extend
    else:
      h = E.height_function()
      roots = g.roots()

      for i in [0,1]:
        x_point = x_val
        y = roots[i][0]
        P = E(x_point, y)
        if P.order() < oo:
          data['height'].append(0)
          data['x_coefficients'].append(tup)
        else:
          # Get the degree of Q(x,y) to multiply by that
          f1 = x_point.minpoly(); f2 = y.minpoly()
          if f1.degree() > 1:
            D.<alpha> = NumberField(f1)
            # f2 is always quadratic, so we just multiply by 2 if irred
            if f2.change_ring(D).is_irreducible():
              multiplier = 2*D.degree()
            else:
              multiplier = D.degree()
          else:
            if f2.degree() == 1:
              multiplier = 1
            else:
              D.<alpha> = NumberField(f2)
              multiplier = D.degree()

          data['height'].append(multiplier*h(E(x_val, roots[i][0])))
          data['x_coefficients'].append(tup)

  # Save height data
  df = pd.DataFrame(data = data)

  # Make heights floats because it is SO much more compatible with pandas
  df['height'] = df['height'].apply(lambda x : float(x))
  df.to_pickle(outfile_name + ".pkl")

  # Save curve data
  curve_df = pd.DataFrame(data = curve_data)
  curve_df.to_pickle(outfile_name + "-curve.pkl")