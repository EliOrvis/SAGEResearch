#### This file contains miscellaneous arithmetic functions that are helpful for working with any code base.

# Function to give a random prime in a particular range with a particular modulo condition. Code from David Lowry-Duda

def random_prime_in_a_mod_b(lowerbound, upperbound, a, b, proof=False):
    """
    Returns a random prime in [lowerbound, upperbound]
    that is congruent to a mod b.
    """
    if not gcd(a, b) == 1:
        raise ValueError("a and b are not coprime")
    if not lowerbound < upperbound:
        raise ValueError("lowerbound is not smaller than upperbound")
    if proof:
        prime_test = is_prime
    else:
        prime_test = is_pseudoprime
    while True:
        n = randint(lowerbound // b, upperbound // b)
        p = n * b + a
        if p < lowerbound or p > upperbound:
            continue
        if prime_test(p):
            return p