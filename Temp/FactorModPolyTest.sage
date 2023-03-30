# This code computes one example: p = 83, l = 3, d_1 = -4 and d_2 = -31
# In this example, we see that the modular polynomial Phi_3 evaluated at 
# j(i) and j(h), where h is a root of the Hilbert Class Polynomial of Q(sqrt(-31))
# is divisible by only primes over 83 that lie over ONE of the primes over 83 in 
# the HCF for Q(sqrt(-31)).

QuadBase = NumberField(x^2 + 31, 't')
HCP = QuadBase.hilbert_class_polynomial()
HCF = QuadBase.extension(HCP, 'h')
h = HCF.gen()
RCP = algdep(elliptic_j(3/2 + 3*sqrt(-31)/2, prec = 100000), 12)
RCP = RCP.change_ring(HCF)
RCP0 = RCP.factor()[0][0]
RCF = HCF.extension(RCP0, 'r')
FinalField = RCF.extension(x^2 + 1, 'i')
R.<X,Y> = ZZ[]
MP3 = X^4 - X^3*Y^3 + 2232*X^3*Y^2 - 1069956*X^3*Y + 36864000*X^3 + 2232*X^2*Y^3 + 2587918086*X^2*Y^2 + 8900222976000*X^2*Y +452984832000000*X^2 - 1069956*X*Y^3 + 8900222976000*X*Y^2 - 770845966336000000*X*Y + 1855425871872000000000*X + Y^4 + 36864000*Y^3 + 452984832000000*Y^2 + 1855425871872000000000*Y
MP3 = MP3.change_ring(FinalField)
ModPolyideal = FinalField.ideal(MP3(1728, h))

# By looking at the norm, I see that 83^{32} exactly divides the norm of ModPolyideal
# Since the field is Galois, and the norm of the ideal 83 is 83^{48}, I think this can only happen
# if the norm of each prime over 83 is 83^2, and there are 5 of them that divide MP3(1728, h).
# but do all 5 of these lie over the same prime in HCF? 

# Actually, I figured out how to check this below:

HCFideal = HCF.ideal(83)
idealfactor0 = HCFideal.factor()[0][0]
idealfactor1 = HCFideal.factor()[1][0]
idealfactor2 = HCFideal.factor()[2][0]

primeFromHCF0 = FinalField.ideal(idealfactor0)
primeFromHCF1 = FinalField.ideal(idealfactor1)
primeFromHCF2 = FinalField.ideal(idealfactor2)

ModPolyideal.is_coprime(primeFromHCF0)
ModPolyideal.is_coprime(primeFromHCF1)
ModPolyideal.is_coprime(primeFromHCF2)

# We can also check that replacing h by a Galois conjugate simply permutes which prime over 83 in the HCF for Q(sqrt(-31))
# is under all the primes over 83 that divide the modular polynomial

otherJs = h.galois_conjugates(HCF)

ModPolyideal0 = FinalField.ideal(MP3(1728, otherJs[0]))

ModPolyideal0.is_coprime(primeFromHCF0)
ModPolyideal0.is_coprime(primeFromHCF1)
ModPolyideal0.is_coprime(primeFromHCF2)

ModPolyideal2 = FinalField.ideal(MP3(1728, otherJs[2]))

ModPolyideal2.is_coprime(primeFromHCF0)
ModPolyideal2.is_coprime(primeFromHCF1)
ModPolyideal2.is_coprime(primeFromHCF2)
