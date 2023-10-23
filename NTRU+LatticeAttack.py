import sys
import math
from math import gcd
from sympy import ZZ, Poly, invert, GF, NotInvertible
from sympy.abc import x
import numpy as np
import random
import time

# NTRU CRYPTOSYSTEM

print()
print("NTRU CRYPTOSYSTEM")
print()

# Parameter Creation

print()
while True:
    keygen = input("Press enter to begin Parameter Creation")
    if keygen.upper() == '':
        break
    else:
        sys.exit()
print()

def is_prime(number):
    if number <= 1:
        return False
    if number == 2:
        return True
    sqrt = int(math.sqrt(number)) + 1
    for divisor in range(2, sqrt):
        if number % divisor == 0:
            return False
    return True

def generate_prime(x,y):
    while True:
        prime = random.randint(x, y)  
        if is_prime(prime):
            return prime
        
print()
print("Enter values for the public parameters.")
print()
print("If you wish to have the value randomly generated, leave it blank by pressing enter.")
print()

N = input("Choose a prime value for N: ")

if N == '':
    N = int(generate_prime(3,20))
else:
    N = int(N)
    while is_prime(N) == False:
        if is_prime(N) == False:
            print("N has to be prime")
        N = int(input("Choose a PRIME value for N: "))
print("Your value of N is", N)
print()

d = input("Choose a value for d: ")

if d == '':
    d = int(random.randint(1, ((-1/2)+(N/2))))
else:
    d = int(d)
print("Your value of d is", d)
print()

p = input("Choose a prime value for p: ")

if p == '':
    p = int(generate_prime(3,20))
else:
    p = int(p)
    while is_prime(p) == False:
        if is_prime(p) == False:
            print("p has to be prime")
        p = int(input("Choose a PRIME value for p: "))
print("Your value of p is", p)
print()

q = input(f"Choose a prime value for q such that gcd({p},q) = gcd({N},q) = 1 and q > {(((6*d)+1)*p)}: ")

if q == '':
    q = int(generate_prime((((6*d)+1)*p) + 1, 1122))
    while gcd(p,q) != 1 or gcd(N,q) != 1 or q <= (((6*d)+1)*p):
        q = int(generate_prime((((6*d)+1)*p) + 1, 1122))
else:
    q = int(q)
    while is_prime(q) == False or gcd(p,q) != 1 or gcd(N,q) != 1 or q <= (((6*d)+1)*p):
        print()
        if is_prime(q) == False:
            print("q has to be prime")
        if gcd(p,q) != 1:
            print(f"The gcd({p},q) has to be 1")
        if gcd(N,q) != 1:
            print(f"The gcd({N},q) has to be 1")
        if q <= (((6*d)+1)*p):
            print(f"q has to be > {(((6*d)+1)*p)}")
        print()
        q = int(input(f"Choose a prime value for q such that gcd({p},q) = gcd({N},q) = 1 and q > {(((6*d)+1)*p)}: "))
print("Your value of q is", q)
print()

print()
print("Your public parameters (N,p,q,d) are:", (N,p,q,d))
print()
        
# Key Generation

print()
while True:
    keygen = input("Press enter to begin Key Generation")
    if keygen.upper() == '':
        break
    else:
        sys.exit()
print()

R = Poly(x ** N - 1, x).set_domain('ZZ')

while True:
    f = Poly(np.random.permutation(np.concatenate((np.ones(d+1),-np.ones(d),np.zeros(N-((2*d)+1))))),x).set_domain('ZZ')
    f = f % R
    
    try:
        pinverse = invert(f, R, domain=GF(p))
    except NotInvertible:
        pass
    
    try:
        qinverse = invert(f, R, domain=GF(q))
    except NotInvertible:
        pass
    else:
        break
     
g = Poly(np.random.permutation(np.concatenate((np.ones(d),-np.ones(d),np.zeros(N-(2*d))))),x).set_domain('ZZ')

f = f % R
g = g % R

print()
print("Your randomly generated private key (f,g) is: ")
print()
print("f =", str(f.as_expr()).replace("**", "^"))
print()
print("g =", str(g.as_expr()).replace("**", "^"))
print()


F_p = Poly([((coeff % p) + p) % p for coeff in invert(f, R, domain=GF(p)).all_coeffs()], x).set_domain('ZZ')

F_q = Poly([((coeff % q) + q) % q for coeff in invert(f, R, domain=GF(q)).all_coeffs()], x).set_domain('ZZ')

h = Poly([((coeff % q) + q) % q for coeff in (((F_q * g) % R).trunc(q)).all_coeffs()], x).set_domain('ZZ')

print("Your public key h is: ")
print("h =", str(h.as_expr()).replace("**", "^"))
print()

# Encryption

print()
while True:
    keygen = input("Press enter to begin Encryption")
    if keygen.upper() == '':
        break
    else:
        sys.exit()
print()

print()
while True:
    choice = input("Would you like your message to be numeric (N) or text (T)? ")
    if choice.upper() == 'N':
        print()
        while True:
            print("You will now enter the coefficients of your polynomial m.")
            print("Separate entries with a space.")
            print(f"m must be a center lift of a polynomial in R{p}.")
            print(f"If you enter more than {N} values, coefficients may change due to reduction in R.")
            print(f"m should have coefficients between {-p/2} and {p/2}.")
            print()
            message = input("Enter your coefficients here: ")
            coeffs = [int(num) for num in message.split()]
            m = (Poly(coeffs,x).set_domain(ZZ) % R)
            def check_coeffs(coeffs):
                k = 0
                for entry in coeffs:
                    if entry <= (-p/2) or entry > (p/2):
                        k += 1
                if k > 0:
                    return False
                elif k == 0:
                    return True
            if check_coeffs(m.all_coeffs()) == True:
                print()
                print("Your message m is: ")
                print("m =", str(m.as_expr()).replace("**", "^"))
                break
            else:
                print()
                print("Your polynomial m =", str(m.as_expr()).replace("**", "^"))
                print(f"is not a center lift of a polynomial in R{p}.")
                print()
        break
    if choice.upper() == 'T':
        print()
        message = input("Enter your text message m: ")
        coeffs = [ord(char) for char in message]
        m = (Poly(coeffs,x).set_domain(ZZ) % R).trunc(p)
        print()
        print("Your message m converted into a polynomial is: ")
        print("m =", str(m.as_expr()).replace("**", "^"))
        break
    else:
        print("Invalid choice. Choose either N or T.")
print()


r = Poly(np.random.permutation(np.concatenate((np.ones(d),-np.ones(d),np.zeros(N-(2*d))))),x).set_domain(ZZ)

e = Poly([((coeff % q) + q) % q for coeff in (((p * r * h + m) % R).trunc(q)).all_coeffs()], x).set_domain('ZZ')

print()
print("Your ciphertext e is: ")
print("e =", str(e.as_expr()).replace("**", "^"))
print()

# LATTICE ATTACK

print()
print("LATTICE ATTACK")
print()

print()
while True:
    keygen = input("Press enter to begin the Lattice Attack")
    if keygen.upper() == '':
        break
    else:
        sys.exit()
print()

start_time = time.time()

if len(h.all_coeffs()) < N:
    hcoeffs = h.all_coeffs()
    num_zeros_to_add = N - len(h.all_coeffs())
    zeros_to_add = [0] * num_zeros_to_add
    hcoeffs = zeros_to_add + hcoeffs
else:
    hcoeffs = h.all_coeffs()

def create_circulant_matrix(lst):
    n = len(lst)
    return np.array([lst[(j-i) % n] for i in range(n) for j in range(n)]).reshape(n, n)

upper_right = create_circulant_matrix(hcoeffs[::-1])
upper_left = np.identity(N)
lower_left = np.zeros((N, N))
lower_right = q * np.identity(N)
M = np.vstack(((np.hstack((upper_left, upper_right))), (np.hstack((lower_left, lower_right))))).astype(int)

print()
print("The NTRU matrix is: ")
print("M =", M)
print()

def GramSchmidt(v):
    (n,m) = v.shape
    u = np.zeros((n, m), dtype=float)
    u[0] = v[0]
    def coef(v1, v2):
        return np.dot(v2,v1)/np.dot(v1,v1)
    def proj(v1, v2):
        return np.dot(coef(v1,v2),v1)
    for i in range(1,n):
        tot = 0
        for j in range(1,i+1):
            tot += proj(u[j-1],v[i])
        val = np.subtract(v[i],tot)
        u[i] = val
    return u

def LLL(v):
    k = 2
    
    (n,m) = v.shape
    
    while k <= n:
        V = GramSchmidt(v)
    
        def mu(i,j):
            return (np.dot(v[i-1],V[j-1]))/(np.linalg.norm(V[j-1])**2)
        
        for j in range(k-1,0,-1):
            m = mu(k,j)
            v[k-1] = v[k-1] - (round(m)*v[j-1])
         
        if (np.linalg.norm(V[k-1])**2) >= ((3/4)-((mu(k,k-1))**2))*(np.linalg.norm(V[k-2])**2):
            k = k + 1 
        else:
            v[[k-2, k-1]] = v[[k-1, k-2]]
            k = max(k-1,2) 
    
    return v

Mred = LLL(M)

print()
print("The LLL reduced NTRU matrix is: ")
print("Mred =", Mred)
print()


def find_smallest_nonzero_vector(matrix):
    smallest_nonzero_vector = None
    smallest_magnitude = None

    for row_vector in matrix:
        nonzero_indices = np.nonzero(row_vector)[0]
        first_half_indices = nonzero_indices[nonzero_indices < len(row_vector)//2]
        second_half_indices = nonzero_indices[nonzero_indices >= len(row_vector)//2]
        
        if len(first_half_indices) > 0 and len(second_half_indices) > 0:
            current_magnitude = np.linalg.norm(row_vector)
            if smallest_nonzero_vector is None or current_magnitude < smallest_magnitude:
                smallest_nonzero_vector = row_vector
                smallest_magnitude = current_magnitude
    
    return smallest_nonzero_vector
    
smallest_vector = find_smallest_nonzero_vector(Mred)

print()
print("The smallest nonzero vector is: ")
print(smallest_vector)
print()

fcoeffs = smallest_vector[:(len(smallest_vector) // 2)]
gcoeffs = smallest_vector[(len(smallest_vector) // 2):]

fprime = Poly(fcoeffs[::-1], x).set_domain('ZZ')
gprime = Poly(gcoeffs[::-1], x).set_domain('ZZ')

print()
print("The recovered keys are: ")
print()
print("f' =", str(fprime.as_expr()).replace("**", "^"))
print()
print("g' =", str(gprime.as_expr()).replace("**", "^"))
print()

Fprime_p = Poly([((coeff % p) + p) % p for coeff in invert(fprime, R, domain=GF(p)).all_coeffs()], x).set_domain('ZZ')
Fprime_q = Poly([((coeff % q) + q) % q for coeff in invert(fprime, R, domain=GF(q)).all_coeffs()], x).set_domain('ZZ')    

aprime = Poly([((coeff % q) + q) % q for coeff in (((e * fprime) % R).trunc(q)).all_coeffs()], x).set_domain('ZZ')
ahatprime = aprime.trunc(q)
bprime = Poly([((coeff % p) + p) % p for coeff in (((ahatprime * Fprime_p) % R).trunc(p)).all_coeffs()], x).set_domain('ZZ')

print()
print("Using f' and g', we were able to recover: ")
print()
print("m =", str(bprime.trunc(p).as_expr()).replace("**", "^"))
print()

if bprime.trunc(p) == m:
    print("The message was recovered and the attack was successful!")
else:
    print("Something went wrong!")

end_time = time.time()

runtime = end_time - start_time

print()
print("The runtime of this attack was", runtime, "seconds")














