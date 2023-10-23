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

# CHOSEN CIPHERTEXT ATTACK

print()
print("CHOSEN CIPHERTEXT ATTACK")
print()

print()
while True:
    keygen = input("Press enter to begin the Chosen Ciphertext Attack")
    if keygen.upper() == '':
        break
    else:
        sys.exit()
print()

start_time = time.time()

def choose_c():
    while True:
        c = random.randint(q//4,q//2)
        if c % p == 0:
            return c
    
c = choose_c()

print()
print("Your c is: ", c)
print()

#df = d+1
#dg = d
#pred_coll = (((2*df)-1)*dg)/(N)

#print()
#print("The predicted number of collisions is", pred_coll)
#print()

i = 0

while True:
    
    print()
    print("------------")
    print("Iteration", i)
    print("------------")
    print()
    
    chosen_ciphertext = c*h + c*(Poly(x**i,x).set_domain(ZZ))

    print("Your chosen ciphertext is: ")
    print(str(chosen_ciphertext.as_expr()).replace("**", "^"))
    print()

    def decryption_alg(ciphertext):
        a = Poly([((coeff % q) + q) % q for coeff in (((ciphertext * f) % R).trunc(q)).all_coeffs()], x).set_domain('ZZ')
        ahat = a.trunc(q)
        b = Poly([((coeff % p) + p) % p for coeff in (((ahat * F_p) % R).trunc(p)).all_coeffs()], x).set_domain('ZZ')
        return b

    decrypted_message = decryption_alg(chosen_ciphertext)

    print("Feeding your chosen ciphertext to the decryption algorithm yields: ")
    print(str(decrypted_message.as_expr()).replace("**", "^"))
    print()
    
    def red_coeffs(coeff):
        while coeff != 0 and coeff != c and coeff != -c and coeff != 2*c and coeff != -2*c:
            coeff = coeff - q
            if coeff < -2*c:
                break
        while coeff != 0 and coeff != c and coeff != -c and coeff != 2*c and coeff != -2*c:
            coeff = coeff + q
            if coeff > 2*c:
                break
        return coeff

    a = ((chosen_ciphertext * f) % R)

    a_red = Poly([red_coeffs(coeff) for coeff in (a).all_coeffs()], x).set_domain('ZZ')

    def kcoeffs(poly):
        coeffs = poly.all_coeffs()
        kcoeffs = []
        for coeff in coeffs:
            if coeff == 2*c or coeff == -2*c:
                kcoeffs.append(coeff/(2*c))
            else:
                kcoeffs.append(0)
        return kcoeffs

    k = Poly(kcoeffs(a_red), x).set_domain('ZZ')

    print("k =", str(k.as_expr()).replace("**", "^"))
    print()
    
    try:
        minverse = invert(decrypted_message, R, domain=GF(p))
    except NotInvertible:
        i += 1
        pass
    else:
        f_prime = (((invert(decrypted_message, R, domain=GF(p))) * (-q*k))%R).trunc(p)

        print()
        print("f' =", str(f_prime.as_expr()).replace("**", "^"))
        print()
        break
    
F_prime_p = (Poly([((coeff % p) + p) % p for coeff in invert(f_prime, R, domain=GF(p)).all_coeffs()], x).set_domain('ZZ')).trunc(p)
a = Poly([((coeff % q) + q) % q for coeff in (((e * f_prime) % R).trunc(q)).all_coeffs()], x).set_domain('ZZ')
ahat = a.trunc(q)
b = Poly([((coeff % p) + p) % p for coeff in (((ahat * F_prime_p) % R).trunc(p)).all_coeffs()], x).set_domain('ZZ')

print()
print("Using f', we were able to recover: ")
print()
print("m =", str(b.trunc(p).as_expr()).replace("**", "^"))
print()

if b.trunc(p) == m:
    print("The message was recovered and the attack was successful!")
    end_time = time.time()
    runtime = end_time - start_time
    print()
    print("The runtime of this attack was", runtime, "seconds")
    sys.exit()
else:
    print("Something went wrong!")




