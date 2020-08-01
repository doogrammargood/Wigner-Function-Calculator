def dig(a,p,i):
    #returns the ith digit of a in base p
    return int(((a%(p**(i+1)) )- (a%(p**(i)) ) ) / p**i)
def to_base(a,base):
    to_return = []
    current = a
    smallest_digit = dig(a,base,0)
    while(current > 0):
        to_return.append(smallest_digit)
        current = (current-smallest_digit)//base
        smallest_digit = dig(current,base,0)
    return to_return

def is_prime(n) : #from https://www.geeksforgeeks.org/python-program-to-check-whether-a-number-is-prime-or-not/
    if (n <= 1) :
        return False
    if (n <= 3) :
        return True
    if (n % 2 == 0 or n % 3 == 0) :
        return False
    i = 5
    while(i * i <= n) :
        if (n % i == 0 or n % (i + 2) == 0) :
            return False
        i = i + 6
    return True

def modinv(a, m):
    a = a % m;
    for x in range(1, m) :
        if ((a * x) % m == 1) :
            return x
    return 1

def next_prime(x): #returns the smallest prime geq than x.
    while(True):
        if is_prime(x):
            return x
        x += 1
def previous_prime(x): #returns the greatest prime leq than x
    while True:
        if is_prime(x):
            return x
        x -= 1
def legendre_symbol(a, n):#credit John D. Cook
    #print(n,a)
    assert(n > a > 0 and n%2 == 1)
    t = 1
    while a != 0:
        while a % 2 == 0:
            a /= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == n % 4 == 3:
            t = -t
        a %= n
    if n == 1:
        return t
    else:
        return 0
def prime_factorization(n): #This function assumes n < 510,510
    current_n = n
    current_prime = 2
    prime_list = []
    to_return = [] #list of pairs (current, current_prime)
    to_return = [0]*7
    while (current_prime <= current_n ** 0.5):
        while(current_n % current_prime == 0):
            if not current_prime == prime_list:
                prime_list.append(current_prime)
                to_return.append(1)
            else:
                to_return[-1] += 1
            current_n = current_n // current_prime
        current_prime = next_prime(current_prime+1)
    to_return = [tr for tr in to_return if tr!=0]
    if current_n != 1: #current_n must be prime itself,
        to_return.append(1)
        prime_list.append(current_n)

    return zip(prime_list, to_return)

def maximal_proper_factors(n):
    return [n//p[0] for p in prime_factorization(n)]
