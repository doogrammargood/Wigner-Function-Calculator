
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
