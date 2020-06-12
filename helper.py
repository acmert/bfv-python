
# Modular inverse of an integer
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('Modular inverse does not exist')
    else:
        return x % m

# GCD of two integers
def gcd(n1, n2):
    a = n1
    b = n2
    while b != 0:
        a, b = b, a % b
    return a

# Bit-Reverse integer
def intReverse(a,n):
    b = ('{:0'+str(n)+'b}').format(a)
    return int(b[::-1],2)

# Bit-Reversed index
def indexReverse(a,r):
    n = len(a)
    b = [0]*n
    for i in range(n):
        rev_idx = intReverse(i,r)
        b[rev_idx] = a[i]
    return b

# Reference Polynomial Multiplication
# with f(x) = x^n + 1
def RefPolMul(A, B, M):
    C = [0] * (2 * len(A))
    D = [0] * (len(A))
    for indexA, elemA in enumerate(A):
        for indexB, elemB in enumerate(B):
            C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % M

    for i in range(len(A)):
        D[i] = (C[i] - C[i + len(A)]) % M
    return D

# Reference Polynomial Multiplication (w/ modulus)
# with f(x) = x^n + 1
def RefPolMulv2(A, B):
    C = [0] * (2 * len(A))
    D = [0] * (len(A))
    for indexA, elemA in enumerate(A):
        for indexB, elemB in enumerate(B):
            C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB)

    for i in range(len(A)):
        D[i] = (C[i] - C[i + len(A)])
    return D

# Check if input is m-th (could be n or 2n) primitive root of unity of q
def isrootofunity(w,m,q):
    if pow(w,m,q) != 1:
        return False
    elif pow(w,m//2,q) != (q-1):
        return False
    else:
        v = w
        for i in range(1,m):
            if v == 1:
                return False
            else:
                v = (v*w) % q
        return True
