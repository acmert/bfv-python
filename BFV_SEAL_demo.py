from generate_prime import *
from BFV import *
from helper import *

from random import randint

# Parameter generation (or enter manually)

# Determine n and bit-size of q, then find a
# q satisfying the condition: q = 1 (mod 2n)
#
# Based on n and q, generate NTT parameters

# Parameter (generated or pre-defined)
PC = 1 # 0: generate -- 1: pre-defined

n   = 0
q   = 0
psi = 0
psiv= 0
w   = 0
wv  = 0

g   = 0
# g (gamma) parameter is a constant 61-bit integer in SEAL
# For a K-bit q, gamma could be a K-bit integer larger than q and co-prime to q, t

# Determine t
t = 4

print("--- Starting BFV (SEAL) Demo")

if PC:
    n = 2048

    q_bit = 37

    # Parameters for q
    q   = 137438691329
    psi = 22157790
    psiv= modinv(psi,q)
    w   = pow(psi,2,q)
    wv  = modinv(w,q)

    g   = 0x1fffffffffc80001 # original integer from SEAL

    print("* q is calculated.")
    print("* n   : {}".format(n))
    print("* q   : {}".format(q))
    print("* Parameters are calculated.")
    print("* w   : {}".format(w))
    print("* wv  : {}".format(wv))
    print("* psi : {}".format(psi))
    print("* psiv: {}".format(psiv))
    print("* g   : {}".format(g))

else:
    n = 4096

    q_bit = 25

    # calculate q and qnp
    wfound = False
    while(not(wfound)):
        q = generate_large_prime(q_bit)
        # check q = 1 (mod 2n)
        while (not ((q % (2*n)) == 1)):
            q = generate_large_prime(q_bit)

        # generate NTT parameters
        for i in range(2,q-1):
            wfound = isrootofunity(i,2*n,q)
            if wfound:
                psi = i
                psiv= modinv(psi,q)
                w   = pow(psi,2,q)
                wv  = modinv(w,q)
                break

    # calculate g
    while(1):
        g = randint(q,2**q_bit - 1)
        if((g>q) and (gcd(g,q) == 1) and (gcd(g,t) == 1)):
            break

    print("* q is calculated.")
    print("* n   : {}".format(n))
    print("* q   : {}".format(q))
    print("* Parameters are calculated.")
    print("* w   : {}".format(w))
    print("* wv  : {}".format(wv))
    print("* psi : {}".format(psi))
    print("* psiv: {}".format(psiv))
    print("* g   : {}".format(g))

# Determine B (bound of distribution X)
sigma = 3.1
B = int(10*sigma)

# Determine T, p (relinearization)
T = 256
p = 16

# Generate tables
w_table    = [1]*n
wv_table   = [1]*n
psi_table  = [1]*n
psiv_table = [1]*n
for i in range(1,n):
    w_table[i]    = ((w_table[i-1]   *w)    % q)
    wv_table[i]   = ((wv_table[i-1]  *wv)   % q)
    psi_table[i]  = ((psi_table[i-1] *psi)  % q)
    psiv_table[i] = ((psiv_table[i-1]*psiv) % q)

qnp = [w_table,wv_table,psi_table,psiv_table]

# Generate BFV evaluator
Evaluator = BFV(n, q, t, B, qnp)

# Generate Keys
Evaluator.SecretKeyGen()
Evaluator.PublicKeyGen()
Evaluator.EvalKeyGenV1(T)
Evaluator.EvalKeyGenV2(p)

# print system parameters
print(Evaluator)

# Generate random message
n1 = randint(0,q-1)

print("--- Random integer n1 is generated.")
print("* n1: {}".format(n1))
print("")

# Encode random messages into plaintext polynomials
m1 = Evaluator.Encode(n1)

print("--- n1 is encoded as polynomial m1(x).")
print("* m1(x): {}".format(m1))
print("")

# Encrypt message
ct1 = Evaluator.EncryptionSEAL(m1)

print("--- m1 is encrypted as ct1.")
print("* ct1[0]: {}".format(ct1[0]))
print("* ct1[1]: {}".format(ct1[1]))
print("")

# Decrypt message
mt = Evaluator.DecryptionSEAL(ct1,g)

print("--- ct1 is decrypted as mt(x).")
print("* mt(x): {}".format(mt))
print("")

# Decode decrypted polynomial into integer
nr = Evaluator.Decode(mt)

print("--- mt(x) is decoded as integer nr.")
print("* nr: {}".format(nr))
print("")

if (m1 == mt) and (nr == n1):
    print("--- Encryption-Decryption works.")
else:
    print("--- Encryption-Decryption does not work.")

#
