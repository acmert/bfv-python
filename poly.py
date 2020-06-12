# Poly

from random import randint
from ntt import *

class Poly:
    def __init__(self, n, q, np=[0,0,0,0]):
        self.n = n
        self.q = q
        self.np= np # NTT parameters: [w,w_inv,psi,psi_inv]
        self.F = [0]*n
        self.inNTT = False
    #
    def randomize(self, B, domain=False):
        self.F = [randint(0, B-1) for i in range(self.n)]
        self.inNTT = domain
    #
    def __str__(self):
        pstr = str(self.F[0])
        tmp = min(self.n,8)

        for i in range(1,tmp):
            pstr = pstr+" + "+str(self.F[i])+"*x^"+str(i)

        if self.n > 8:
            pstr = pstr + " + ..."
        return pstr
    #
    def __add__(self, b):
        if self.inNTT != b.inNTT:
            raise Exception("Polynomial Addiditon: Inputs must be in the same domain.")
        elif self.q != b.q:
            raise Exception("Polynomial Addiditon: Inputs must have the same modulus")
        else:
            c = Poly(self.n, self.q, self.np)
            c.F = [(x+y)%self.q for x,y in zip(self.F,b.F)]
            c.inNTT = self.inNTT
            return c
    #
    def __sub__(self, b):
        if self.inNTT != b.inNTT:
            raise Exception("Polynomial Subtraction: Inputs must be in the same domain.")
        elif self.q != b.q:
            raise Exception("Polynomial Subtraction: Inputs must have the same modulus")
        else:
            c = Poly(self.n, self.q, self.np)
            c.F = [(x-y)%self.q for x,y in zip(self.F,b.F)]
            c.inNTT = self.inNTT
            return c
    #
    def __mul__(self, b):
        if self.inNTT != b.inNTT:
            raise Exception("Polynomial Multiplication: Inputs must be in the same domain.")
        elif self.q != b.q:
            raise Exception("Polynomial Multiplication: Inputs must have the same modulus")
        else:
            """
            Assuming both inputs in POL/NTT domain
            If in NTT domain --> Coeff-wise multiplication
            If in POL domain --> Full polynomial multiplication
            """
            c = Poly(self.n, self.q, self.np)
            if self.inNTT == True and b.inNTT == True:
                c.F = [((x*y)%self.q) for x,y in zip(self.F,b.F)]
                c.inNTT = True
            else:
                # x1=self*psi, x2=b*psi
                # x1n = NTT(x1,w), x2n = NTT(x2,w)
                # x3n = x1n*x2n
                # x3 = INTT(x3n,w_inv)
                # c = x3*psi_inv

                w_table    = self.np[0]
                wv_table   = self.np[1]
                psi_table  = self.np[2]
                psiv_table = self.np[3]

                s_p = [(x*psi_table[pwr])%self.q for pwr,x in enumerate(self.F)]
                b_p = [(x*psi_table[pwr])%self.q for pwr,x in enumerate(b.F)]
                s_n = NTT(s_p,w_table,self.q)
                b_n = NTT(b_p,w_table,self.q)
                sb_n= [(x*y)%self.q for x,y in zip(s_n,b_n)]
                sb_p= INTT(sb_n,wv_table,self.q)
                sb  = [(x*psiv_table[pwr])%self.q for pwr,x in enumerate(sb_p)]

                c.F = sb
                c.inNTT = False
            return c
    #
    def __mod__(self,base):
        b = Poly(self.n, self.q, self.np)
        b.F = [(x%base) for x in self.F]
        b.inNTT = self.inNTT
        return b
    #
    def __round__(self):
        b = Poly(self.n, self.q, self.np)
        b.F = [round(x) for x in self.F]
        b.inNTT = self.inNTT
        return b
    #
    def __eq__(self, b):
        if self.n != b.n:
            return False
        elif self.q != b.q:
            return False
        else:
            for i,j in zip(self.F,b.F):
                if i != j:
                    return False
            return True
    #
    def __neg__(self):
        b = Poly(self.n, self.q, self.np)
        b.F = [((-x) % self.q) for x in self.F]
        b.inNTT = self.inNTT
        return b
    #
    def toNTT(self):
        b = Poly(self.n, self.q, self.np)
        if self.inNTT == False:
            b.F = NTT(self.F,self.np[0],self.q)
            b.inNTT = True
        else:
            b.F = [x for x in self.F]
            b.inNTT = True
        return b
    #
    def toPOL(self):
        b = Poly(self.n, self.q, self.np)
        if self.inNTT == False:
            b.F = [x for x in self.F]
            b.inNTT = False
        else:
            b.F = INTT(self.F,self.np[1],self.q)
            b.inNTT = False
        return b
#
