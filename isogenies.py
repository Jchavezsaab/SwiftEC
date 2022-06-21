from decimal import DivisionByZero
from Xencoding import *

""" 
This file includes the parameters and isogeny implementations for cases
where we want to compose the encoding with an isogeny, which allows us
to reach even more curves. This should never be needed for a=0 curves.
"""

if curve == "P25519":

    # Codomain curve:
    A = F(486662)
    A2 = (A+2)/4
    E1 = EllipticCurve([0, A, 0, 1, 0])
    # Index of the subgroup
    index = 8
    # Precomputed parameters for 2-isogeny:
    x2 = F(57896044618658097711785492504343953926634992332820282019728792003956553140061)
    t = F(306944513303040)
    # Precomputed parameters for isomorphism:
    u2 = F(23721296059033526145801000401085370011607392691919421105305546723843314752618)
    r = F(38597363079105398474523661669562635951089994888546854679819194669304376384412)
    dual_u2 = F(36)
    dual_r = F(5839944)

    def isogeny(x,z):

        # Apply a 2-isogeny: (x:z) -> ( (x(x-zx2)+tz^2) : z(x-x2z) )
        tz2 = fp_sqr(z)
        tz2 = fp_mul(tz2, t)
        zprime = fp_mul(x2, z)
        zprime = fp_sub(x, zprime)
        xprime = fp_mul(zprime, x)
        xprime = fp_add(xprime, tz2)
        zprime = fp_mul(zprime, z)

        # Apply an isomorphism: (x:z) -> ( u^2x+r*z : z)
        tz2 = fp_mul(r, zprime)
        xprime = fp_mul(xprime, u2)
        xprime = fp_add(xprime, tz2)

        # Multiplication by 4
     
        t0 = fp_sub(xprime, zprime)
        t1 = fp_add(xprime, zprime)
        zprime = fp_sqr(t0)
        t1 = fp_sqr(t1)
        xprime = fp_mul(zprime, t1)
        t1 = fp_sub(t1, zprime)
        t0 = fp_mul(A2, t1)
        zprime = fp_add(zprime, t0)
        zprime = fp_mul(zprime, t1)
     
        t0 = fp_sub(xprime, zprime)
        t1 = fp_add(xprime, zprime)
        zprime = fp_sqr(t0)
        t1 = fp_sqr(t1)
        xprime = fp_mul(zprime, t1)
        t1 = fp_sub(t1, zprime)
        t0 = fp_mul(A2, t1)
        zprime = fp_add(zprime, t0)
        zprime = fp_mul(zprime, t1)

        return xprime,zprime

    def dual_isogeny(x,z):

        # Apply a 2-isogeny: (x:z) -> ( (x^2+1) : zx )
        t0 = fp_sqr(z)
        xprime = fp_sqr(x)
        xprime = fp_add(xprime, t0)
        zprime = fp_mul(x, z)
        
        # Apply an isomorphism: (x:z) -> ( u^2x+r*z : z)
        t0 = fp_mul(dual_r, zprime)
        xprime = fp_mul(xprime, dual_u2)
        xprime = fp_add(xprime, t0)

        # Add a random 8-torsion poin
        return xprime,zprime


    def get_y2(x):
        y = fp_add(x, A)
        y = fp_mul(y, x)
        y = fp_mul(y, x)
        y = fp_add(y, x)
        return y

# elif curve == 'MyCurve':
#     """ For any new curve, define the objects
#         -E1
#         -index
#         -isogeny(x,y)
#         -dual_isogeny(x,y)
#         -get_y2(x)
#     """

else:
    raise NotImplementedError("Isogeny for curve "+curve+" must be defined in isogenies.py")


### The remaining functions are common to all curves (no need to edit)

def decisog(u, t, s):
    x,z = Xdecode(u,t)
    x,z = isogeny(x,z)
    try:
        z = fp_inv(z)
    except DivisionByZero:
        raise PointAtInfinity
    x = fp_mul(x, z)
    y = get_y2(x)
    y = fp_sqrt(y)
    negy = fp_neg(y)
    c = (int(y) & 1) ^ (int(s) & 1)
    y, negy = fp_cswap(c, y, negy)
    return x, y

def encisog(x, y):
    z = F(1)
    x,z = dual_isogeny(x, z)
    s = int(y) & 1
    u,t =  Xencode(x,z)
    return u,t,s

def Xdecisog(u, t):
    x,z = Xdecode(u, t)
    x,z = isogeny(x, z)
    return x,z

def Xencisog(x, z):
    x,z = dual_isogeny(x, z)
    u,t =  Xencode(x,z)
    return u,t