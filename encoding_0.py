from fp import *
from random import random

def decode(u, t, s):

    #Evalaute point in conic X(u,t),Y(u,t)
    X1 = fp_sqr(u)
    X1 = fp_mul(X1, u)
    Y1 = fp_sqr(t)
    X1 = fp_add(X1, b)
    X1 = fp_sub(X1, Y1)
    Y1 = fp_add(Y1, Y1)
    Y1 = fp_add(Y1, X1)
    Z1 = fp_mul(u, X[0])   # X[0] is just sqrt(-3)
    X1 = fp_mul(X1, Z1)
    Z1 = fp_mul(Z1, t)
    Z1 = fp_add(Z1, Z1)

    #print("Test X1,Y1,Z1",X1**2+(3*u**2+4*a)*Y1**2+(u**3+a*u+b)*Z1**2==0)

    # Compute projective point in surface S
    #   y = (2Y1)^2
    #   v = X1*Z1 - u*Y1*Z1
    #   w = 2*Y1*Z1
    y = fp_add(Y1, Y1)
    y = fp_sqr(y)
    v = fp_mul(Y1, u)
    v = fp_sub(X1, v)
    v = fp_mul(v, Z1)
    w = fp_mul(Y1, Z1)
    w = fp_add(w, w)

    #print("Tests v,y,w", y**2*(w**2*u**2 + w*v*u + v**2 + w**2*a) == -w**4*(u**3 + a*u + b))



    # Compute affine point in V
    #   x1 = v/w
    #   x2 = -u-v/w
    #   x3 = u + y^2/w^2
    try:
        w = fp_inv(w)
    except ZeroDivisionError:
        raise PointAtInfinity
    x1 = fp_mul(v, w)
    x2 = fp_add(u, x1)
    x2 = fp_neg(x2)
    x3 = fp_mul(y, w)
    x3 = fp_sqr(x3)
    x3 = fp_add(u, x3)

    # Compute g(x_i)
    y21 = fp_sqr(x1)
    y21 = fp_add(y21, a)
    y21 = fp_mul(y21, x1)
    y21 = fp_add(y21, b)

    y22 = fp_sqr(x2)
    y22 = fp_add(y22, a)
    y22 = fp_mul(y22, x2)
    y22 = fp_add(y22, b)

    y23 = fp_sqr(x3)
    y23 = fp_add(y23, a)
    y23 = fp_mul(y23, x3)
    y23 = fp_add(y23, b)

    # Find the square
    c2 = fp_jacobi(y22)
    c3 = fp_jacobi(y23)
    x1, x2 = fp_cswap(c2, x1, x2)
    y21, y22 = fp_cswap(c2, y21, y22)
    x1, x3 = fp_cswap(c3, x1, x3)
    y21, y23 = fp_cswap(c3, y21, y23)

    # Find the square-root and choose sign
    y21 = fp_sqrt(y21)
    y22 = fp_neg(y21)
    c1 = ((int(y21) % 2) ^ (int(s) % 2))
    y21, y22 = fp_cswap(c1, y21, y22)
    return x1, y21


def encode(x, y):
    u = F.random_element()
    case = ceil(4*random())
    
    # When x is x1
    if case == 1:
        
        # Check that x_2 doesnt yield a square
        v = fp_neg(x)
        v = fp_sub(v, u)
        y2 = fp_sqr(v)
        y2 = fp_add(y2, a)
        y2 = fp_mul(y2, v)
        y2 = fp_add(y2, b)
        if fp_jacobi(y2):
            return encode(x,y)
        
        v = x
        # Compute coordinate y^2 in S
        y1 = fp_add(v, u)
        y1 = fp_mul(y1, v)
        y2 = fp_sqr(u)
        y1 = fp_add(y1, y2)
        y1 = fp_add(y1, a)
        y2 = fp_add(y2, a)
        y2 = fp_mul(y2, u)
        y2 = fp_add(y2, b)
        y2 = fp_neg(y2)
        y1 = fp_inv(y1)
        y2 = fp_mul(y2, y1)
        if not fp_jacobi(y2):
            return encode(x,y)

    # When x is x2    
    elif case == 2:
        v = fp_add(u, x)
        v = fp_neg(v)

        # Check that x1 did not yield a square
        y1 = fp_sqr(v)
        y1 = fp_add(y1, a)
        y1 = fp_mul(y1, v)
        y1 = fp_add(y1, b)
        if fp_jacobi(y1):
            return encode(x,y)

        # Compute coordinate y^2 in S
        y1 = fp_add(v, u)
        y1 = fp_mul(y1, v)
        y2 = fp_sqr(u)
        y1 = fp_add(y1, y2)
        y1 = fp_add(y1, a)
        y2 = fp_add(y2, a)
        y2 = fp_mul(y2, u)
        y2 = fp_add(y2, b)
        y2 = fp_neg(y2)
        y1 = fp_inv(y1)
        y2 = fp_mul(y2, y1)
        if not fp_jacobi(y2):
            return encode(x,y)


    # When x is x3 
    else:
        y2 = fp_sub(x,u)    #y2 = x-u
        if not fp_jacobi(y2):
            return encode(x,y)
        y1 = fp_sqr(u)      #y1 = u^2
        v = fp_mul(y1, y2)  #v = u^2*y2
        y1 = fp_add(y1, a)  #y1 = u^2+a
        y1 = fp_add(y1, y1)
        y1 = fp_add(y1, y1) #y1 = 4(u^2+a)
        y1 = fp_mul(y1, y2) #y1 = 4y2(u^2+a)
        v = fp_sub(y1, v)   #v = y2(4a+3u^2)
        v = fp_mul(v, y2)   #v = y2^2(4a+3u^2)
        y1 = fp_mul(y1, u)
        v = fp_add(v, y1)   #v = y2^2(4a+3u^2) + 4y2(u^3+au)
        y1 = fp_mul(b, y2)
        y1 = fp_add(y1, y1)
        y1 = fp_add(y1, y1)
        v = fp_add(v, y1)   # v = y2^2(4a+3u^2) + 4y2(u^3+au+b)
        v = fp_neg(v)
        if fp_jacobi(v):
            v = fp_sqrt(v)
        else:
            return encode(x,y)
        if case == 4:
            v = fp_neg(v)
        y1 = fp_mul(u, y2)
        v = fp_sub(v, y1)
        w = fp_add(y2, y2)
        

    # Compute conic coordinates

    Y1 = fp_sqrt(y2)

    if case > 2:
        X1 = fp_add(v, v)
        Z1 = fp_mul(u, w)
        X1 = fp_add(X1, Z1)
        X1 = fp_mul(X1, Y1)
        Y1 = fp_mul(Y1, w)
        Z1 = fp_add(w, w)
        
        t = fp_mul(Y1, u)
        t = fp_mul(t, X[0])   #X[0] is just sqrt(-3)
        t = fp_sub(t, X1)
        Z1 = fp_inv(Z1)
        t = fp_mul(t, Z1)
    
    else:
        X1 = fp_add(v, v)
        X1 = fp_add(X1, u)
        X1 = fp_mul(X1, Y1)

        t = fp_mul(Y1, u)
        t = fp_mul(t, X[0])   #X[0] is just sqrt(-3)
        t = fp_sub(t, X1)
        t = fp_mul(t, X[1])   #X[1] is just 2^-1 - this saves us from doing any inverses in cases 1,2

    s = int(y) % 2

    return u,t,s