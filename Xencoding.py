from encoding import *
from random import random

def Xdecode(u, t):

    #Evalaute initial point in conic X(u),Y(u)
    X0 = fp_mul(X[2], u)
    X0 = fp_add(X0, X[1])
    X0 = fp_mul(X0, u)
    X0 = fp_add(X0, X[0])
    Y0 = fp_mul(Y[1], u)
    Y0 = fp_add(Y0, Y[0])

    #print("Test X0,Y0,Z0",X0**2+(3*u**2+4*a)*Y0**2+(u**3+a*u+b)==0)

    # Evaluate f(u)=3u^2+4a and g(u)=u^3+au+b
    f = fp_sqr(u)
    g = fp_add(f, a)
    g = fp_mul(g, u)
    g = fp_add(g, b)
    Z1 = fp_add(f, f)
    f = fp_add(Z1, f)
    f = fp_add(f, ax4)

    #print("Test f,g",X0**2+f*Y0**2+g==0)

    # Compute new point in conic
    #   X1 = f*(Y0-t*X0)^2 + g
    #   Z1 = X0(1 + f*t^2)
    #   Y1 = Z1*Y0 + t*(X - Z*X0)
    Z1 = fp_mul(t, X0)
    Y1 = fp_sub(Y0, Z1)
    X1 = fp_sqr(Y1)
    X1 = fp_mul(X1, f)
    X1 = fp_add(X1, g)
    Z1 = fp_mul(Z1, t)
    Z1 = fp_mul(Z1, f)
    Z1 = fp_add(Z1, X0)
    Y1 = fp_mul(Y1, Z1)
    tX = fp_mul(t, X1)
    Y1 = fp_add(Y1, tX)
    #assert(not Z1 == 0)

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

    # Compute projective point in V
    #   x1 = v*y*w^6
    #   x2 = -y*(w*u+v)*w^6
    #   x3 = (w^2*u+y^2)*y*w^5
    #   YY  = ( (w^2*u + y^2)^3 + a*w^4*(w^2*u + y^2) + w^6*b ) * ( w^2*u + u*v*w + v^2 + w^2*a )
    #   ZZ  = y*w^7
    ZZ = fp_sqr(w)
    x2 = fp_mul(u, w)
    x1 = fp_mul(x2, w)
    x2 = fp_add(x2, v)
    x2 = fp_mul(x2, w)
    x2 = fp_neg(x2)
    x3 = fp_sqr(y)
    x3 = fp_add(x3, x1)
    x1 = fp_mul(v, w)

    # Compute ZZ*xi^3+A*ZZ^3*xi+B*ZZ^4
    w = fp_sqr(ZZ)
    t2 = fp_sqr(w)
    t1 = fp_mul(w, a)
    t2 = fp_mul(t2, b)
    
    y23 = fp_sqr(x3)
    y23 = fp_add(y23, t1)
    y23 = fp_mul(y23, ZZ)
    y23 = fp_mul(y23, x3)
    y23 = fp_add(y23, t2)
    
    y22 = fp_sqr(x2)
    y22 = fp_add(y22, t1)
    y22 = fp_mul(y22, ZZ)
    y22 = fp_mul(y22, x2)
    y22 = fp_add(y22, t2)

    # Find the square
    c2 = fp_jacobi(y22)
    c3 = fp_jacobi(y23)
    x1, x2 = fp_cswap(c2, x1, x2)
    x1, x3 = fp_cswap(c3, x1, x3)
    
    return x1, ZZ


def Xencode(x, z):
    u = F.random_element()
    case = ceil(4*random())
    
    # When x is x1
    if case == 1:
        
        # Check that x_2 doesnt yield a square
        v = fp_neg(u)
        v = fp_mul(v, z)
        v = fp_sub(v, x)
        y1 = fp_sqr(v)
        y2 = fp_sqr(z)
        Y1 = fp_mul(y2, a)     #w^2a
        Y1 = fp_add(Y1, y1)    #w^2a+v^2
        Y1 = fp_mul(Y1, v)     #w^2va+v^3
        Y1 = fp_mul(Y1, z)     #w^3va + wv^3
        y2 = fp_sqr(y2)
        y2 = fp_mul(y2, b)     #w^4b
        Y1 = fp_add(Y1, y2)    #wv^3 + w^3va + w^4b
        if fp_jacobi(Y1):
            return Xencode(x, z)

        v = x
        w = z
        # Compute coordinate y^2 in S
        y2 = fp_mul(u, w)
        y1 = fp_add(y2, v)
        y1 = fp_mul(y1, y2)
        y2 = fp_sqr(v)
        y1 = fp_add(y2, y1)
        Y1 = fp_sqr(w)
        y2 = fp_mul(Y1, a)
        y1 = fp_add(y1, y2)
        y2 = fp_sqr(u)
        y2 = fp_add(y2, a)
        y2 = fp_mul(y2, u)
        y2 = fp_add(y2, b)
        y2 = fp_mul(y2, Y1)
        y2 = fp_neg(y2)
        y1 = fp_inv(y1)
        y2 = fp_mul(y2, y1)
        if not fp_jacobi(y2):
            return Xencode(x, z)

    # When x is x2    
    elif case == 2:


        v = fp_mul(z, u)
        v = fp_add(v, x)
        v = fp_neg(v)
        w = z

        # Check that x1 did not yield a square
        y1 = fp_sqr(v)
        y2 = fp_sqr(w)
        Y1 = fp_mul(y2, a)     #w^2a
        Y1 = fp_add(Y1, y1)    #w^2a+v^2
        Y1 = fp_mul(Y1, v)     #w^2va+v^3
        Y1 = fp_mul(Y1, w)     #w^3va + wv^3
        y2 = fp_sqr(y2)
        y2 = fp_mul(y2, b)     #w^4b
        Y1 = fp_add(Y1, y2)    #wv^3 + w^3va + w^4b
        if fp_jacobi(Y1):
            return Xencode(x, z)

        # Compute coordinate y^2 in S
        y2 = fp_mul(u, w)
        y1 = fp_add(y2, v)
        y1 = fp_mul(y1, y2)
        y2 = fp_sqr(v)
        y1 = fp_add(y2, y1)
        Y1 = fp_sqr(w)
        y2 = fp_mul(Y1, a)
        y1 = fp_add(y1, y2)
        y2 = fp_sqr(u)
        y2 = fp_add(y2, a)
        y2 = fp_mul(y2, u)
        y2 = fp_add(y2, b)
        y2 = fp_mul(y2, Y1)
        y2 = fp_neg(y2)
        y1 = fp_inv(y1)
        y2 = fp_mul(y2, y1)
        if not fp_jacobi(y2):
            return Xencode(x, z)


    # When x is x3 
    else:
        X1 = fp_inv(z)
        X1 = fp_mul(X1, x)
        y2 = fp_sub(X1,u)    #y2 = x-u
        if not fp_jacobi(y2):
            return Xencode(x, z)
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
            return Xencode(x,z)
        if case == 4:
            v = fp_neg(v)
        y1 = fp_mul(u, y2)
        v = fp_sub(v, y1)
        w = fp_add(y2, y2)


    # Compute points in conic
    Y1 = fp_sqrt(y2)

    X0 = fp_mul(X[2], u)
    X0 = fp_add(X0, X[1])
    X0 = fp_mul(X0, u)
    X0 = fp_add(X0, X[0])
    Y0 = fp_mul(Y[1], u)
    Y0 = fp_add(Y0, Y[0])
    X0 = fp_add(X0, X0)
    Y0 = fp_add(Y0, Y0)

    X1 = fp_add(v, v)
    y2 = fp_mul(u, w)
    X1 = fp_add(X1, y2)
    X1 = fp_mul(X1, Y1)
    Y1 = fp_mul(Y1, w)
    X0 = fp_mul(X0, w)
    Y0 = fp_mul(Y0, w)

    t = fp_sub(Y1, Y0)
    y1 = fp_sub(X1, X0)
    y1 = fp_inv(y1)
    t = fp_mul(t, y1)

    return u,t