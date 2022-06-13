from parameters import *

add_counter = 0
mul_counter = 0
sqr_counter = 0
jac_counter = 0
inv_counter = 0
sqrt_counter = 0

def init_counters():
    global add_counter
    global mul_counter
    global sqr_counter
    global jac_counter
    global inv_counter
    global sqrt_counter
    add_counter = 0
    mul_counter = 0
    sqr_counter = 0
    jac_counter = 0
    inv_counter = 0
    sqrt_counter = 0
    return

def print_counters():
    global add_counter
    global mul_counter
    global sqr_counter
    global jac_counter
    global inv_counter
    global sqrt_counter
    print(str(add_counter)+"A + "+str(mul_counter)+"M + "+str(sqr_counter)+"S + "+str(jac_counter)+"J + "+str(inv_counter)+"I + "+str(sqrt_counter)+"R")

def fp_add(a, b):
    global add_counter
    add_counter += 1
    return a+b

def fp_sub(a, b):
    global add_counter
    add_counter += 1
    return a-b

def fp_neg(a):
    global add_counter
    add_counter += 1
    return -a

def fp_mul(a, b):
    global mul_counter
    mul_counter += 1
    return a*b

def fp_sqr(a):
    global sqr_counter
    sqr_counter += 1
    return a**2

def fp_inv(a):
    global inv_counter
    inv_counter += 1
    return 1/a

def fp_sqrt(a):
    global sqrt_counter
    sqrt_counter += 1
    return a.nth_root(2)

def fp_jacobi(a):
    global jac_counter
    jac_counter += 1
    return a.is_square()

def fp_cswap(c, a, b):
    # Change this to constant time!
    if c:
        x = a
        a = b
        b = x
    return a,b

def fp_isodd(a):
    return int(a)%2 == 1
