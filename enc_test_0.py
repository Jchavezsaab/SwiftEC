from encoding_0 import *

def main():

    if a != 0:
        print("Error: for curves with a!=0, call enc_test.py instead.")
        exit(-1)

    E = EllipticCurve([0, 0, 0, a, b])
    x,y = E.random_point().xy()

    print("+++PARAMETERS+++")
    print("CURVE:\t"+sys.argv[2]) 
    print("a=\t"+str(hex(int(a))))
    print("b=\t"+str(hex(int(b))))
    print("\n+++INITIAL POINT+++")
    print("x=\t"+str(hex(int(x))))
    print("y=\t"+str(hex(int(y))))

    init_counters()
    u,t,s = encode(x,y)
    print("\n+++ENCODING+++")
    print("u=\t"+str(hex(int(u))))
    print("t=\t"+str(hex(int(t))))
    print("s=\t"+str(s))
    print("Encoding Cost:")
    print_counters()

    init_counters()
    x1,y1 = decode(u,t,s)
    print("\n+++DECODING+++")
    print("x=\t"+str(hex(int(x1))))
    print("y=\t"+str(hex(int(y1))))
    print("Decoding Cost:")
    print_counters()
    
    assert(x == x1)
    assert(y == y1)

if __name__ == "__main__":
    main()

