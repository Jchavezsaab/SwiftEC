from Xencoding_0 import *

def main():

    if a != 0:
        print("Error: for curves with a!=0, call Xenc_test.py instead.")
        exit(-1)

    E = EllipticCurve([0, 0, 0, a, b])
    x,y = E.random_point().xy()
    z = F.random_element()
    x = x*z

    print("+++PARAMETERS+++")
    print("CURVE:\t"+sys.argv[2]) 
    print("a=\t"+str(hex(int(a))))
    print("b=\t"+str(hex(int(b))))
    print("\n+++INITIAL POINT+++")
    print("x=\t"+str(hex(int(x/z))))

    init_counters()
    u,t = Xencode(x,z)
    print("\n+++ENCODING+++")
    print("u=\t"+str(hex(int(u))))
    print("t=\t"+str(hex(int(t))))
    print("Encoding Cost:")
    print_counters()

    init_counters()
    x1,z1 = Xdecode(u,t)
    print("\n+++DECODING+++")
    print("x=\t"+str(hex(int(x1/z1))))
    print("Decoding Cost:")
    print_counters()
    
    assert(x/z == x1/z1)

if __name__ == "__main__":
    main()

