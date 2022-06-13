from isogenies import *

def main():

    P = E1.random_point()
    x,y = P.xy()
    z = F.random_element()
    x = x*z

    print("+++PARAMETERS+++")
    print("CURVE:\t"+sys.argv[2]) 
    print("a=\t"+str(hex(int(a))))
    print("b=\t"+str(hex(int(b))))
    print("\n+++INITIAL POINT+++")
    print("x=\t"+str(hex(int(x/z))))

    init_counters()
    u,t = Xencisog(x,z)
    print("\n+++ENCODING+++")
    print("u=\t"+str(hex(int(u))))
    print("t=\t"+str(hex(int(t))))
    print("Encoding Cost:")
    print_counters()

    init_counters()
    x1,z1 = Xdecisog(u,t)
    print("\n+++DECODING+++")
    print("x=\t"+str(hex(int(x1/z1))))
    print("Decoding Cost:")
    print_counters()
    
    assert((index*P).xy()[0] == x1/z1)

if __name__ == "__main__":
    main()

