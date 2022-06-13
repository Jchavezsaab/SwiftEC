from isogenies import *

def main():

    P = E1.random_point()
    x,y = P.xy()

    print("+++PARAMETERS+++")
    print("CURVE:\t"+sys.argv[2]) 
    print("a=\t"+str(hex(int(a))))
    print("b=\t"+str(hex(int(b))))
    print("\n+++INITIAL POINT+++")
    print("x=\t"+str(hex(int(x))))
    print("y=\t"+str(hex(int(y))))

    init_counters()
    u,t,s = encisog(x,y)
    print("\n+++ENCODING+++")
    print("u=\t"+str(hex(int(u))))
    print("t=\t"+str(hex(int(t))))
    print("s=\t"+str(s))
    print("Encoding Cost:")
    print_counters()

    init_counters()
    x1,y1 = decisog(u,t,s)
    print("\n+++DECODING+++")
    print("x=\t"+str(hex(int(x1))))
    print("y=\t"+str(hex(int(y1))))
    print("Decoding Cost:")
    print_counters()
    
    assert((index*P).xy()[0] == x1)
    assert((index*P).xy()[1] in [y1, -y1])

if __name__ == "__main__":
    main()

