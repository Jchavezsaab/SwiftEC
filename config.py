import sys
try:
    from sage.all import *
except:
    print("Failed to import sage module.")
    print("Try \"sage --python "+sys.argv[0]+" -p PRIME\"")
    print(" or \"python "+sys.argv[0]+" -p PRIME\"")
    exit()
    
class MissingFileError(Exception):
    pass

class PointAtInfinity(Exception):
    pass

if not (len(sys.argv) == 3 and sys.argv[1]=='-p'):
    print("USAGE: \"sage "+sys.argv[0]+" -p PRIME\"")
    exit()

curve = sys.argv[2]

try:

    f = open('./curves/' + curve)
    p = int(f.readline())
    F = GF(p)
    a = F(int(f.readline()))
    b = F(int(f.readline()))
    f.close()

except IOError:

    raise MissingFileError("Curve \""+curve+"\" not found, try adding it to curves/"+curve)
