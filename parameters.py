from config import *

X = []
Y = []

try:

    f = open('./parameters/' + curve)
    X.append(F(f.readline()))
    X.append(F(f.readline()))
    X.append(F(f.readline()))
    Y.append(F(f.readline()))
    Y.append(F(f.readline()))
    ax4 = F(int(f.readline()))
    pm1d2 = F(int(f.readline()))
    pbits = int(f.readline())
    f.close()

except IOError:

    raise MissingFileError("Precomputation not found, try running \"python generate_parameters.py -p "+curve+"\" first.")
    