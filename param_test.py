from generate_parameters import *

def print_polys(polys):
	print("X=",polys[0][0],"+",polys[0][1],"* t +",polys[0][2],"* t^2")
	print("Y=",polys[1][0],"+",polys[1][1],"* t")
	print("Z=",polys[2][0])

def in_conic(point, a, b, t):
	X = point[0]
	Y = point[1]
	Z = point[2]
	return X**2 + (3*t**2+4*a)*Y**2 + (t**3 + a*t + b)*Z**2 == 0

def poly_eval(polys, u):
	"""Evaluates X(u),Y(u),Z(u)"""
	point = []
	for i in range(len(polys)):
		point.append(0)
		for d in range(len(polys[i])):
			point[i] += polys[i][d]*u**d
	return point

def in_S(point, a, b, t):
	v = point[0]
	y = point[1]
	w = point[2]
	return y**2*(t**2*w**2 + v*t*w + v**2 + a*w**2) == -w**4*(t**3+a*t+b)

def poly_multieval(polys, t, u):
	"""Evaluates alpha(t,u),beta(t,u),gamma(t,u)"""
	point = []
	for i in range(len(polys)):
		point.append(0)
		for d1 in range(len(polys[i])):
			for d2 in range(len(polys[i][d1])):
				point[i] += polys[i][d1][d2]*t**d1*u**d2
	return point

def main(iterations):

	print("\033[KProgress: 0%",end="\r")
	for iteration in range(iterations):
		print("\033[KProgress: "+str(math.floor(iteration/iterations*100))+"%",end="\r")
		# Any choice of signs will work and produce a different parametrization
		signs = [1-2*int(round(random())),1-2*int(round(random())),1-2*int(round(random())),1-2*int(round(random())),1-2*int(round(random()))]
		X0u,Y0u,Z0u,case = conic_point(a,b,p,signs)
		X1tu,Y1tu,Z1tu = conic_parametrization(X0u,Y0u,Z0u, a, b)
		u = F.random_element()
		t = F.random_element()
		X0,Y0,Z0 = poly_eval((X0u,Y0u,Z0u), u)
		X1,Y1,Z1 = poly_multieval((X1tu,Y1tu,Z1tu), t, u)
		assert(in_conic((X0,Y0,Z0), a, b, u))
		#assert(in_conic((X1,Y1,Z1), a, b, u))
	print("All tests passed")


if __name__ == "__main__":
    main(1024)


#Case 1
#a=Mod(19338309632770419086156185150203927566229393627978547513184166475234416919874,P256)
#b=Mod(70666588922232729606665262871975404666340305607818632618971496586132780213941,P256)
#Case 2
#a=Mod(49002211989758983106857941188968855211451296442363798735770804588539430089297,P256)
#b=Mod(67481648871508048090344741784856540682703176942379260294063979557761568307167,P256)
#Case 3
#a=Mod(32651201986284816203485320389763519797417674197317260843268792113488522340136,P256)
#b=Mod(52929209091472368122750625035971472539267788555382426132109478834668639862094,P256)
#Case 4
#a=Mod(27604794116331584293826600666414899828856610780898312398942711365269289829480,P256)
#b=Mod(52439542961005085591279224534662415077119473180968504672412193213438792230303,P256)
#Case Random
#a = F.random_element()
#b = F.random_element()
#print(a,b)