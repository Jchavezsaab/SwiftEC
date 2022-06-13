from encoding import *
import matplotlib.pyplot as plt
import numpy as np

def main():
    counts = {}
    inf = -1
    counts[inf] = 0
    for i in range(p):
        x = F(i)
        y2 = x**3+a*x+b
        if y2.is_square():
            counts[i] = 0

    print("Progress: 0%",end="\r")
    for i in range(p):
        print("\033[KProgress: "+str(math.floor(i/p*100))+"%",end="\r")
        t = F(i)
        for j in range(p):
            u = F(j)
            try:
                x,y = decode(t, u, 0)
                counts[int(x)] += 1
            except PointAtInfinity:
                counts[inf] += 1

    count_mean = float(mean(counts.values()))
    count_min = min(counts.values())
    count_max = max(counts.values())
    count_std = float(std(counts.values()))

    y_pos = np.arange(len(counts.values()))
    plt.bar(y_pos, counts.values(), 1, align='center', alpha=0.5)
    plt.ylabel('Counts')
    plt.xlabel('X coord')
    plt.title("Number of preimages for "+curve)
    
    plt.annotate("Mean counts: "+str(count_mean), [.8*y_pos[-1],.9*count_max])
    plt.annotate("Min counts: "+str(count_min), [.8*y_pos[-1],.8*count_max])
    plt.annotate("Max counts: "+str(count_max), [.8*y_pos[-1],.7*count_max])
    plt.annotate("STD: "+str(count_std)+" ("+str(math.ceil(count_std/count_mean*100))+"%)", [.8*y_pos[-1],.6*count_max])

    plt.show()

    print("MEAN: "+str(count_mean))
    print("MIN: "+str(count_min))
    print("MAX: "+str(count_max))
    print("STD: "+str(count_std)+" ("+str(math.ceil(count_std/count_mean*100))+"%)")

if __name__ == "__main__":
    main()

