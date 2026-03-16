# Author Nick M.
# transient analysis via 
from numpy import log
import matplotlib.pyplot as plt

#function = integral( ln(r)/sqrt(r))dr from 0->1

# using substitution r=x^2, dr=2xdx we get integral( ln(x^2)/sqrt(x^2))2xdx
# simplifying---> integral( 2ln(x)/x)2xdx = integral(4ln(x))dx

# bounds of integration do not change, 1^2 = 1, 0^2 = 0

a = 10**-12
b = 1 

n = int(input("Number of sub-intervals(multiple of 3):"))
h = (b-a)/n


def f(r):
    height = 4*log(r) #using modified function
    return height

def simp38(n):
    h = (b-a)/n
    area = 0

    for i in range(0, n, 3): # n/3 subdivision groupings:
        x0 = a + i*h    #start at 1st point in grouping of 3, i*h being the grouping # 
        x1 = x0 + h     #next pts based on the first of the interval
        x2 = x0 + 2*h
        x3 = x0 + 3*h
        area+= f(x0) + 3*f(x1) + 3*f(x2) + 1*f(x3)
    

    area*= 3*h/8 #Multiply by the 

        
    return area

def plot_integrals(aprxs, nlist):
    plt.figure()
    plt.plot(nlist, aprxs, marker='*')
    plt.xlabel("n Value")
    plt.ylabel("Approximation Value")
    plt.title("Approximation Value vs n Value")
    plt.show()

approx = simp38(n)
print(f'integral = {approx}')

#justify b)


aprxs = []
nlist = []
n=0
#C) Graphing every 

for i in range(1,300,10):
    n = 3*i
    aprxs += [simp38(n)]
    nlist += [n]

#print(aprxs)
plot_integrals(aprxs, nlist)

