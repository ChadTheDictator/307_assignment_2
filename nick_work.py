# Author Nick M.
# transient analysis via 
from math import log
import matplotlib.pyplot as plt

#function = integral( ln(r)/sqrt(r))dr from 0->1

# using substitution r=x^2, dr=2xdx we get integral( ln(x^2)/sqrt(x^2))2xdx
# simplifying---> integral( 2ln(x)/x)2xdx = integral(4ln(x))dx

# bounds of integration do not change, 1^2 = 1, 0^2 = 0

a = 10**-15
b = 1 

n = int(input("Number of sub-intervals(multiple of 3):"))
h = (b-a)/n
 
gweights = [0.56888888, 0.4786286704993,0.4786286704993,0.236926885056189,0.236926885056189]
nodes = [0, -0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640]



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

def plot_integrals(approxes, nlist):
    plt.figure()
    plt.plot(nlist, approxes, marker='*')
    plt.xlabel("Value of i, n = 10^i")
    plt.ylabel("Approximation Value")
    plt.title("Approximation Value vs Value of n")
    plt.show()

approx = simp38(n)
print(f'integral = {simp38(n)}')

approxes = []
nlist = []

for i in range(1, 9):
    approxes+= [simp38(10**i)]
    nlist += [i]
plot_integrals(approxes, nlist)

print(approxes)