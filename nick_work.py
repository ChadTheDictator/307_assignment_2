# Author Nick M.

from math import log
import matplotlib.pyplot as plt

#function = integral( ln(r)/sqrt(r))dr from 0->1

# using substitution r=x^2, dr=2xdx we get integral( ln(x^2)/sqrt(x^2))2xdx
# simplifying---> integral( 2ln(x)/x)2xdx = integral(4ln(x))dx

# bounds of integration do not change, 1^2 = 1, 0^2 = 0

a = 10**-15
b = 1 

# Node/weight data for 5 pt quadrature
nodes = [0, -0.538469, 0.538469, -0.906180, 0.906180]
weights = [0.568889, 0.478629, 0.478629, 0.236927, 0.236927]

    
 
def f(r):
    height = 4*log(r) #using modified function
    return height

def gauss5pt(n):
 
    area = 0 
    h = (b-a)/n

    for i in range(n): #for each subinterval
        left_end = a +i*h
        right_end = left_end+h # right boundary, a + width of subinterval
        midpt = (left_end+right_end)/2 # finds location of middle point of subinterval
        half = (right_end-left_end)/2  # halflength of subinterval


        for j in range(5):      
            x = midpt + half*nodes[j]
            area +=weights[j]*f(x)*half

    return area


# Main code

n = int(input("Number of sub-intervals:"))
approx = gauss5pt(n)
print(f'integral = {gauss5pt(n)}')

approxes = []
nlist = []

print(f"n\t\t|\tApproximation value")
for i in range(1,21):
    approxes+= [gauss5pt(i)]
    nlist += [i]

for i in range(len(approxes)):
    print(f"{nlist[i]}\t\t|\t {round(approxes[i], 5)}")

