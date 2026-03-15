# Author Nick M.
# transient analysis via 
from numpy import log
#function = ln(r)/r
a = 0
b = 1 
n = int(input("Number of sub-intervals(multiple of 3):"))

h = (b-a)/n


def simp38(n):
    area = 0
    heights = []
    for i in range(n):
        if  i*h <= 10^-12:
            heights += log((10^-12))/(10^-12)
        else:
            heights += log((i*h))/(i*h)

    area+=heights[0]        #add endpoints
    area+=heights[-1]

    for i in range(1,n-1):
        area+=3*heights[i]
    area*=3/8


        
    return
