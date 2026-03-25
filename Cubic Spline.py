import csv
from helper_functions.matrix_solving import gauss_elim_p_piv

file  = "SV&SH_data.csv"

def read_csv(file):
    with open(file,'r') as file:    #
        reader = csv.reader(file)
        sv = []
        sh = []
        for SV,SH in reader: 
            sv.append(float(SV))
            sh.append(float(SH))
    return sv, sh

def S(x_value, x, h):
    interval_num = 0
    for i in range(n):
        if x[i] <= x_value <= x[i+1]:
            interval_num = i
            break
    
    hi = h[i]

    #term_1 = 
    return

def spline(file):
    sv, sh = read_csv(file)
    n = len(sv)-1

    z=[]

    h = []
    u = []
    v = []
    b = []
    y = sh
    x = sv
    A = []
    for i in range(n):
        h += [(x[i+1] - x[i])]
    
    for i in range(n):
        b += [(y[i+1]-y[i])/h[i]]
    
    for i in range(1,n):
        u += [2*(h[i-1] + h[i])]
        v += [6*(b[i]-b[i-1])]

    for i in range(n-1):
        row = [0]*(n-1)

        # Main diagonal
        row[i] = u[i]

        # Lower diagonal
        if i > 0:
            row[i-1] = h[i]

        # Upper diagonal
        if i < n-2:
            row[i+1] = h[i+1]

        A.append(row)

       
        v.append(6*(b[i+1] - b[i]))

    return A, v
    #print(h,u,b,v,y,x)
   # v+=[h[i-1]*z[i-1]+u[i]*z[i]+h[i]*z[i+1]]
sv, sh = read_csv(file)
n = len(sv)-1
    
A,v = spline(file)
#vars = []
#for i in range(n):
#    vars+=["z"]

z = gauss_elim_p_piv(A,v,vars)

print(z)


#print(sv)
#print(sh)
