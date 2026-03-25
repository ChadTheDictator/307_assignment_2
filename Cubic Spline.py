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

def create_S(x_value, zlist, sv, sh, h):
    interval_num = 0
    n = len(sv) - 1

    for i in range(n):
        if sv[i] <= float(x_value) <= sv[i+1]:
            interval_num = i
            break

    i = interval_num
    hi = h[i]    
        
    term_1 = zlist[i] * (sv[i+1] - x_value)**3 / (6 * hi)
    term_2 = zlist[i+1] * (x_value - sv[i])**3 / (6 * hi)
    term_3 = (sh[i] - zlist[i]*hi**2/6) * (sv[i+1] - x_value) / hi
    term_4 = (sh[i+1] - zlist[i+1]*hi**2/6) * (x_value - sv[i]) / hi
        
    return (term_1 + term_2 + term_3 + term_4)

def get_spline_matrix(file):
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
        A+=[row]


    return A, v, h 
#------------Main Code----------

x_value = float(input("value of sv to estimate sh at:"))
sv, sh = read_csv(file)
n = len(sv)-1
    
A,v,h = get_spline_matrix(file)
#print(A)
#print(v)

vars = []
for i in range(n-1):
    vars+=["z"+str(i+1)]

#print(vars)
zlist = []
z, nums = gauss_elim_p_piv(A, vars, v) #via helper function for partial pivoting gaussian elim.

for i in range(len(nums)):  #convert into list instead of dict.
    zlist += [z.get(vars[i])]
zlist =[0] + zlist + [0] 

a = create_S(x_value, zlist, sv, sh, h)

print("Spline value at", x_value, "=", a)

