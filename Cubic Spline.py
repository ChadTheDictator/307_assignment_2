import csv
import matplotlib.pyplot as plt

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
    '''
    gets system of equations to create spline equations
    '''
    
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

def print_spline_fxs(sv, sh, zlist, h):
    n = len(sv)-1

    print("\nCubic Spline Functions (Expanded Form):\n")
    
    a = []  # Coefficients ax3 +bx2 +cx +d
    b = []
    c = []
    d = []

    for i in range(n):
        xi = sv[i]
        xi1 = sv[i+1]
        hi = h[i]

        zi = zlist[i]
        zi1 = zlist[i+1]

        # Coefficients from standard spline form
        A = zi/(6*hi)
        B = zi1/(6*hi)
        C = (sh[i] - zi * hi**2 / 6)/hi
        D = (sh[i+1] - zi1 * hi**2 / 6)/hi

        # (xi1 - x)^3 = -x^3 + 3xi1 x^2 - 3xi1^2 x + xi1^3
        # (x - xi)^3 = x^3 - 3xi x^2 + 3xi^2 x - xi^3
        
        ai = -A + B
        bi = 3*A*xi1 - 3*B*xi
        ci = -3*A*(xi1**2) + 3*B*(xi**2) - C + D
        di = A*(xi1**3) - B*(xi**3) + C*xi1 - D*xi

        print(f"Interval [{xi}, {xi1}]:")
        print(f"S_{i}(x) = ({ai:.4e})X^3 + ({bi:.4e})X^2 + ({ci:.4e})X + ({di:.4e})\n")
        
        a += [ai]   # Add to list for later use
        b += [bi]
        c += [ci]
        d += [di]
    return a, b, c, d

def plot_splines(sv, sh, a, b, c, d):
    x_plot = []
    y_plot = []

    n = len(sv) - 1

    for j in range(n):
        xi = sv[j]
        xi1 = sv[j+1]

        num_points = 100
        step = (xi1 - xi) / num_points

        x = xi

        for k in range(num_points + 1):
            y = a[j]*x**3 + b[j]*x**2 + c[j]*x + d[j]

            x_plot.append(x)
            y_plot.append(y)

            x += step

    # Table for data
    '''
    print("SV\t|\tSH")
    for i in range(len(x_plot)):
        print(f"{round(x_plot[i], 3)}\t{round(y_plot[i], 5)}")
    '''
    # plot spline
    plt.plot(x_plot, y_plot)

    # plot original data
    plt.scatter(sv, sh)

    plt.xlabel("Specific Volume (m^3/kg)")
    plt.ylabel("Specific Enthalpy (kJ/kg)")
    plt.title("Specific Enthalpy vs. Specific Volume")
    plt.grid()

    plt.show()


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

givenpt = create_S(x_value, zlist, sv, sh, h)

print("Spline value at", x_value, "=", givenpt)

print("\n\nSPLINE FUNCTIONS")
a, b, c, d = print_spline_fxs(sv, sh, zlist, h)

plot_splines(sv, sh, a, b, c, d)


