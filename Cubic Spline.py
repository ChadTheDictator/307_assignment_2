# Author = Nick Mrazek
'''Create a set of cubic spline functions based on given dataset'''

import csv
import matplotlib.pyplot as plt

file  = "SV&SH_data.csv"

# SPP gaussian elimination function created by Garrett Posehn
def get_pivot(
    matrix:list[list[float|int]],
    index:int,
    pivoted:list[int]
  ) -> int:
  '''
    Gets the index of the pivot row for gauss elim. with partial pivoting.

    Args:
      matrix: The matrix of size n*n being solved
      index: The current column index being worked on
      pivoted: A list of the already pivoted row indices

    Returns:
      pivot_row: The index of the pivot row
  '''

  # Defining scale vector
  lst = [max(abs(x) for x in row) for row in matrix]
  pivot = []

  # Determining pivot row
  for i in range(len(matrix)):
    if i in pivoted:
      pivot.append(float('-inf'))
    else:
      pivot.append((abs(matrix[i][index]))/(abs(lst[i])))
  pivot_row = pivot.index(max(pivot))
  return(pivot_row)

def back_sub(
    matrix:list[list[int|float]],
    outs:list[int|float],
    order:list[int]|None=None
    )-> list[int|float]:
  """
  Performs backwards substitution on any n*n matrix.
  @Args:,
    matrix: The matrix being solved
    outs: A list of expected outputs
    order(optional): The order in which the matrix has been pivoted.
      assuming it was not re-arranged to be lower triangular

  @Returns:
    returns result of the back substitution in order:
      [x1,x2,x3,...,xn]
  """
  if order == None:
    order = list(range(len(matrix)-1,-1,-1))

  var_locs = order[::-1]

  for i in order:
    cd_index = (len(order)-(order.index(i))-1) # Center Diagonal Index

    # Back substituting anything not in a CD spot
    total = 0
    for n in range(cd_index+1,len(order)):
      outs[i] -= (matrix[i][n]*outs[var_locs[n]])

    # Normalizing
    if matrix[i][cd_index] != 0:
      outs[i] = outs[i]/matrix[i][cd_index]
      matrix[i] = list(map(lambda x: x/matrix[i][cd_index],matrix[i]))

  x = []
  for n in var_locs:
    x.append(outs[n])
  return(x)

def gauss_elim_p_piv(
    matrix:list[list[float|int]],
    var:list[str],
    outs:list[int|float]
  ):
  '''
    Performs gaussian elimination with partial pivoting on provided
    matrix and outputs using provided variable names.

    Args:
      matrix: The matrix of size n*n being solved
      var: A list of variable names
      outs: A list of anticipated outputs for each variable

    returns:
      answer: A dictionary of the answers where the variable names are
        mapped to their values.

  '''

  order_pivot = []
  rows = len(matrix)
  columns = len(matrix[0])
  if rows == columns:

    # Partial Pivoting
    for column_index in range(columns-1):
      # Getting Pivot index
      pivot_index = get_pivot(matrix,column_index,order_pivot)
      order_pivot.append(pivot_index)

      # Performing Row Order Operations
      for row_index in range(columns):

        # if the row index isn't on main diagonal, and has not been pivoted already
        if row_index != pivot_index and row_index not in order_pivot:

          scaler = matrix[row_index][column_index]/matrix[pivot_index][column_index]
          for i in range(len(matrix[row_index])):
            matrix[row_index][i] = matrix[row_index][i] - (scaler*matrix[pivot_index][i])
            if abs(matrix[row_index][i]) < (1*(10**-6)):
              matrix[row_index][i] = 0
          outs[row_index] = outs[row_index]-(outs[pivot_index]*scaler)
          
        if column_index == columns-2 and row_index not in order_pivot:
          order_pivot.append(row_index)

  # Reverse Substitution
  reverse_order = order_pivot[::-1]
  outs = back_sub(matrix,outs,reverse_order)
  answer = {}

  for i in range(len(outs)):
    answer[var[i]] = (outs[i])

  # Uncomment below to print out pivot orders
  # print(f'Pivot Order: {order_pivot}')
  return(answer,(reverse_order[::-1]))


# ----- Main Program Functions -----

def read_csv(file):
    with open(file,'r') as file:    #read csv, retrieve sv/sh data. 
        reader = csv.reader(file)
        sv = []
        sh = []
        for SV,SH in reader: 
            sv.append(float(SV))
            sh.append(float(SH))
    return sv, sh

def create_S(x_value, zlist, sv, sh, h):
    #Use spline function for a given interval to estimate a value of sh given a value for sv.

    interval_num = 0
    n = len(sv) - 1

    for i in range(n):          # Find correct interval
        if sv[i] <= float(x_value) <= sv[i+1]:
            interval_num = i
            break

    i = interval_num
    hi = h[i]    

    # Use spline function in expanded form to estimate. 
    term_1 = zlist[i] * (sv[i+1] - x_value)**3 / (6 * hi)
    term_2 = zlist[i+1] * (x_value - sv[i])**3 / (6 * hi)
    term_3 = (sh[i] - zlist[i]*hi**2/6) * (sv[i+1] - x_value) / hi
    term_4 = (sh[i+1] - zlist[i+1]*hi**2/6) * (x_value - sv[i]) / hi
        
    return (term_1 + term_2 + term_3 + term_4)

def get_spline_matrix(file):
    # Creates a system of equations to create spline equations
   

    sv, sh = read_csv(file)
    n = len(sv)-1

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

        # Simplified into ax3+bx2+cx+d form.
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

    # Print Data in table
    '''
    print("SV\t|\tSH")
    for i in range(len(x_plot)):
        print(f"{round(x_plot[i], 3)}\t{round(y_plot[i], 5)}")
    '''
    # plot spline using "y" and "x" values
    plt.plot(x_plot, y_plot)

    # plot original data
    plt.scatter(sv, sh)

    plt.xlabel("Specific Volume (m^3/kg)")
    plt.ylabel("Specific Enthalpy (kJ/kg)")
    plt.title("Specific Enthalpy vs. Specific Volume")
    plt.grid()

    plt.show()

#------------Main Code----------

sv, sh = read_csv(file)
n = len(sv)-1

A,v,h = get_spline_matrix(file)
#print(A)
#print(v)

vars = []       # create list of variable names
for i in range(n-1):
    vars+=["z"+str(i+1)]

#print(vars)
zlist = []  # Zlist = list of second derivatives.
z, nums = gauss_elim_p_piv(A, vars, v) # A = matrix, vars = variable names, v = solution vector

for i in range(len(nums)):  #convert into list instead of dict.
    zlist += [z.get(vars[i])]
zlist =[0] + zlist + [0] 



# OUTPUTS

# Estimate sh at a given sv
x_value = float(input("value of sv to estimate sh at:"))
givenpt = create_S(x_value, zlist, sv, sh, h)   
print("Spline value at", x_value, "=", givenpt)

# Print out spline functions
print("\n\nSPLINE FUNCTIONS")
a, b, c, d = print_spline_fxs(sv, sh, zlist, h)

# Plot splines
plot_splines(sv, sh, a, b, c, d)