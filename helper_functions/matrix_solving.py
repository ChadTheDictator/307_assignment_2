# Author: Garrett Posehn

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


# ----- Main Program Functions defined below -----

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
  return(answer)

