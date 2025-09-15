import math
import sys
import numpy as np
import pandas as pd
import symnmf_c as mf

# Use the same RNG seed as the previous reference implementation to ensure
# deterministic, comparable initialization across versions.
np.random.seed(0)

#this function calculates and returns the average of all entries of the normalized matrix W,
#which has num_of_elements rows and num_of_elements columns.
def calculate_average(W, num_of_elements):
    sum = 0.0
    for i in range(num_of_elements):
        for j in range(num_of_elements):
            sum += W[i][j]
    average = sum/(num_of_elements**2)
    return average

#this function randomly initializes the matrix H for the symNMF algorythm 
# from the normalized matrix W. k is the number of required clusters.
# W has num_of_elements rows and num_of_elements columns.
# the returned matrix has num_of_elements rows and k columns.

def init_H0(W, num_of_elements, k):
    average = calculate_average(W, num_of_elements)
    return np.random.uniform(low=0.0, high=2*math.sqrt(average/k), size=(num_of_elements, k))

def initialize_H(W, num_of_elements, k):
    H0 = init_H0(W, num_of_elements, k)
    return H0

def is_float(n):
    try:
        float(n)
        return True
    except:
        return False

#this function prints the elements of the given matrix, which has num_of_rows rows 
#and num_of_cols columns, separated by commas    
def printMatrix(matrix, num_of_rows, num_of_cols):
    for i in range(num_of_rows):
        for j in range(num_of_cols):
            print('%.4f' % matrix[i][j], end='')
            if j<num_of_cols-1:
                print(",", end='')
        print('')
    return

def main():
    k = 0  # default; only meaningful for the symnmf goal
    file_name = ""
    elements = []
    num_of_elements = 0
    d = 0

    # Tester always invokes: python3 symnmf.py <k> <goal> <filename>
    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        exit(0)

    k_arg = sys.argv[1]
    goal = str(sys.argv[2])
    file_name = str(sys.argv[3])

    # Only parse & validate k if the goal actually requires it (symnmf)
    if goal == "symnmf":
        if not k_arg.isdigit():
            print("An Error Has Occurred")
            exit(0)
        k = int(k_arg)
    else:
        # For goals sym / ddg / norm ignore the k argument (tester passes 0)
        try:
            # Accept numeric (possibly 0); ignore value
            float(k_arg)
        except ValueError:
            print("An Error Has Occurred")
            exit(0)

    file = pd.read_csv(file_name, header=None)
    elements = file.values  # numpy ndarray of datapoints
    num_of_elements = len(elements)
    d = len(elements[0]) if num_of_elements > 0 else 0

    if num_of_elements == 0 or d == 0:
        print("An Error Has Occurred")
        exit(0)

    # Validate k only for symnmf goal
    if goal == "symnmf" and (k < 1 or k > num_of_elements):
        print("An Error Has Occurred")
        exit(0)

    if goal == "sym":
        printMatrix(mf.sym(elements.tolist(), num_of_elements, d), num_of_elements, num_of_elements)
    elif goal == "ddg":
        printMatrix(mf.ddg(elements.tolist(), num_of_elements, d), num_of_elements, num_of_elements)
    elif goal == "norm":
        printMatrix(mf.norm(elements.tolist(), num_of_elements, d), num_of_elements, num_of_elements)
    elif goal == "symnmf":
        W = mf.norm(elements.tolist(), num_of_elements, d)
        H = initialize_H(W, num_of_elements, k)
        printMatrix(mf.symnmf(H.tolist(), W, k, num_of_elements), num_of_elements, k)
    else:
        print("An Error Has Occurred")
        exit(0)

if __name__ == "__main__":
    main()