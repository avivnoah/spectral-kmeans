import sys
import importlib.util
import os
import math

sys.path.insert(0, os.getcwd())
# load current built module
import symnmf_c as curr

# load previous built module from Prev_final_100 if available
prev_path = os.path.join(os.getcwd(), 'Prev_final_100')
prev = None
if os.path.isdir(prev_path):
    sys.path.insert(0, prev_path)
    for name in ('symnmfmodule', 'symnmf_c'):
        try:
            prev = __import__(name)
            break
        except Exception:
            prev = None

# read input
file_name = 'test_input_small4.txt'
import pandas as pd
file = pd.read_csv(file_name, header=None)
elements = file.values
num = len(elements)
d = len(elements[0])
K = 2

# compute W using current norm
W_curr = curr.norm(elements.tolist(), num, d)
import numpy as np
W_curr_arr = np.array(W_curr)
average = float(np.mean(W_curr_arr))

# create identical H using same init_H0 logic as symnmf.py
np.random.seed(1234)
H0 = np.random.uniform(low=0.0, high=2*math.sqrt(average/K), size=(num,K))

print('Running current symnmf_c.symnmf with default signature...')
H_out_curr = curr.symnmf(H0.tolist(), W_curr, K, num)
for row in H_out_curr:
    print(','.join('%.4f' % x for x in row))

if prev:
    print('\nRunning previous symnmfmodule.symnmf with same inputs...')
    W_prev = prev.norm(elements.tolist(), num, d)
    H_out_prev = prev.symnmf(H0.tolist(), W_prev, K, num)
    for row in H_out_prev:
        print(','.join('%.4f' % x for x in row))
else:
    print('Previous module not available')
