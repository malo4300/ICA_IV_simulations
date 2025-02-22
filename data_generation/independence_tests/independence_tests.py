import numpy as np
from fcit import fcit
import pandas as pd
from tqdm import tqdm

def p_vals_conditional(data, signals, n, J):
    T = data[:,-2]
    controlls = data[:, 0:(J-3)]
    p_val = []
    for j in range(J):
        p_val.append(fcit.test(T.reshape(n,1), signals[:,j].reshape(n,1),controlls))
    return p_val

def p_vals_unconditional(data, signals, n, J):
    T = data[:,-2]
    p_val = []
    for j in range(J):
        p_val.append(fcit.test(T.reshape(n,1), signals[:,j].reshape(n,1), ))
    return p_val