import numpy as np
import IMASW

filename = 'data.txt'

n = 4
vs = [50,100, 150, 250]
x_initial=  vs + [0.35]*n + [2000]*n + [1,5,10]  
lb = [40]*n  + [0.01]*n  + [1200]*n + [0.5]*(n-1)  
ub = [500]*n + [0.45]*n + [8000]*n + [30]*(n-1)   

test = IMASW.masw(filename, np.array(x_initial), np.array(lb), np.array(ub))
test.inverse()
