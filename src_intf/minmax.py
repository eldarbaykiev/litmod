
#LITMOD 4.0
#gmt minmax substitute for interface
#Eldar Baykiev, 2020

import numpy as np
import sys

if len(sys.argv) != 2:
    print('ERROR')
    exit(-1)
    
filename = sys.argv[1]
import os
if not (os.path.isfile(filename)):
    print('ERROR')
    exit(-1)
    
try:
    grd = np.loadtxt(filename)
except:
    print('ERROR')
    exit(-1)
    
col_n = len(grd[0, :])
for i in range(len(grd[0, :])):
    print(str(np.min(grd[:, i])) + '\t' + str(np.max(grd[:, i])), end='\t')
    
print('\n', end='')



