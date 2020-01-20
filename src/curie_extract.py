import sys
print "**************************************************"
print "Calculation of Curie depth temperature map"

CURIETEMP = float(sys.argv[1])
print "Curie temperature is set to " + str(CURIETEMP) + " degrees"
print "This is the name of the script: ", sys.argv[0]
if(len(sys.argv)!=7):
    exit(-1)


import numpy as np
import os
import subprocess
from scipy.interpolate import interp1d


lon_min=-13
lon_max=4
lat_min=49
lat_max=60


lon_min=int(sys.argv[2])
lon_max=int(sys.argv[3])
lat_min=int(sys.argv[4])
lat_max=int(sys.argv[5])


proj=str(sys.argv[6])

print proj

#CURIETEMP = 585

#print "Curie temperature: " + str(CURIETEMP)

column_T = []
column_Z = []
curr_X = 0.0
curr_Y = 0.0

column_X = []
column_Y = []
column_curiedepth = []

X_max = 0
Y_max = 0
Z_max = 0
with open('OUT_xyz/Trhopvels.xyz') as f:
    for line in f:
        current = line.split()
        if Y_max < float(current[1]):
	    Y_max = float(current[1])
	if X_max < float(current[0]):
	    X_max = float(current[0])
	if Z_max < float(current[2]):
	    Z_max = float(current[2])

print "X_max: " + str(X_max) + " Y_max: " + str(Y_max) + " Z_max: " + str(Z_max)

firstline = 1
with open('OUT_xyz/Trhopvels.xyz') as f:

    for line in f:

        current = line.split()

        X = float(current[0])
        Y = Y_max-float(current[1])
        Z = float(current[2])
        T = float(current[3])
        rho = float(current[4])
        p = float(current[5])

        if Z == Z_max:

            if firstline == 1:
                firstline = 0
                continue
            else:
                f = interp1d(column_T, column_Z)
                column_X.append(curr_X)
                column_Y.append(curr_Y)
                column_curiedepth.append(f(CURIETEMP))

            column_T = []
            column_Z = []
            curr_X = X
            curr_Y = Y


            column_T.append(T)
            column_Z.append(Z)



        else:

            column_T.append(T)
            column_Z.append(Z)

    else:
        f = interp1d(column_T, column_Z)
        column_X.append(curr_X)
        column_Y.append(curr_Y)
        column_curiedepth.append(f(CURIETEMP))

        column_T = []
        column_Z = []
        curr_X = X
        curr_Y = Y


        column_T.append(T)
        column_Z.append(Z)

#column_X, column_Y, column_curiedepth = zip(*sorted(zip(column_X, column_Y, column_curiedepth), key=lambda x: -x[1]))


X_nodes = len(np.unique(column_X))
Y_nodes = len(np.unique(column_Y))

print "X_nodes: " + str(X_nodes) + " Y_nodes: " + str(Y_nodes)

h = open('OUT_xyz/curie.xyz', 'w')
for i in range(len(column_curiedepth)):
    h.write(str(column_X[i]/1000.0) + ' ' + str(column_Y[i]/1000.0) + ' ' + str(column_curiedepth[i]/1000.0) +'\n')

h.close()

os.system('rm OUT_xyz/curie_geo.xyz\n')
os.system('GMT surface OUT_xyz/curie.xyz -R0/' + str(X_max/1000.0) + '/0/' + str(Y_max/1000.0) + ' -I' + str(X_nodes)+ '+/' + str(Y_nodes)+'+ -GOUT_xyz/curie_cart.grd\n')
os.system('GMT grdproject OUT_xyz/curie_cart.grd -I -Rd' + str(lon_min) + '/' + str(lon_max) + '/' + str(lat_min) + '/' + str(lat_max) + ' -J' + proj + ' -Fk -GOUT_xyz/curie_geo.grd -V\n')
os.system('GMT grd2xyz OUT_xyz/curie_geo.grd > OUT_xyz/curie_geo.xyz\n')

import os.path

if os.path.isfile("OUT_xyz/curie_geo.xyz"):
    print "Curie depth map created"

print "**************************************************"
