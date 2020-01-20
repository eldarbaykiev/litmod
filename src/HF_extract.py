import numpy as np
import os
import subprocess
from scipy.interpolate import interp1d

import sys
print "**************************************************"
print "Calculation of surface heat flow map"

lon_min=int(sys.argv[1])
lon_max=int(sys.argv[2])
lat_min=int(sys.argv[3])
lat_max=int(sys.argv[4])
proj=str(sys.argv[5])


X_max = 0
Y_max = 0
Z_max = 0

with open('OUT_xyz/HF.xyz') as f:
    for line in f:
        current = line.split()
        if Y_max < float(current[1]):
	    Y_max = float(current[1])
	if X_max < float(current[0]):
	    X_max = float(current[0])
	if Z_max < float(current[2]):
	    Z_max = float(current[2])

print "X_max: " + str(X_max) + " Y_max: " + str(Y_max) + " Z_max: " + str(Z_max)


column_X = []
column_Y = []
column_HF = []


firstline = 1
with open('OUT_xyz/HF.xyz') as f:
    for line in f:
        l = line.split(' ')
        l_r = list(filter(lambda x : x != '', l))
        #print l_r

        column_X.append(float(l_r[0]))
        column_Y.append(float(l_r[1]))
        column_HF.append(float(l_r[2][:-2]))



X_nodes = len(np.unique(column_X))
Y_nodes = len(np.unique(column_Y))

print "X_nodes: " + str(X_nodes) + " Y_nodes: " + str(Y_nodes)

h = open('OUT_xyz/HF_out.xyz', 'w')
for i in range(len(column_HF)):
    h.write(str(column_X[i]/1000.0) + ' ' + str(Y_max/1000.0-column_Y[i]/1000.0) + ' ' + str(column_HF[i]) +'\n')

h.close()


os.system('GMT makecpt -Cjet -T' + str(min(column_HF)) + '/' + str(max(column_HF)) + '/0.01 > OUT_xyz/HF.cpt\n')
os.system('GMT surface OUT_xyz/HF_out.xyz -R0/' + str(X_max/1000.0) + '/0/' + str(Y_max/1000.0) + ' -I' + str(X_nodes)+ '+/' + str(Y_nodes)+'+ -GOUT_xyz/HF_cart.grd\n')
os.system('GMT grdproject OUT_xyz/HF_cart.grd -I -Rd' + str(lon_min) + '/' + str(lon_max) + '/' + str(lat_min) + '/' + str(lat_max) + ' -J' + proj + ' -Fk -GOUT_xyz/HF_geo.grd -V\n')
#os.system(GMT + ' grdimage HF_geo.grd -Rd' + str(lon_min) + '/' + str(lon_max) + '/' + str(lat_min) + '/' + str(lat_max) + ' -J' + proj + ' -CHF.cpt -P -K > HF_geo.ps\n')
#os.system(GMT + ' pscoast -Dl -Rd' + str(lon_min) + '/' + str(lon_max) + '/' + str(lat_min) + '/' + str(lat_max) + ' -J' + proj + ' -Ba5g5f1/a5g5f1WESN -V -W1/0.5 -O -K >> HF_geo.ps\n')
#os.system(GMT + ' psscale -D5.3i/5i/5c/0.7c -CHF.cpt -I -B0.01:"HF":/:"w": -O >> HF_geo.ps\n')
#os.system(GMT + ' ps2raster HF_geo.ps\n')
os.system('GMT grd2xyz OUT_xyz/HF_geo.grd > OUT_xyz/HF_geo.xyz\n')

import os.path

if os.path.isfile("OUT_xyz/HF_geo.xyz"):
    print "Surface heatflow map created"
print "**************************************************"
