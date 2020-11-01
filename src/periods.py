
#LITMOD 4.0
#script to plot MINEOS grid output
#Eldar Baykiev, 2020

import numpy as np
import sys

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return a.flat[idx], idx


if len(sys.argv) != 2:
    print('ERROR')
    exit(-1)
    
PERIOD = float(sys.argv[1])

print(PERIOD)

try:
    grd_lonlat = np.loadtxt('RET.lonlat')
except:
    print('ERROR')
    exit(-1)
    
lon = grd_lonlat[:, 0]
lat = grd_lonlat[:, 1]
ray = lon * 0.0
love = lon * 0.0

for i in range(0, len(lon)):
    try:
        out_ray = np.loadtxt('%06d' % (i+1) + '_out_phase_ray.dat')
        out_love = np.loadtxt('%06d' % (i+1) + '_out_phase_love.dat')
    except:
        print('ERROR')
        exit(-1)
                
    dummy, ind_PERIOD = find_nearest(out_ray[:, 0], PERIOD)
            
    ray[i] = out_ray[ind_PERIOD, 1]
    love[i] = out_love[ind_PERIOD, 1]
            
ray_avg = np.average(ray)
love_avg = np.average(love)

ray_anom = ((ray - ray_avg) / ray) * 100.0
love_anom = ((love - love_avg) / love) * 100.0
        
with open('MINEOS_out_grid_ray_' + str(PERIOD) + '.xyz', 'w') as f_ray:
    with open('MINEOS_out_grid_love_' + str(PERIOD) + '.xyz', 'w') as f_love:
        for i in range(0, len(lon)):
            f_ray.write(str(lon[i]) + ' ' + str(lat[i]) + ' ' + str(ray[i]) + '\n')
            f_love.write(str(lon[i]) + ' ' + str(lat[i]) + ' ' + str(love[i]) + '\n')
            
with open('MINEOS_out_grid_ray_' + str(PERIOD) + '_anom.xyz', 'w') as f_ray:
    with open('MINEOS_out_grid_love_' + str(PERIOD) + '_anom.xyz', 'w') as f_love:
        for i in range(0, len(lon)):
            f_ray.write(str(lon[i]) + ' ' + str(lat[i]) + ' ' + str(ray_anom[i]) + '\n')
            f_love.write(str(lon[i]) + ' ' + str(lat[i]) + ' ' + str(love_anom[i]) + '\n')
