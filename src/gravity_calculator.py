import numpy as np
import time
import multiprocessing
import subprocess
import os

import sys
import glob

systime = time.time()

#print "**************************************************"
#print "Calculation of Gravity Potentials"
#print ''

tsplit = time.time()
filename_grav_input = "grav_input.dat"

n_elem = int(subprocess.check_output("sed -n '$=' " + filename_grav_input, shell=True))
print("n_elem:"+str(n_elem))


n_proc = int(multiprocessing.cpu_count())
print("n_proc:"+str(n_proc))

n_lines = int(n_elem/n_proc)

if (float(n_elem)/float(n_proc)>0):
    n_lines = n_lines + (n_elem%n_proc)

print("n_lines:"+str(n_lines))



for filename in glob.glob('gchunk*'):
    os.remove(filename)

if n_proc > 26:
    sufflen = 2
else:
    sufflen = 1

os.system('split -a ' + str(sufflen) + ' -l ' + str(n_lines) + ' ' + filename_grav_input + ' gchunk ')


chunks = glob.glob('gchunk*')
for filename in chunks:
    os.rename(filename, filename + '.dat')

filenames = glob.glob('gchunk*')
print(filenames)
tsplit = time.time()-tsplit


#forming a command for computing
command_for_shell = ''
for i in range(0,n_proc):
    command_for_shell = command_for_shell + './gravcalc_parallel ' + filenames[i] + ' | '

command_for_shell = command_for_shell[0:len(command_for_shell)-2] + '\n'
print(command_for_shell)

os.system(command_for_shell)


print('Summation of the output grids...')
tsum = time.time()
GEOID = np.loadtxt(chunks[0] + '_GEOID.dat')
FA    = np.loadtxt(chunks[0] + '_FA.dat')
BGA   = np.loadtxt(chunks[0] + '_BGA.dat')

Uxx = np.loadtxt(chunks[0] + '_Uxx.dat')
Uyy = np.loadtxt(chunks[0] + '_Uyy.dat')
Uzz = np.loadtxt(chunks[0] + '_Uzz.dat')
Uzx = np.loadtxt(chunks[0] + '_Uzx.dat')
Uzy = np.loadtxt(chunks[0] + '_Uzy.dat')
Uxy = np.loadtxt(chunks[0] + '_Uxy.dat')

GEOID_cr = np.loadtxt(chunks[0] + '_GEOID_cr.dat')
FA_cr    = np.loadtxt(chunks[0] + '_FA_cr.dat')
BGA_cr   = np.loadtxt(chunks[0] + '_BGA_cr.dat')

Uxx_cr = np.loadtxt(chunks[0] + '_Uxx_cr.dat')
Uyy_cr = np.loadtxt(chunks[0] + '_Uyy_cr.dat')
Uzz_cr = np.loadtxt(chunks[0] + '_Uzz_cr.dat')
Uzx_cr = np.loadtxt(chunks[0] + '_Uzx_cr.dat')
Uzy_cr = np.loadtxt(chunks[0] + '_Uzy_cr.dat')
Uxy_cr = np.loadtxt(chunks[0] + '_Uxy_cr.dat')

for i in range(1,len(chunks)):
    GEOID = GEOID + np.loadtxt(chunks[i] + '_GEOID.dat')
    FA    = FA    + np.loadtxt(chunks[i] + '_FA.dat')
    BGA   = BGA   + np.loadtxt(chunks[i] + '_BGA.dat')

    Uxx = Uxx + np.loadtxt(chunks[i] + '_Uxx.dat')
    Uyy = Uyy + np.loadtxt(chunks[i] + '_Uyy.dat')
    Uzz = Uzz + np.loadtxt(chunks[i] + '_Uzz.dat')
    Uzx = Uzx + np.loadtxt(chunks[i] + '_Uzx.dat')
    Uzy = Uzy + np.loadtxt(chunks[i] + '_Uzy.dat')
    Uxy = Uxy + np.loadtxt(chunks[i] + '_Uxy.dat')

    GEOID_cr = np.loadtxt(chunks[0] + '_GEOID_cr.dat')
    FA_cr    = np.loadtxt(chunks[0] + '_FA_cr.dat')
    BGA_cr   = np.loadtxt(chunks[0] + '_BGA_cr.dat')

    Uxx_cr = Uxx_cr + np.loadtxt(chunks[i] + '_Uxx_cr.dat')
    Uyy_cr = Uyy_cr + np.loadtxt(chunks[i] + '_Uyy_cr.dat')
    Uzz_cr = Uzz_cr + np.loadtxt(chunks[i] + '_Uzz_cr.dat')
    Uzx_cr = Uzx_cr + np.loadtxt(chunks[i] + '_Uzx_cr.dat')
    Uzy_cr = Uzy_cr + np.loadtxt(chunks[i] + '_Uzy_cr.dat')
    Uxy_cr = Uxy_cr + np.loadtxt(chunks[i] + '_Uxy_cr.dat')

np.savetxt('grav_input_GEOID.dat', GEOID)
np.savetxt('grav_input_FA.dat', FA)
np.savetxt('grav_input_BGA.dat', BGA)

np.savetxt('grav_input_Uxx.dat', Uxx)
np.savetxt('grav_input_Uyy.dat', Uyy)
np.savetxt('grav_input_Uzz.dat', Uzz)
np.savetxt('grav_input_Uzx.dat', Uzx)
np.savetxt('grav_input_Uzy.dat', Uzy)
np.savetxt('grav_input_Uxy.dat', Uxy)

np.savetxt('grav_input_GEOID_cr.dat', GEOID_cr)
np.savetxt('grav_input_FA_cr.dat', FA_cr)
np.savetxt('grav_input_BGA_cr.dat', BGA_cr)

np.savetxt('grav_input_Uxx_cr.dat', Uxx_cr)
np.savetxt('grav_input_Uyy_cr.dat', Uyy_cr)
np.savetxt('grav_input_Uzz_cr.dat', Uzz_cr)
np.savetxt('grav_input_Uzx_cr.dat', Uzx_cr)
np.savetxt('grav_input_Uzy_cr.dat', Uzy_cr)
np.savetxt('grav_input_Uxy_cr.dat', Uxy_cr)

for filename in glob.glob('gchunk*'):
    os.remove(filename)


print("Time spent of operations with files: {:.2f} sec".format(((time.time()-tsum)+tsplit)/1.0))

print("Total execution time: {:.1f} min".format((time.time()-systime)/60.0))
