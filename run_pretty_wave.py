import os
from numpy import *
programName = 'pretty_wave_2d_plug_fredag.py'
#programName = 'copy_pretty.py'

dt = [1, 0.5, 0.2, 0.1, 0.02, 0.01]#, 0.005]
#w = linspace(4.50800, 4.51600, 20)
#mx = linspace(4.9,5.377777,10)

for i in range(len(dt)):
    cmd = 'python %s -standing -T 10 -Nx 100 -remove -dt %g ' % (programName, dt[i])
    #print cmd
    os.system(cmd)
