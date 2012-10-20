from mayavi import mlab
from numpy import *
from scitools.std import *
import os, glob, time, argparse
from mayavi.api import OffScreenEngine

parser = argparse.ArgumentParser()
parser.add_argument("-rm",action="store_true") 
args = parser.parse_args()

infile = open('initial.txt','r')
varlist = [line for line in infile]
infile.close()
Nx=0;Ny=0;Lx=0;Ly=0
for i in varlist:
	if i.startswith('Nx'):
		Nx = int(i.split('=')[-1])
	elif i.startswith('Ny'):
		Ny = int(i.split('=')[-1])
	elif i.startswith('Lx'):
		Lx = int(i.split('=')[-1])
	elif i.startswith('Ly'):
		Ly = int(i.split('=')[-1])
	else:
		print "hummm huuuumm... something is not right.."
		
dx = Lx/float(Ny+1)
dy = Ly/float(Ny+1)
files = glob.glob('texttmp*.txt')
sorted_files = sorted(files)
length = len(files)
f = open('u0.txt','r')
u0 = [ map(float,line.split()) for line in f ]
u0 = array(u0)
f.close()
i=1
print "Making picture files!"
mlab.figure()
#mlab.view(50, -100)
mlab.options.offscreen = True
X,Y = meshgrid(linspace(0,Lx,Nx),linspace(0,Ly,Ny))
s = mlab.mesh(X, Y, u0, color=(0.0,0.75,1.0))
mlab.savefig("wtmp0000.png")
mlab.clf()
for matrix in sorted_files:
	f = open(matrix,'r');l=[]
	l = [ map(float,line.split()) for line in f ]
	l=array(l)
	f.close()
	mlab.options.offscreen = True
	#mlab.figure()
	#mlab.view(50, -100)
	mlab.mesh(X, Y, l, color=(0.0,0.75,1.0))
	mlab.savefig("wtmp%.4d.png" %i)
	mlab.clf()
	i+=1
	print"picture %d out of %d done" %(i, length)

print "Making movie!"
movie("wtmp*.png", encoder='convert', fps=25, output_file='movie.gif')
'''
if args.rm:
	time.sleep(1)
	os.remove('initial.txt');os.remove('u0.txt')
	for i in glob.glob("wtmp*.png"):
		os.remove(i)
	for i in xrange(length):
		os.remove(files[i])
'''
