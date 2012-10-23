import scitools.std as sci
import glob,os,sys, argparse,random
from numpy import *
from mayavi import mlab

V_f_glob = lambda x,y: x*0
q_f_glob = lambda x,y,c=1.0: ones(shape(x))*c
f_f_glob = lambda x,y,t: x*0

def beach(x,y):
    q = x*(Lx - x)*y*(Ly - y)
    q = 0.5*x**2 + 0.5*y
    q += (q.min() + 0.1) 
    q /= (q.max() + 0.1)
    s = mlab.mesh(x, y, q)
    mlab.savefig("beach.png")
    return q

def box(x,y):
    """Returns the initial condition in the shape of a square plug with amplitude 2 """
    I = ones(shape(x))*0.5
    I_plot = I.copy()
    a = shape(x)
    sigma = 0.4
    #g=9.81
    for i in range(len(x[0])):
        for j in range(len(x[0])):
            I[j,i]-=0.3 if (abs(x[i,j]-Lx/2.0)<sigma) else 0
            I_plot[j,i] +=0.3 if (abs(x[i,j]-Lx/2.0)<sigma) else 0
    s = mlab.mesh(x, y, I_plot)
    mlab.savefig("box.png")
    return I#*g

def gauss_hill(x,y):
    sigma_x = 0.35
    sigma_y = 0.35
    x_0 = Lx/4.0
    y_0 = Ly/4.0
    q_plot = exp(-(((x-x_0)/(2*sigma_x))**2+((y-y_0)/(2*sigma_y))**2))*0.5
    q = (1-exp(-(((x-x_0)/(2*sigma_x))**2+((y-y_0)/(2*sigma_y))**2)))*0.5
    s = mlab.mesh(x, y,q_plot )
    mlab.savefig("mountain.png")
    return q
    

"""
def plog (x,y):
	""Returns the initial condition in the shape of a square plug with amplitude 2 ""
	I = zeros(shape(x))
	a = shape(x)
	sigma = 0.2
	#print x
	for i in range(len(x[0])):
	    for j in range(len(x[0])):
	    	I[j,i]=2 if (abs(x[i,j]-Lx/2.0)<sigma) else 0
	return I
"""
	
def const(x,y):
	I = zeros(shape(x))
	for i in range(len(x[0])):
	    for j in range(len(x[0])):
	    	I[i,j]=1.3 
	return I
	
def zero(x,y):
	I = zeros(shape(x))
	for i in range(len(x[0])):
	    for j in range(len(x[0])):
	    	I[i,j]=0
	return I

class wave_2d:

    def __init__(self,Lx, Ly, T, Nx, Ny, dt, b = 0,I_f=None, V_f=V_f_glob, q_f=q_f_glob, f_f=f_f_glob,exact=None, gauss=False, standing=False,plug=False,constant =False, rain=False):
        """
        I_f is the initial state-function, x,y
        V_f is the initial velocity-function, x,y
        q_f is the velocity of the wave at position x,y
        f_f is the source function , x,y,t
        gauss = True => I_f is set to a gauss function
        standing = True => I_f, V_f, q_f, f_f is set to fit a standing wave solution.
        """
        self.exact = exact
        self.standing = standing
        self.plug = plug

        if I_f==None and gauss==False and standing==False and plug==False and constant == False and rain==False:
            print "You failed to give wave properties"
            sys.exit(1)
        elif(gauss):
            self.I_f = self.gauss_f
            self.f_f = f_f
            self.V_f = V_f
            self.q_f = q_f
        elif(standing):
            self.I_f = self.standing_exact
            self.f_f = self.standing_source
            self.V_f = self.standing_V
            self.q_f = q_f
            self.exact = self.standing_exact
            self.w = 3.777775
            #self.w = 4.51179
            #self.mx = 5.1
            #self.my = 5.1
            self.mx = 5.03777
            self.my = 5.03777
            self.standing = standing
        elif(plug):
            self.I_f = self.plog
            #self.dx = self.dt
            #self.dy = self.dt
            self.b = 0
       	elif(constant):
            self.I_f = const
            self.V_f = zero
        elif(rain):
            self.I_f = self.gauss_f
            #self.I_f = lambda x,y: zeros(shape(x))
            self.f_f = self.rain
        else:
            self.I_f = I_f

		
        self.Lx = Lx
        self.Ly = Ly
        self.Nx = Nx
        self.Ny = Ny
        self.T = T
        self.dt = dt
        self.b = b
        #self.f_f = f_f
        #print "w = ", w, " mx = my = ", mx 
        self.x = linspace(0,Lx,Nx)
        self.y = linspace(0,Ly,Ny)
        self.dx = self.x[1]-self.x[0]
        self.dy = self.y[1]-self.y[0]
        self.Nt = int(T/dt)
        #self.T = Nt*dt
        self.t = linspace(0,self.T,self.Nt+1)
        #print (self.t[1]-self.t[0])
        #print self.t
        self.X1,self.Y1 = meshgrid(linspace(0,Lx,Nx),linspace(0,Ly,Ny))

        self.X,self.Y = meshgrid(linspace(0,Lx,Nx+2),linspace(0,Ly,Ny+2))
        X = self.X
        Y = self.Y
        self.V = V_f(X,Y)
        self.I = self.I_f(X,Y)
        print shape(X), shape(self.I)
        #print self.I
        self.q = q_f(X,Y)
        self.f = f_f(X,Y,0)
        self.C_y = dt/self.dy
        self.C_x = dt/self.dx
        print self.C_x, " C_x", self.dt, " dt", self.dx, " dx"
        
      

        
    def plog (self,x,y):
	"""Returns the initial condition in the shape of a square plug with amplitude 2 """
	I = zeros(shape(x))
	a = shape(x)
	sigma = 0.2
	
	for i in range(len(x[0])):
	    for j in range(len(x[0])):
	    	I[j,i]=2 if (abs(x[i,j]-Lx/2.0)<sigma) else 0
	return I

    def gauss_f(self,x,y):
    	"""Returns a gaussian shaped initial condition """
        sigma_x = 0.15
        sigma_y = 0.15
        x_0 = self.Lx/2.0
        y_0 = self.Ly/2.0
        return exp(-(((x-x_0)/(2*sigma_x))**2+((y-y_0)/(2*sigma_y))**2))

    def gauss_rain(self,x,y,sigma_x,sigma_y,x_0,y_0):
        """Returns a gaussian shaped initial condition """
        #sigma_x = 0.15
        #sigma_y = 0.15
        #x_0 = self.Lx/2.0
        #y_0 = self.Ly/2.0
        return exp(-(((x-x_0)/(2*sigma_x))**2+((y-y_0)/(2*sigma_y))**2))

    def standing_source(self,x,y,t):
    	"""Returns the source term for a standing wave solution of the equation. """
        w = self.w
        mx = self.mx
        my = self.my
        b = self.b
        Lx = self.Lx
        Ly = self.Ly
        q = self.q
        return (w*b*tan(w*t) - w**2 + q*pi**2*((mx/Ly)**2 + (my/Ly)**2))*self.standing_exact(x,y,t)
    
    def standing_V(self,x,y):
    	"""Returns the velocity distribution for the standing wave solution. """
        return -self.b*self.standing_exact(x,y,0)
    
    def standing_exact(self,x,y,t=0):
    	"""Returns the exact solution for the standing wave case at time level t."""
        w = self.w
        mx = self.mx
        my = self.my
        return exp(-self.b*t)*cos(w*t)*cos(mx*x*pi/self.Lx)*cos(my*y*pi/self.Ly)


    def make_exact(self):
        if(self.exact!=None):
            sol = self.exact(self.X,self.Y,0)
            savetxt('u0.txt',sol[1:-1,1:-1])
            for k in xrange(self.Nt):
                if (k+1)%3 == 0:
                    sol = self.exact(self.X, self.Y, self.t[k])
                    savetxt('texttmp%.4d.txt'%k,sol[1:-1,1:-1])
                
            
    def plot_exact(self):
        u_exact = self.standing_exact(self.X,self.Y)
        mlab.mesh(self.X1,self.Y1,u_exact[1:-1,1:-1])
        mlab.savefig("u_exact.png")
        savetxt("u_exact.txt",u_exact[1:-1,1:-1])
    
    def rain(self,x,y,t):
    	#Returns the source term which gives small raindroplets every 25th timestep """
        f = zeros((self.Nx+2,self.Ny+2))
        i = int(t/dt)
        if (i+1)%25 == 0:
            x0 = random.uniform(1,self.Lx-1)
            y0 = random.uniform(1,self.Ly-1)
            sigma_x = random.uniform(0.05,0.15)
            sigma_y = sigma_x
            f = self.gauss_rain(self.X, self.Y, x0,y0,sigma_x,sigma_y)
            """
            for o in xrange(1,self.Nx-1):
                for p in xrange(1,self.Ny-1):
                    f[o,p] = 0.1*self.gauss_rain(self.X[o,p],self.Y[o,p],x0,y0,sigma_x,sigma_y);
            """
        return f 


    def solve_num(self):
        Nx = self.Nx
        Ny = self.Ny
        Nt = self.Nt
        I = self.I
        q = self.q
        f = self.f
        dt = self.dt
        V = self.V
        b = self.b
        C_x = self.C_x
        C_y = self.C_y
        f_f = self.f_f
        X = self.X
        Y = self.Y
        t = self.t
        
        u = zeros((Nx+2,Ny+2), float)
        up = u.copy()
        upp = up.copy()
        if self.plug:
        	#dx = dt = 0.01
        	#q *=sqrt(2)
        	#C_y = 0; C_x = 1.0
        	b = 0
        #print shape(I)
        up[1:-1,1:-1] = I[1:-1,1:-1].copy()

        savetxt('u0.txt',up[1:-1,1:-1])
        
        """
        if self.standing:
            savetxt('u0_e.txt',up[1:-1,1:-1])
        """
           

        up[0,:] = up[2,:]
        up[:,0] = up[:,2]
        up[-1,:] = up[-3,:]
        up[:,-1] = up[:,-3]
        up[0][0] = up[2][2]
        up[0][-1] = up[2][-3]
        up[-1][-1] = up[-3][-3]
        up[-1][0] = up[-3][2]
        #making u^1
        if self.exact != None:
            if self.dt==1:
                errorfil = open('error.txt','w')
            errorfil = open('error.txt', 'a')

        for i in xrange(1,Nx+1):
            for j in xrange(1,Ny+1):
                x_para = ((q[i][j] + q[i+1][j])*(up[i+1][j] - up[i][j]) - (q[i-1][j] + q[i][j])*(up[i][j] - up[i-1][j]))
                y_para = ((q[i][j] + q[i][j+1])*(up[i][j+1] - up[i][j]) - (q[i][j-1] + q[i][j])*(up[i][j] - up[i][j-1]))
                rest = f[i][j] + 4*up[i][j] + 2*dt*V[i][j]*(b*dt-2)
                u[i][j] = 1.0/(4.0)*(C_x**2*x_para + C_y**2*y_para + rest)


        """


        x_para = (q[1:-1,1:-1] + q[2:,1:-1])*(up[2:,1:-1] - up[1:-1,1:-1]) - (q[2:,1:-1] + q[1:-1,1:-1])*(up[1:-1,1:-1] - up[:-2,1:-1])
        y_para = (q[1:-1,1:-1] + q[1:-1,2:])*(up[1:-1,2:] - up[1:-1,1:-1]) - (q[1:-1,:-2] + q[1:-1,1:-1])*(up[1:-1,1:-1] - up[1:-1,:-2])
        rest = 2*f[1:-1,1:-1] + 4*up[1:-1,1:-1] + 2*dt*V[1:-1,1:-1]*(b*dt-2)

        u[1:-1, 1:-1] = 1.0/4.0*(C_x**2*x_para +C_y**2*y_para + rest)
        """
     
        u[0,:] = u[2,:]
        u[:,0] = u[:,2]
        u[-1,:] = u[-3,:]
        u[:,-1] = u[:,-3]
        u[0][0] = u[2][2]
        u[0][-1] = u[2][-3]
        u[-1][-1] = u[-3][-3]
        u[-1][0] = u[-3][2]
        #making u^-1
        upp = 2*dt*V + u
        err = zeros(Nt)
        #vectorized:


        for k in xrange(Nt):
            
            x_para = (q[1:-1,1:-1] + q[2:,1:-1])*(up[2:,1:-1] - up[1:-1,1:-1]) - (q[:-2, 1:-1] + q[1:-1,1:-1])*(up[1:-1,1:-1] - up[:-2,1:-1])
            y_para = (q[1:-1,1:-1] + q[1:-1,2:])*(up[1:-1,2:] - up[1:-1,1:-1]) - (q[1:-1,:-2] + q[1:-1,1:-1])*(up[1:-1,1:-1] - up[1:-1,:-2])
            f = f_f(X,Y,t[k])
            rest = f[1:-1,1:-1] + 4*up[1:-1,1:-1] + upp[1:-1,1:-1]*(b*dt-2)
            
            u[1:-1,1:-1] = 1.0/(2+b*dt)*(C_x**2*x_para + C_y**2*y_para + rest)
            
            u[0,:] = u[2,:]
            u[:,0] = u[:,2]
            u[-1,:] = u[-3,:]
            u[:,-1] = u[:,-3]
            u[0][0] = u[2][2]
            u[0][-1] = u[2][-3]
            u[-1][-1] = u[-3][-3]
            u[-1][0] = u[-3][2]
            if self.exact != None:
                u_e = self.exact(X,Y,t[k])
                err[k] = sqrt(sum((u_e-u)**2)/(Nx*Ny))
            
            if (k+1)%3 == 0:    # and self.standing:
                #savetxt('u%.4d.txt'%k,u[1:-1,1:-1])
                #savetxt('u_e%.4d.txt'%k,u_e[1:-1,1:-1])
                #savetxt('texttmp%.4d.txt'%k, u_e[1:-1,1:-1]-u[1:-1,1:-1])
                #savetxt('texttmp%.4d.txt'%k,u_e[1:-1,1:-1])
                savetxt('texttmp%.4d.txt'%k,u[1:-1,1:-1])
            """
            elif (k+1)%3 ==0 and not self.standing:
                savetxt('texttmp%.4d.txt'%k,u[1:-1,1:-1])
            """
            
            upp = up.copy()
            up = u.copy()
        
        """
        a = open("max(err).txt", "a")
        a.write("%.10f\n" %max(err))
        a.close()
        """
        if self.exact != None:
            errorfil.write("error for dt=%4.2f: %6.4f \n"%(dt, max(err)))
            rel_error =  max(err)/dt**2
            errorfil.write("E/h**2: %6.4f \n"%rel_error )
            errorfil.close()
	

    

#--------Initialization---------------
parser = argparse.ArgumentParser()
parser.add_argument("-s",action="store_true", help="Scalar solver with damping and subsea geometry. If -s and -b Scalar solver without damping and subsea geometry") 
parser.add_argument("-standing",action="store_true", help="Choosing standing wave") 
parser.add_argument("-plug",action="store_true", help="Verify that the program can reproduce a square plug exactly in 1D")
parser.add_argument("-beach",action="store_true", help="Choose the sea bottom to look like a beach")
parser.add_argument("-box",action="store_true", help="Choose the sea bottom to look like a box")
parser.add_argument("-mountain",action="store_true", help="Choose the sea bottom to look like a mountain")
parser.add_argument("-gauss",action="store_true", help="Choosing gauss initiat condition")
parser.add_argument("-rain",action="store_true", help="Choosing the rain simulator") 
parser.add_argument("-b", type = float, dest="b", help="Damping coeff")
parser.add_argument("-Lx", type = int, dest="Lx", help="Size of area in x direction")
parser.add_argument("-Ly", type = int, dest="Ly", help="Size of area in y direction")
parser.add_argument("-T", type = int, dest="T", help="Number of timesteps")
parser.add_argument("-Nx", type = int, dest="Nx", help="Number of gridpoints in x direction")
parser.add_argument("-Ny", type = int, dest="Ny", help="Number of gridpoints i y direction")
parser.add_argument("-dt", type = float, dest="dt", help="timestep")
parser.add_argument("-movie", action="store_true",help="") 
parser.add_argument("-remove", action="store_true",help="remove the picture and text files produced from plotting")
parser.add_argument("-constant", action="store_true",help="verify that the program can produce a constant solution")
args = parser.parse_args()

b = args.b if args.b != None else 0.0
Lx = args.Lx if args.Lx != None else 5
Ly = args.Ly if args.Ly != None else Lx
T = args.T if args.T != None else 100
Nx = args.Nx if args.Nx != None else 25
Ny = args.Ny if args.Ny != None else Nx
#g = beach if args.beach==True elif box if args.box==True else q_f_glob
if args.beach:
    g = beach
elif args.box:
    g = box
elif args.mountain:
    g = gauss_hill




dx = Lx/float(Nx+1); dy = Ly/float(Ny+1);
buf=0.8
if args.dt <=0 or args.dt==None:
    dt = dx/sqrt(2)
elif args.dt != None:
    dx = args.dt*sqrt(2)
    dy = args.dt*sqrt(2)
    Nx = int(Lx/dx) + 1
    Ny = int(Ly/dy) + 1
    dt = args.dt
    c = 1. #sqrt(2)
if args.plug:
    """
    if args.dt!=None:
        dx = args.dt
        dy = args.dt
        Nx = int(Lx/dx) +1
        Ny = int(Ly/dy) +1
        dt = args.dt
    else:
    """
    dt = dx
    dy = dx
print "dt = ", dt, "dx = ", dx, "C_x = ", dt/dx
    
#print "Nx = ", Nx, "Ny = ", Ny, "dx = ",dx, "dy = ",dy, "dt = ",dt, "T = ",T
#print dt
#--------Write initial values to file--

init = ["Nx = "+str(Nx),"Lx = "+str(Lx),"Ny = "+str(Ny),"Ly = "+str(Lx)]
outfile = open('initial.txt','w')
for i in init:
    outfile.write(i);outfile.write(chr(10))
outfile.close()

#--------END Write initial values to file--

w = wave_2d(Ly,Lx,T,Nx,Ny,dt,b,standing=args.standing,gauss=args.gauss,plug=args.plug,constant = args.constant,q_f=g, rain=args.rain)


w.solve_num()
#w.make_exact()


if args.movie and not args.remove:
	#must install xvfb to run this
	cmd = 'xvfb-run --server-args="-screen 0 1024x768x24" python plot_wave.py'
	failure = os.system(cmd)
	if failure:
		print 'Execution of "%s failed\n' %cmd

if args.movie and args.remove:
	cmd = 'xvfb-run --server-args="-screen 0 1024x768x24" python plot_wave.py'
	failure = os.system(cmd)
	if failure:
		print 'Execution of "%s failed\n' %cmd
	text_files = glob.glob('texttmp*.txt')
	img_files = glob.glob("wtmp*.png")
	length_text_files = len(text_files)
	length_img_files = len(img_files)
	print "removing files"
	os.remove('initial.txt');os.remove('u0.txt')
	for i in xrange(length_img_files):
		os.remove(img_files[i])
	for i in xrange(length_text_files):
		os.remove(text_files[i])

if args.remove and not args.movie:
    text_files = glob.glob('texttmp*.txt')
    img_files = glob.glob("wtmp*.png")
    length_text_files = len(text_files)
    length_img_files = len(img_files)
    print "removing files"
    os.remove('initial.txt');os.remove('u0.txt')
    for i in xrange(length_img_files):
        os.remove(img_files[i])
    for i in xrange(length_text_files):
        os.remove(text_files[i])
