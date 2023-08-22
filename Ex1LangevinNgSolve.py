from ngsolve import *
from netgen.geom2d import unit_square
import numpy as np
from netgen.geom2d import SplineGeometry
import mylibmaterials as libmat


pi = np.pi
mu0 = 4*pi*1E-7
sigma = 5E6

h = 1e-3 # Height
b = 10e-3 # Length

freq = 60.
nper = 1   # Number of Periods
ptspp = 200 # Number of Points per Period
dt = (1/freq)/ptspp
Hpeak = 600.
t=-dt
a= 10
Ms=1.6e6

geo = SplineGeometry()
pnts =[(0,0),(b,0),(b,h),(0,h)]

p1,p2,p3,p4 = [geo.AppendPoint(*pnt) for pnt in pnts]

curves = [[["line",p1,p2],"bottom"],
          [["line",p2,p3],"right"],
          [["line",p3,p4],"top"],
          [["line",p4,p1],"left"]]

[geo.Append(c,bc=bc) for c,bc in curves]

omega = Mesh(geo.GenerateMesh(maxh=0.05*1e-3))

# H1-conforming finite element space
omegaP1 = H1(omega, order=1, dirichlet="right|top|left")


# Define Grid-Functions

H = GridFunction(omegaP1)
npH = H.vec.FV().NumPy()
npH[:] = Hpeak

M = GridFunction(omegaP1)
npM = M.vec.FV().NumPy()

Hold = GridFunction(omegaP1)
npHold = Hold.vec.FV().NumPy()

dBdH = GridFunction(omegaP1)
npdBdH = dBdH.vec.FV().NumPy()

# define trial- and test-functions
u,v = omegaP1.TnT()


f = LinearForm(omegaP1)
f += (sigma/dt*dBdH*v*Hold)*dx

Mm = BilinearForm(omegaP1)
Mm += (Grad(u)*Grad(v)+sigma/dt*dBdH*u*v)*dx 


vtk = VTKOutput(ma=omega,
                coefs=[H],
                names = ["H-field"],
                filename="gitignore/res",
                subdivision=0)
vtk.Do()

for n in range(nper*ptspp):
	npHold = npH[:]
	t = t + dt
	npH[:]=0
	H.Set(Hpeak*np.cos(2*pi*freq*t), BND)
	libmat.lib.pyLangevin(npM,npdBdH,npHold,1,len(npH),Ms,a)
	for i in range(len(npdBdH)):
		npdBdH[i] = mu0*(1+npdBdH[i])
	Mm.Assemble()
	f.Assemble()
	r = f.vec.CreateVector()
	r.data = f.vec - Mm.mat * H.vec
	H.vec.data += Mm.mat.Inverse(omegaP1.FreeDofs()) * r
	if (n%20 == 0):
		vtk.Do()

print("Finished")

