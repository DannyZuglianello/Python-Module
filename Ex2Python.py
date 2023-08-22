import numpy as np
import matplotlib.pyplot as plt
import mylibmaterials as libmat

model_name = "Langevin" # Pick Langevin or Interpolation
dataset = 1 # Pick Interpolation Dataset 0 or 1

fig1, ax1 = plt.subplots()
ax1.set_title(model_name)
ax1.set(xlabel="H [A/m]",ylabel="dBdH [H/m]")
fig2, ax2 = plt.subplots()
ax2.set_title(model_name)
ax2.set(xlabel="H [A/m]",ylabel="B [T]")
# Data
mu0 = 4*np.pi*1e-7
Hpeak = 500

#Langevin Data
a=[1,10,100]
Ms=1.6e6
lena = len(a)
#Interpolation Data
if(dataset != 0):
	bhH=np.array([0.0, 8.06, 10.85, 13.12, 15.07, 16.79, 18.42, 19.93, 21.31, 22.65, 23.95, 25.24, 26.53, 28.25, 32.66, 43.00, 70.15, 182.84, 799.18, 842.93, 2503.30, 4861.50, 10077.00],dtype = np.float64)
	bhB=np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.8965, 1.9000, 1.9527, 1.9673, 1.9755],dtype = np.float64)
else:
	bhH=np.array([0,100,250,500,1000,3000],dtype = np.float64)
	bhB=np.array([0,0.7,1.2,1.45,1.55,1.7],dtype = np.float64)
ndata=len(bhH)

freq = 60
npts = 400
direction = np.array([1,1,1], dtype = np.float64)
direction = direction/np.linalg.norm(direction)
dim = len(direction)

eyes = np.array([libmat.eye(dim) for j in range(npts)], dtype = np.float64)
if(model_name.lower()=="langevin"):
	H =np.array([[[Hpeak*np.cos(2*np.pi*freq*t)*direction[i] for i in range(dim)] for t in np.linspace(0,1/freq,num=npts)] for j in range(lena)], dtype = np.float64)
	B = np.zeros((lena,npts,dim))
	M = np.zeros((lena,npts,dim))
	dMdH = np.zeros((lena,npts,int(dim*(dim+1)/2)))
	dBdH = np.zeros((lena,npts,int(dim*(dim+1)/2)))
else:
	H =np.array([[Hpeak*np.cos(2*np.pi*freq*t)*direction[i] for i in range(dim)] for t in np.linspace(0,1/freq,num=npts)], dtype = np.float64)
	B = np.zeros((npts,dim))
	M = np.zeros((npts,dim))
	dMdH = np.zeros((npts,int(dim*(dim+1)/2)))
	dBdH = np.zeros((npts,int(dim*(dim+1)/2)))
if(model_name.lower()=="langevin"):
    for i in range(lena):
        tH = H[i].reshape((dim*npts,))
        tM = M[i].reshape((dim*npts,))
        tdMdH = dMdH[i].reshape((int(dim*(dim+1)/2)*npts,))

        libmat.lib.pyLangevin(tM,tdMdH,tH,dim,npts,Ms,a[i])

        H[i] = tH.reshape((npts,dim))
        M[i] = tM.reshape((npts,dim))
        dMdH[i] = tdMdH.reshape((npts,int(dim*(dim+1)/2)))


        B[i] = mu0*(H[i]+M[i])
        dBdH[i] = mu0*(eyes+dMdH[i])

        Hh = [np.sum(H[i][j]*direction) for j in range(npts)]

        Bh = [np.sum(B[i][j]*direction) for j in range(npts)]

        dBdHh = [np.sum(libmat.prod(dBdH[i][j],direction)*direction) for j in range(npts)]

        ax1.plot(Hh,dBdHh,label="a = "+str(a[i]))
        ax2.plot(Hh,Bh,label="a = "+str(a[i]))
elif(model_name.lower()=="interpolation"):
    tH = H.reshape((dim*npts,))
    tB = B.reshape((dim*npts,))
    tdBdH = dBdH.reshape((int(dim*(dim+1)/2)*npts,))

    libmat.lib.pyInterpolation(tB,tdBdH,tH,bhB,bhH,ndata,dim,npts)

    H = tH.reshape((npts,dim))
    B = tB.reshape((npts,dim))
    dBdH = tdBdH.reshape((npts,int(dim*(dim+1)/2)))

    Hh = [np.sum(H[j]*direction) for j in range(npts)]

    Bh = [np.sum(B[j]*direction) for j in range(npts)]

    dBdHh = [np.sum(libmat.prod(dBdH[j],direction)*direction) for j in range(npts)]

    ax1.plot(Hh,dBdHh)
    ax2.plot(Hh,Bh)
else:
    print("Model " + model_name + " is not yet available. Pick between Langevin and Interpolation")

if(model_name.lower()=="langevin"):
	ax1.legend()
	ax2.legend()
if(model_name.lower()=="interpolation"):
	model_name = model_name + "data" + str(dataset)

fig1.savefig("./gitignore/dBdHH_"+model_name+".png")
fig2.savefig("./gitignore/BH_"+model_name+".png")
