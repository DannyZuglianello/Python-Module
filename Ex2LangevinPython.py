import numpy as np
import matplotlib.pyplot as plt
import mylibmaterials as libmat

fig1, ax1 = plt.subplots()
ax1.set_title("Langevin")
ax1.set(xlabel="H [A/m]",ylabel="dBdH [H/m]")
fig2, ax2 = plt.subplots()
ax2.set_title("Langevin")
ax2.set(xlabel="H [A/m]",ylabel="B [T]")
# Data
mu0 = 4*np.pi*1e-7
Hpeak = 500
a=[1,10,100]
Ms=1.6e6
freq = 60
npts = 400
direction = np.array([1,1,1], dtype = np.float64)
direction = direction/np.linalg.norm(direction)
dim = len(direction)

eyes = np.array([libmat.eye(dim) for j in range(npts)], dtype = np.float64)

H =np.array([[[Hpeak*np.cos(2*np.pi*freq*t)*direction[i] for i in range(dim)] for t in np.linspace(0,1/freq,num=npts)] for j in range(len(a))], dtype = np.float64)
B = np.zeros((len(a),npts,dim))
M = np.zeros((len(a),npts,dim))
dMdH = np.zeros((len(a),npts,int(dim*(dim+1)/2)))
dBdH = np.zeros((len(a),npts,int(dim*(dim+1)/2)))

for i in range(len(a)):

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

ax1.legend()
ax2.legend()
fig1.savefig("./gitignore/dBdHH.png")
fig2.savefig("./gitignore/BH.png")
