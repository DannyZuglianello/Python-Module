import numpy as np
import ctypes
import pathlib

# Read Langevin Module

lib = ctypes.cdll.LoadLibrary('./gitignore/pylibmaterials.so')

lib.pyLangevin.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double]

lib.pyInterpolation = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_int]
    
# prod defines the matrix product of a symmetric matrix A of size 1x1, 2x2 or 3x3 in Voigt Form and a vector.

# Note dBdH = mu0*eye(n) + dMdH
def prod(A,b):
    temp = np.zeros(len(b), dtype = np.float64)
    for i in range(len(b)):
        temp[i] = temp[i] + A[i]*b[i]
        for j in range(len(b)):
            if i != j:
#                temp[j] = temp[j] + A[(int(len(b)*(len(b)+1)/2))-i-j]*b[i]
                temp[j] = temp[j] + A[len(A)-i-j]*b[i]
    return temp
    
def eye(n):
	temp = np.zeros(int((n**2+n)/2),dtype = np.float64)
	for i in range(n):
		temp[i] = 1.
	return temp
		
