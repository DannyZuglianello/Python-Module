#ifndef PYMODULE_H
#define PYMODULE_H


#ifdef __cplusplus
   extern "C" {
#endif

void pyLangevin(double* M, double* dMdH, double* H, int dim, int nH ,double Ms, double a);
void pyInterpolation(double *Y,double *dYdX, double *X,double *Ydata,double *Xdata,int nXdata, int dimX,int nX);

#ifdef __cplusplus
   }
#endif

#endif



