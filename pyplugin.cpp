#include "pylibmaterials.h"
#include "../libmaterials/AuxiliaryMath.cpp"
#include "../libmaterials/ModelLangevin.cpp" 
#include "../libmaterials/TensorProduct.cpp" 
#include "../libmaterials/ModelInterpolation.cpp" 

void pyLangevin(double *M, double *dMdH, double *H,int dim,int nH, double Ms, double a)
{
   int* pdim = new int;
   *pdim = dim;
   int* pnH = new int;
   *pnH = nH;
   double* pMs = new double;
   *pMs = Ms;
   double* pa = new double;
   *pa = a;



   ModelLangevin model(pMs, pa);
   model.Model_Langevin(dMdH,M, H, pdim, pnH);


   delete pdim;
   delete pnH;
   delete pMs;
   delete pa;
}

void pyInterpolation(double *Y,double *dYdX, double *X, int dimX,int nX)
{
   int* pdimX = new int;
   *pdimX = dimX;
   int* pnX = new int;
   *pnX = nX;


   ModelInterpolation model(Y, X, pdimX);
   model.Model_Interpolation(dYdX, Y, X, pdimX, pnX);


   delete pdimX;
   delete pnX;

}
