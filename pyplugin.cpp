#include "pylibmaterials.h"
#include "../libmaterials/AuxiliaryMath.cpp"
#include "../libmaterials/ModelLangevin.cpp" 
#include "../libmaterials/TensorProduct.cpp" 

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

