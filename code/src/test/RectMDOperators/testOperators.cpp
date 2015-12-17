#include <iostream>
#include <cmath>
#include "RectMDOperators.H"
using namespace std;

int main()
{
  RectMDOperators opSet;
  
  int N = 100;
  cout << "Please input the number of grid points in each direction:";
  cin >> N;
  Real a_h = 1.0/N;
  Point lCorner(getZeros());
  Point hCorner(getOnes()*N);
  Box b(lCorner, hCorner);
  Box bInner(b.grow(-1));
//  cout << '[' << b.getLowCorner()[0] <<',' << b.getLowCorner()[1] << "],[" << b.getHighCorner()[0] <<',' << b.getHighCorner()[1] << ']' <<endl;
//  cout << '[' << bInner.getLowCorner()[0] <<',' << bInner.getLowCorner()[1] << "],[" << bInner.getHighCorner()[0] <<',' << bInner.getHighCorner()[1] << ']' <<endl;
  
  // Create RectMDArrays for test
  RectMDArray<Real> field(b);
  RectMDArray<array<Real,DIM> > vecField(b);
  for(Point pIter = b.getLowCorner(); b.notDone(pIter); b.increment(pIter))
  {
    field[pIter] = 1;
    for(int i=0; i<DIM; i++)
    {
      field[pIter] *= sin(a_h*pIter[i]);
    }
    //cout << field[pIter] <<", ";

    for(int i=0; i<DIM; i++)
    {
      vecField[pIter][i] = sin(a_h*pIter[i]);
    }
  }
  
  // Calculate field results
  RectMDArray<Real> rLap(bInner); //result Laplacian
  RectMDArray<array<Real,DIM> > rGrad(bInner); //result Gradient
  RectMDArray<Real> rDiv(bInner); //result of divergence
  opSet.getLaplacianField(field, rLap, a_h);
  opSet.getGradientField(field, rGrad, a_h);
  opSet.getDivergenceField(vecField, rDiv, a_h);
  
  // Calculate point results
  RectMDArray<Real> rLapP(bInner); //result Laplacian 
  RectMDArray<array<Real,DIM> > rGradP(bInner); //result Gradient
  RectMDArray<Real> rDivP(bInner); //result of divergence
  for(Point pIter = bInner.getLowCorner(); bInner.notDone(pIter); bInner.increment(pIter))
  {
    rLapP[pIter] = opSet.getLaplacian(field, pIter, a_h);
    rGradP[pIter] = opSet.getGradient(field, pIter, a_h);
    rDivP[pIter] = opSet.getDivergence(vecField, pIter, a_h);
  }
  
  // Create correct results
  RectMDArray<Real> cLap(bInner); //corret Laplacian
  RectMDArray<array<Real,DIM> > cGrad(bInner); //corret gradient
  RectMDArray<Real> cDiv(bInner); //correct divergence
  for(Point pIter = bInner.getLowCorner(); bInner.notDone(pIter); bInner.increment(pIter))
  {
    cLap[pIter] = 1;
    for(int i=0; i<DIM; i++)
    {
      cLap[pIter] *= sin(a_h*pIter[i]);
    }
    cLap[pIter] *= -1*DIM;
    //cout << cLap[pIter] <<", " <<rLap[pIter] <<endl;
    
    for(int i=0; i<DIM; i++)
    {
      cGrad[pIter][i] = 1;
      for(int j=0; j<DIM; j++)
      {
        if(j==i)
        {
          cGrad[pIter][i] *= cos(a_h*pIter[j]);
        }
        else
        {
          cGrad[pIter][i] *= sin(a_h*pIter[j]);
        }
      }
    }
    
    cDiv[pIter] = 0;
    for(int i=0; i<DIM; i++)
    {
      cDiv[pIter] += cos(a_h*pIter[i]);
    }
  }
  
  // Get the error
  Real maxLapE = 0;
  Real maxLapPE = 0;
  Real maxGradE = 0;
  Real maxGradPE = 0;
  Real maxDivE = 0;
  Real maxDivPE = 0;
  for(Point pIter = bInner.getLowCorner(); bInner.notDone(pIter); bInner.increment(pIter))
  {
    if(abs(cLap[pIter]-rLap[pIter])>maxLapE)
    {
      maxLapE = abs(cLap[pIter]-rLap[pIter]);
    }
    if(abs(cLap[pIter]-rLapP[pIter])>maxLapPE)
    {
      maxLapPE = abs(cLap[pIter]-rLapP[pIter]);
    }
    
    for(int i=0; i<DIM; i++)
    {
      if(abs(cGrad[pIter][i]-rGrad[pIter][i])>maxGradE)
      {
        maxGradE = abs(cGrad[pIter][i]-rGrad[pIter][i]);
      }
      if(abs(cGrad[pIter][i]-rGradP[pIter][i])>maxGradPE)
      {
        maxGradPE = abs(cGrad[pIter][i]-rGradP[pIter][i]);
      }
    }
    
    if(abs(cDiv[pIter]-rDiv[pIter])>maxDivE)
    {
      maxDivE = abs(cDiv[pIter]-rDiv[pIter]);
    }
    if(abs(cDiv[pIter]-rDivP[pIter])>maxDivPE)
    {
      maxDivPE = abs(cDiv[pIter]-rDivP[pIter]);
    }
  }
  cout << "The max error of Laplacian field calculation when N=" << N << " is: " << maxLapE <<endl;
  cout << "The max error of Laplacian point calculation when N=" << N << " is: " << maxLapPE <<endl;
  cout << "The max error of gradient field calculation when N=" << N << " is: " << maxGradE <<endl;
  cout << "The max error of gradient point calculation when N=" << N << " is: " << maxGradPE <<endl;
  cout << "The max error of divergence field calculation when N=" << N << " is: " << maxDivE <<endl;
  cout << "The max error of divergence point calculation when N=" << N << " is: " << maxDivPE <<endl;

  return 1;
}
