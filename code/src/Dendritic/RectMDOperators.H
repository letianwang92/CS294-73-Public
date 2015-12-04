#ifndef _MDOperator_H_
#define _MDOperator_H_
#include "RectMDArray.H"
using namespace std;

class RectMDOperators
{
/// This class is a set of operators that act on RectMDArray
public:
  // Laplacian operator
  void getLaplacian(const RectMDArray<Real>& a_field, RectMDArray<Real>& a_laplacian, Real a_h) 
  // Gradient operator
  void getGradient(const RectMDArray<Real>& a_field, RectMDArray<array<Real, 2> >& a_gradient, Real a_h) 
  // Divergence operator
  void getDivergence(const RectMDArray<array<Real, 2> >& a_vecField, RectMDArray<Real>& a_divergence, Real a_h) 

};
#endif
