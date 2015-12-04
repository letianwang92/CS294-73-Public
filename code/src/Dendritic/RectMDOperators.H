#ifndef _MDOperator_H_
#define _MDOperator_H_
#include "RectMDArray.H"
using namespace std;

class RectMDOperators
{
/// This class is a set of operators that act on RectMDArray and a given point
public:
  // Laplacian operator
  Real getLaplacian(const RectMDArray<Real>& a_field, Point a_p, Real a_h) 
  // Gradient operator
  array<Real,2> getGradient(const RectMDArray<Real>& a_field, Point a_p, Real a_h) 
  // Divergence operator
  Real getDivergence(const RectMDArray<array<Real, 2> >& a_vecField, Point a_p, Real a_h) 

};
#endif
