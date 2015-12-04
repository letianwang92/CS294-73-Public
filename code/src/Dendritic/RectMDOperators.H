#ifndef _MDOperator_H_
#define _MDOperator_H_
#include "RectMDArray.H"
using namespace std;

class RectMDOperators
{
/// This class is a set of operators that act on RectMDArray and a given point
public:
  // Laplacian operator on a single point
  Real getLaplacian(const RectMDArray<Real>& a_field, const Point& a_p, Real a_h) 
  // Gradient operator on a single point
  array<Real,2> getGradient(const RectMDArray<Real>& a_field, const Point& a_p, Real a_h) 
  // Divergence operator on a single point
  Real getDivergence(const RectMDArray<array<Real, 2> >& a_vecField, const Point& a_p, Real a_h) 
  
  // Laplacian operator of the whole array
  void getLaplacianField(const RectMDArray<Real>& a_field, RectMDArray<Real>& a_lapField, Real a_h) 
  // Gradient operator of the whole array
  void getGradient(const RectMDArray<Real>& a_field, RectMDArray<array<Real,2> >& a_gradField, Real a_h) 
  // Divergence operator of the whole array
  void getDivergence(const RectMDArray<array<Real, 2> >& a_vecField, RectMDArray<Real>& a_divField, Real a_h) 

};
#endif
