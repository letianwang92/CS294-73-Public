#ifndef _MDOperator_H_
#define _MDOperator_H_
#include <array>
#include "Real.H"
#include "RectMDArray.H"
using namespace std;

/// This class is a set of operators that act on RectMDArray and a given point. All of the operators are dimension independent.
class RectMDOperators
{
public:
  // Laplacian operator on a single point
  Real getLaplacian(const RectMDArray<Real>& a_field, const Point& a_p, Real a_h);
  // Gradient operator on a single point
  array<Real,DIM> getGradient(const RectMDArray<Real>& a_field, const Point& a_p, Real a_h);
  // Divergence operator on a single point
  Real getDivergence(const RectMDArray<array<Real, DIM> >& a_vecField, const Point& a_p, Real a_h);
  
  // Laplacian operator of the whole array. Only values on interior grids will be calculated.
  void getLaplacianField(const RectMDArray<Real>& a_field, RectMDArray<Real>& a_lapField, Real a_h);
  // Gradient operator of the whole array. Only values on interior grids will be calculated.
  void getGradientField(const RectMDArray<Real>& a_field, RectMDArray<array<Real,DIM> >& a_gradField, Real a_h);
  // Divergence operator of the whole array. Only values on interior grids will be calculated.
  void getDivergenceField(const RectMDArray<array<Real, DIM> >& a_vecField, RectMDArray<Real>& a_divField, Real a_h);

};
#endif
