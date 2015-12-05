#include "RectMDOperators.H"
#include <assert.h>

Real RectMDOperators::getLaplacian(const RectMDArray<Real>& a_field, const Point& a_p, Real a_h)
{
  //Check if the point is in the interior of the RectMDArray. Could be removed for performance issue
  Box innerBox(a_field.getBox().grow(-1));
  assert(innerBox.contains(a_p));
  
  Real result = 0;
  Point unitv;
  for(int i=0; i < DIM; i++)
  {
    unitv = getUnitv(i);
    result += a_field[a_p+unitv] + a_field[a_p-unitv];
  }
  result -= 2*DIM*a_field[a_p];
  result /= a_h*a_h;
  
  return result;
};

array<Real,DIM> RectMDOperators::getGradient(const RectMDArray<Real>& a_field, const Point& a_p, Real a_h)
{
  //Check if the point is in the interior of the RectMDArray. Could be removed for performance issue
  Box innerBox(a_field.getBox().grow(-1));
  assert(innerBox.contains(a_p));
  
  array<Real,DIM> result;
  Point unitv;
  for(int i=0; i<DIM; i++)
  {
    unitv = getUnitv(i);
    result[i] = (a_field[a_p+unitv] - a_field[a_p-unitv]) / (2*a_h);
  }
  
  return result;
};

Real RectMDOperators::getDivergence(const RectMDArray<array<Real, DIM> >& a_vecField, const Point& a_p, Real a_h)
{
  //Check if the point is in the interior of the RectMDArray. Could be removed for performance issue
  Box innerBox(a_vecField.getBox().grow(-1));
  assert(innerBox.contains(a_p));
  
  Real result=0;
  Point unitv;
  for(int i=0; i < DIM; i++)
  {
    unitv = getUnitv(i);
    result += a_vecField[a_p+unitv][i] - a_vecField[a_p-unitv][i];
  }
  result /= 2*a_h;
  
  return result;
};

void RectMDOperators::getLaplacianField(const RectMDArray<Real>& a_field, RectMDArray<Real>& a_lapField, Real a_h)
{
  
  Box innerBox(a_field.getBox().grow(-1));
  assert(innerBox == a_lapField.getBox());
  
  array<Point, DIM> unitvs;
  for(int i=0; i<DIM; i++)
  {
    unitvs[i] = getUnitv(i);
  }
  
  for(Point pIter=innerBox.getLowCorner(); innerBox.notDone(pIter); innerBox.increment(pIter) )
  {
    a_lapField[pIter] = 0;
    for(int i=0; i < DIM; i++)
    {
      a_lapField[pIter] += a_field[pIter+unitvs[i]] + a_field[pIter-unitvs[i]];
    }
    a_lapField[pIter] -= 2*DIM*a_field[pIter];
    a_lapField[pIter] /= a_h*a_h;
  }
};

void RectMDOperators::getGradientField(const RectMDArray<Real>& a_field, RectMDArray<array<Real,DIM> >& a_gradField, Real a_h)
{
  Box innerBox(a_field.getBox().grow(-1));
  assert(innerBox == a_gradField.getBox());
  
  array<Point, DIM> unitvs;
  for(int i=0; i<DIM; i++)
  {
    unitvs[i] = getUnitv(i);
  }
  
  for(Point pIter=innerBox.getLowCorner(); innerBox.notDone(pIter); innerBox.increment(pIter) )
  {
    for(int i=0; i < DIM; i++)
    {
      a_gradField[pIter][i] = (a_field[pIter+unitvs[i]] - a_field[pIter-unitvs[i]]) / (2*a_h);
    }
  }
};

void RectMDOperators::getDivergenceField(const RectMDArray<array<Real, DIM> >& a_vecField, RectMDArray<Real>& a_divField, Real a_h)
{
  Box innerBox(a_vecField.getBox().grow(-1));
  assert(innerBox == a_divField.getBox());
  
  array<Point, DIM> unitvs;
  for(int i=0; i<DIM; i++)
  {
    unitvs[i] = getUnitv(i);
  }
  
  for(Point pIter=innerBox.getLowCorner(); innerBox.notDone(pIter); innerBox.increment(pIter) )
  {
    a_divField[pIter] = 0;
    for(int i=0; i < DIM; i++)
    {
      a_divField[pIter] += a_vecField[pIter+unitvs[i]][i] - a_vecField[pIter-unitvs[i]][i];
    }
    a_divField[pIter] /= 2*a_h;
  }
};

