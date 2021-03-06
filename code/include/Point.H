#ifndef POINT_H
#define POINT_H //what is this for ?
#include <assert.h>
class Point
{
  public:
  /// Default Constructor
  inline Point();
  /// Constructor
  inline Point(const int a_tuple[DIM]);
  /// Copy constructor
  inline Point(const Point& a_pt);
  /// Adds two points
  inline Point operator+(const Point& a_rhsPoint) const;
  /// Subtracts two points
  inline Point operator-(const Point& a_rhsPoint) const;
  /// division by an integer - rounds down.
  inline Point operator/(int a_nref) const;
  /// Componentwise division by a Point - rounds down.
  inline Point operator/(const Point& a_pt) const;
  /// multiplication by an integer.
  inline Point operator*(int a_nref) const;
  /// Componentwise multiplication by a Point.
  inline Point operator*(const Point& a_pt) const;
  /// Update in place arithmetic operations.
  inline void operator*=(const Point& a_pt);
  inline void operator+=(const Point& a_pt);
  inline void operator-=(const Point& a_pt);
  inline void operator/=(const Point& a_pt);
   /// Update in place aritmetic operations. Adds integer to every component.
  inline void operator*=(int a_n);
  inline void operator+=(int a_n);
  inline void operator-=(int a_n);
  inline void operator/=(int a_n);
  /// are two points equal.
  inline bool operator==(const Point& a_pt) const;
  /// Get a component of a point.
  inline const int&  operator[](const int& a_index) const
{
  //assert(a_index < DIM);
  //assert(a_index >= 0);
  return m_tuple[a_index];
};
 private:
  int m_tuple[DIM];
};
inline Point getZeros();
inline Point getOnes();
inline Point getUnitv(int a_dir)
{
  //assert((a_dir < DIM) && (a_dir >= 0));
  int unitv[DIM];
  for (int k = 0 ; k < DIM; k++)
    {
      unitv[k] = 0;
    }
 unitv[a_dir] = 1;
 return Point(unitv);
};
#include "PointImplem.H"
#endif
