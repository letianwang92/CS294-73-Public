#ifndef _BOX_H_
#define _BOX_H_
#include "Point.H"
class Box
{
public:
  inline Box();
  /// constructor for nontrivial Box. 
  inline Box(const Point& a_lowCorner,const Point& a_highCorner);
  /// copy constructor.
  inline Box(const Box& a_Box);
  ///  Computes the intersection of the Box with rightBox.
  inline Box operator*(const Box& a_rightBox) const;
  /// Computes shifted Box by a_offset in direction a_direction.
  inline Box shift(int a_direction, int a_offset);
  /// Computes shifted Box by a_pt.
  inline Box shift(const Point& a_pt);
  /// Grow in all of the coordinate directions by a_numpoints.
  inline Box grow(int a_numpoints);
  /// Grow in each coordinate direction by an amount given by the component of a_pt.
  inline Box grow(const Point& a_pt);
  /// coarsen in all of the coordinate directions by a_numpoints.
  inline Box coarsen(int a_numpoints);
  /// coarsen in each coordinate direction by an amount given by the component of a_pt.
  inline Box coarsen(const Point& a_pt);
 /// refine in all of the coordinate directions by a_numpoints.
  inline Box refine(int a_numpoints);
  /// refine in each coordinate direction by an amount given by the component of a_pt.
  inline Box refine(const Point& a_pt);
  /// Access functions for low, high bounds of Box.
  inline const Point& getLowCorner() const {return m_lowCorner;};
  inline const Point& getHighCorner() const {return m_highCorner;}; 
  /// Computes a linear index corresponding to the input array. 
  inline int linearIndex(const Point& a_pt) const;
  /// Computes the size of the Box.
  inline const int& sizeOf() const{return m_size;};
  /// When are two boxes equal.
  inline bool operator==(const Box& a_rhsBox) const
  {return ((m_lowCorner == a_rhsBox.m_lowCorner) && (m_highCorner == a_rhsBox.m_highCorner));};
  /// Contains the Point a_pt.
  inline bool contains(const Point& a_pt) const;  
  /// Get linear index of a Point in a Box.
  inline int getIndex(const Point& a_pt) const
{
  int factor = 1;
  int linIndex = a_pt[0] - m_lowCorner[0];
  for (int i = 1;i < DIM;i++)
    { 
      factor = factor*(m_highCorner[i-1] - m_lowCorner[i-1]+1);
      linIndex = linIndex + (a_pt[i] - m_lowCorner[i])*factor;
    }
  return linIndex;
};
  /// Says that you have not yet iterated through all ofthe points of the Box.
  inline bool notDone(const Point& a_pt) const;
  /// iteration through the points in a Box. a_pt is incremented to the next point
  /// in the Box.
  inline void increment(Point& a_pt) const;
  /// Get Point corresponding to a linear index in [0, ... sizeOf()-1] inside the Box
  inline Point getPoint(int k) const;
private:  
  Point m_lowCorner,m_highCorner;  
  int m_size;
};
#include "BoxImplem.H"
#endif
