#ifndef _RECTMDARRAY_H_
#define _RECTMDARRAY_H_
#include "Box.H"
#include <assert.h>
template <class T> class RectMDArray 
{
public:
  /// Default constructor.
  RectMDArray();
  /// Constructs an RectMDArray over the Box a_box;
  RectMDArray(const Box& a_box);
  /// Defines a default-constructed RectMDArray.  
  void define(const Box& a_box);
/// Copy constructor.
  RectMDArray(const RectMDArray<T>& a_srcArray);
  /// Destructor.
  ~RectMDArray(); 
  /// Setting an RectMDArray to a constant value.
  void setVal(const T& a_val);
  /// Assignment operator.
  RectMDArray& operator=(const RectMDArray<T>& a_srcArray);
  /// Gets box over which array is defined.
  Box getBox() const {return m_box;};
  /// Indexing operator.
  inline T& operator[](const Point& a_iv) const;
  /// get raw pointer.
  inline T* getPointer() const;
  /// Index with integer.
  inline T& operator[](int a_index) const 
  { return m_data[a_index];}

  T& val(int a_index) 
  {return m_data[a_index];}
private:
  T* m_data;
  Box m_box;
  bool m_isDefined;
};
#include "RectMDArrayImplem.H"  
#endif
