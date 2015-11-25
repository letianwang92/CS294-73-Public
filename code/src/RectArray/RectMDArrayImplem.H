#include "CH_Timer.H"
template <class T> RectMDArray<T>::RectMDArray()
{
  m_isDefined = false;
};
template <class T> void RectMDArray<T>::define(const Box& a_box)
{
  CH_TIME("define");
  if (!m_isDefined)
    {
      m_box=a_box;
      m_data = new T[a_box.sizeOf()];
      m_isDefined=true;
    }
  else
    {
      m_box = a_box;
      delete [] m_data;
      m_data = new T[a_box.sizeOf()];
      m_isDefined=true;
    }
};
template <class T> RectMDArray<T>::RectMDArray(const Box& a_box)
{
  CH_TIME("construct");
  m_box=a_box;
  m_data = new T[a_box.sizeOf()];
  m_isDefined=true;
};
template <class T> RectMDArray<T>::~RectMDArray()
{
  if (m_isDefined)
    {
      delete [] m_data;
    }
};
template <class T> RectMDArray<T>::RectMDArray(const RectMDArray<T>& a_srcArray)
{
  if (a_srcArray.m_isDefined)
    {
      m_box = a_srcArray.m_box;
      m_isDefined = true;
      m_data = new T[m_box.sizeOf()];
      for (int k = 0;k < m_box.sizeOf() ;k++)
        {
          m_data[k] = a_srcArray.m_data[k];
        }
    }
}; 
template <class T> RectMDArray<T>& RectMDArray<T>::operator=(const RectMDArray<T>& a_srcArray)
{
  if (m_isDefined)
    {
      m_isDefined = false;
      delete [] m_data;
      m_box = Box();
    }
  if (a_srcArray.m_isDefined)
    {
      m_box = a_srcArray.m_box;
      m_isDefined = true;
      m_data = new T[m_box.sizeOf()];
      for (int k = 0;k < m_box.sizeOf() ;k++)
        {
          m_data[k] = a_srcArray.m_data[k];
        }
    }
  return *this;
};       
template <class T> void RectMDArray<T>::setVal(const T& a_val)
{
  CH_TIME("setval");
  if (m_isDefined)
    {
      for (int k = 0; k < m_box.sizeOf();k++)
        {
          m_data[k] = a_val;
        }
    }
};
template <class T> T& RectMDArray<T>::operator[](const Point& a_iv) const
{
  // assert(m_isDefined);
  // CH_TIME("index");
  int k = m_box.getIndex(a_iv);
  return m_data[k];
};
template <class T> T* RectMDArray<T>::getPointer() const
{
  assert(m_isDefined);
  CH_TIME("index");
  return m_data;
};
