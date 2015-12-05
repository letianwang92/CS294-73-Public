#include Dendritic.H

void DendriticShift::init(const Dendritic& a_dendritic)
{
  assert(a_dendritic.m_isInitilized==1);
  m_box=a_dendritic.m_box;
  m_phiShift.define(m_box);
  m_uShift.define(m_box);
}

void DendriticShift::increment(double a_scale, const DendriticShift&  a_rhs);
{
  assert(m_box.sizeOf() ==a_rhs.m_box.sizeOf());
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt);m_box.increment(pt))
    {
     m_phiShift[pt]= m_phiShift[pt]+a_scale*a_rhs.m_phi[pt];
     m_uShift[pt]= m_uShift[pt]+a_scale*a_rhs.m_u[pt];
    }
}
void DendriticShift::operator*=(double a_scale)
{
  for(Point pt = m_box.getLowCorner(); m_box.notDone(pt);m_box.increment(pt))
    {
      m_phiShift[pt]=m_phiShift[pt]*a_scale;
      m_uShift[pt]=m_phiShift[pt]*a_scale;
    }
}

void DendriticShift::setToZero()
{

  for(Point pt = m_box.getLowCorner(); m_box.notDone(pt);m_box.increment(pt))
    {
      m_phiShift[pt]=0;
      m_uShift[pt]=0;
    }
 
}

// Implement the Dendritic Class
 Dendritic::Dendritic()
{
}

Dendritic::Dentritic(Box& a_box,RectMDArray<Real>& phi_int,RectMDArray<Real>& u_int)
{
  m_box=a_box;
  m_phi=phi_int;
  m_u=u_int;
  m_isInitialized=1;
}


void  Dendritic::increment(const DendriticShift& a_shift);
{
  assert(a_dendritic.m_isInitilized==1);
  assert(m_box.sizeOf() ==a_shift.m_box.sizeOf());
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt);m_box.increment(pt))
    {
     m_phi[pt]= m_phi[pt]+a_scale*a_shift.m_phi[pt];
     m_u[pt]= m_u[pt]+a_scale*a_shift.m_u[pt];
    }
}
