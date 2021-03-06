#include "Dendritic.H"
#include <cassert>
#include <iostream>

void DendriticShift::init(const Dendritic& a_dendritic)
{
  assert(a_dendritic.m_isInitialized==1);
  m_box=a_dendritic.m_box;
  m_phiShift.define(m_box);
  m_uShift.define(m_box);
  setToZero();
}

void DendriticShift::increment(double a_scale, const DendriticShift&  a_rhs)
{
  assert(m_box.sizeOf() ==a_rhs.m_box.sizeOf());
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt);m_box.increment(pt))
  {
    m_phiShift[pt]= m_phiShift[pt]+a_scale*a_rhs.m_phiShift[pt];
    m_uShift[pt]= m_uShift[pt]+a_scale*a_rhs.m_uShift[pt];
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
  m_h = 0;
  m_D = 0;
  m_tau = 0;
  m_beta = 0;
  m_eta = 0;
  m_um = 0;
  m_W0 = 0;
  m_mu = 0;
  m_a0 = 0;
  m_theta0 = 0;
  m_L = 0;
}

Dendritic::Dendritic(Box& a_box, RectMDArray<Real>& phi_int, RectMDArray<Real>& u_int)
{
  m_box=a_box;
  m_phi=phi_int;
  m_u=u_int;
  m_isInitialized=1;
}


void  Dendritic::increment(const DendriticShift& a_shift)
{
  assert(m_isInitialized==1);
  assert(m_box.sizeOf() ==a_shift.m_box.sizeOf());
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt);m_box.increment(pt))
  {
    m_phi[pt]= m_phi[pt]+a_shift.m_phiShift[pt];
    m_u[pt]= m_u[pt]+a_shift.m_uShift[pt];
  }
}
