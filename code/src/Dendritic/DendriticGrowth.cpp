
#include "DendriticGworth.H"
#include <math.h>
using namespace std;
Real dot(array<Real,2> a, array<Real,2> b)
{
  Real ret=0;
  for( i=0;i<2;1++)
    {
      ret+=a[i]+b[i];
    }
  return ret;
}

DentriticGrowth::DentriticGrowth()
{
}


void DendriticGrowth::()(DendriticShift& a_k, 
                     const Real& a_time, const Real& dt, 
                     Dendritic& a_state);
{
  Real pi = M_PI;
  Real tau_rev=1/m_tau;
  Box domain=a_state.m_phi.getBox();
  RectMDArray <Real> W_square(domain);
  RectMDArray <Real> WX1(domain);
  RectMDArray <Real> WX2(domain);
  RectMDArray <Real> W2_LOfPhi(domain);
  RectMDArray <Real> LOfU(domain);
  RectMDArray <Real> GOfPhi(domain);

    for (Point pt = domain.getLowCorner(); domain.notDone(pt);domain.increment(pt))
      {
	g_phi=getGradient(m_phi,pt,a_h);
	LOfU(pt)=getLaplacian(m_u,pt,a_h); //actually whole field
	LOfPhi=getLaplacian(m_phi,pt,a_h); 
	GOfPhi(pt)=g_phi;
	int sign=0;
	if g_phi_x>0
	{ sign=1;}
	m_theta=atan(g_phi[1]/g_phi[0])+pi*(1-sign);
	// g_phi is designed for faster access data. If not faster, remove it.
	W=m_W0*(1+m_mu*cos(m_a0*(m_theta-m_theta0)));
	W_prime=-m_W0*m_mu*m_a0*sin(m_a0*(m_theta-m_theta0));
	WX1(pt)=W*W_prime*g_phi[1];
	WX2(pt)=W*W_prime*g_phi[0];
	W_square(pt)=W*W;
	W2_LOfPhi(pt)=W_square(pt)*LOfPhi;
      }
  // here we calculated the phi gradient twice. which is time consuming for the 
  for (Point pt = domain.getLowCorner(); domain.notDone(pt);domain.increment(pt))
    {
      array g_W2=getGradient(W_square(pt));
      nu=m_beta/pi*atan(m_mu(m_um-m_u(pt)));
      m_pRHS=tau_rev*(m_phi(pt)*(1-m_phi(pt))*(m_phi(pt)-1/2+nu)
			    +getGradient(WX1(pt))[0]-getGradient(WX2(pt))[1]
			    +dot(GOfPhi(pt),g_W2)+W2_LOfPhi(pt));
      m_uRHS=0.5*m_pRHS+m_D*getLaplacian(m_phi,pt,a_h);
    }
   
    
}
  
RectMDArray<Real> a_Shift DendriticGrowth::RK4(RectMDArray<Real> a_RHS,Real& dt)
{
  k1
  a_Shift=dt/6*(k1+2*k2+3*k3+k4)
}

