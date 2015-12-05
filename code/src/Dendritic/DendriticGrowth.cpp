#include "DendriticGworth.H"
#include <math.h>
M_Pi
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


void DentriticGrowth::Initialize(Point lo, Point hi,RectMDArray<Real> phi_int,RectMDArray<Real> u_int)
{
 
  Point ptmax;
  Point lowCorner(lo);
  Point highCorner(hi);
  Real h = 1./(hi[0]);
  Box domainWithGhost(lo,hi);
  Box domain = domainWithGhost.grow(-1);
  Real* phiPtr = new double[domainWithGhost.sizeOf()];
  m_phi.setVal(0.);
  
  for (Point pt=lowCorner; domainWithGhost.notDone(pt); domainWithGhost.increment(pt))
    {
      m_phi(pt) = phi_int(pt);
      m_u(pt) = u_int(pt);    
    }
  m_isInitialized=1;
}

void DendriticGrowth:: stepForward(const Real& a_time, const Real& dt)
//looks like the a_time is not used again.
{
  m_Shift=RK4(m_RHS,dt);
  m
		
  

}

void DendriticGrowth::calcRHS()
{
  Real pi = M_PI;
  Real tau_rev=1/m_tau;
  RectMDArray <Real> W_square(domain);
  RectMDArray <Real> WX1(domain);
  RectMDArray <Real> WX2(domain);
  RectMDArray <Real> W2_LOfPhi(domain);
  RectMDArray <Real> LOfU(domain);
  RectMDArray <Real> GOfPhi(domain);

  LOfU(pt)=getLapa(pt)
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
      m_RHS(pt)[0]=tau_rev*(m_phi(pt)*(1-m_phi(pt))*(m_phi(pt)-1/2+nu)
			    +getGradient(WX1(pt))[0]-getGradient(WX2(pt))[1]
			    +dot(GOfPhi(pt),g_W2)+W2_LOfPhi(pt));
	m_RHS(pt)[1]=0.5*m_RHS(pt)[0]+m_D*getLaplacian(m_phi,pt,a_h);

	}
   
    
}
  
void DendriticGrowth::RK4(Real& dt)
