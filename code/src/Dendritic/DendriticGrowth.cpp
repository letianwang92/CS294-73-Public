
#include "DendriticGrowth.H"
#include <math.h>

using namespace std;
Real dot(array<Real,DIM> a, array<Real,DIM> b)
{
  Real ret=0;
  for(int i=0;i<DIM;i++)
    {
      ret+=a[i]+b[i];
    }
  return ret;
}

// this is default constructor, however, the parameter will actually require main program to access this class to define.
DendriticGrowth::DendriticGrowth()
{
}



void DendriticGrowth::operator()(DendriticShift& a_k, 
                     const Real& a_time, const Real& dt, 
                     Dendritic& a_state)
{
  RectMDOperators Operator;
  Real pi = M_PI;
  Real m_theta;
  Real W;
  Real W_prime;
  Real LOfPhi;
  Real nu;
  array<Real,DIM> g_phi;
  // setup constant
  Real m_h=a_state.m_h;
  Real m_D=a_state.m_D;
  Real m_tau=a_state.m_tau;
  Real m_beta=a_state.m_beta;
  Real m_eta=a_state.m_eta;
  Real m_um=a_state.m_um;
  Real m_W0=a_state.m_W0;
  Real m_mu=a_state.m_mu;
  Real m_a0=a_state.m_a0;
  Real m_theta0=a_state.m_theta0;
  Real tau_rev=1/m_tau;// the unit cell distance
  // Setup calculation variables
  Box domain(a_state.m_box);
  RectMDArray <Real> a_phi(a_state.m_phi);
  RectMDArray <Real> a_u(a_state.m_u);
  RectMDArray <Real> W_square(domain);
  RectMDArray <Real> WX1(domain);
  RectMDArray <Real> WX2(domain);
  RectMDArray <Real> W2_LOfPhi(domain);
  RectMDArray <Real> LOfU(domain);
  RectMDArray <array <Real, DIM>> GOfPhi(domain);

    for (Point pt = domain.getLowCorner(); domain.notDone(pt);domain.increment(pt))
      {
	a_phi[pt]+=a_k.m_phiShift[pt];
	a_u[pt]+=a_k.m_uShift[pt];
	g_phi=Operator.getGradient(a_phi,pt,m_h);
	LOfU[pt]=Operator.getLaplacian(a_u,pt,m_h); //actually whole field
	LOfPhi=Operator.getLaplacian(a_phi,pt,m_h); 
	GOfPhi[pt]=g_phi;
	int sign=0;
	if (g_phi[0]>0)
	{ sign=1;}
	m_theta=atan(g_phi[1]/g_phi[0])+pi*(1-sign);
	// g_phi is designed for faster access data. If not faster, remove it.
	W=m_W0*(1+m_mu*cos(m_a0*(m_theta-m_theta0)));
	W_prime=-m_W0*m_mu*m_a0*sin(m_a0*(m_theta-m_theta0));
	WX1[pt]=W*W_prime*g_phi[1];
	WX2[pt]=W*W_prime*g_phi[0];
	W_square[pt]=W*W;
	W2_LOfPhi[pt]=W_square[pt]*LOfPhi;
      }
  // here we calculated the phi gradient twice. which is time consuming for the 
  for (Point pt = domain.getLowCorner(); domain.notDone(pt);domain.increment(pt))
    {
      array <Real,DIM> g_W2=Operator.getGradient(W_square,pt,m_h);
      nu=m_beta/pi*atan(m_eta*(m_um-a_u[pt]));
      a_k.m_phiShift[pt]=dt*(tau_rev*(a_phi[pt]*(1-a_phi[pt])*(a_phi[pt]-1/2+nu)
				      +Operator.getGradient(WX1,pt,m_h)[0]-Operator.getGradient(WX2,pt,m_h)[1]
				  +dot(GOfPhi[pt],g_W2)+W2_LOfPhi[pt]));
      a_k.m_uShift[pt]=0.5*a_k.m_phiShift[pt]+dt*m_D*Operator.getLaplacian(a_phi,pt,m_h);
    }
   
    
}
