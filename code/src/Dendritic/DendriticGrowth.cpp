#include "DendriticGrowth.H"
#include <math.h>

using namespace std;
Real dot(array<Real,DIM> a, array<Real,DIM> b)
{
  Real ret=0;
  for(int i=0;i<DIM;i++)
    {
      ret+=a[i]*b[i];
    }
  return ret;
}

DendriticGrowth::DendriticGrowth(){}

void DendriticGrowth::operator()(DendriticShift& a_k, 
                     const Real& a_time, const Real& dt, 
                     Dendritic& a_state)
{
  // Parameter initialization
  RectMDOperators Operator;
  Real pi = M_PI;
  Real theta;
  Real W;
  Real Wp;
  Real nu;
  Real h=a_state.m_h;
  Real D=a_state.m_D;
  Real tau=a_state.m_tau;
  Real beta=a_state.m_beta;
  Real eta=a_state.m_eta;
  Real um=a_state.m_um;
  Real W0=a_state.m_W0;
  Real mu=a_state.m_mu;
  Real a0=a_state.m_a0;
  Real theta0=a_state.m_theta0;
  Real lh = a_state.m_L;
  array<Real,DIM> g_phi;
  Box domain(a_state.m_box);
  RectMDArray <Real> a_phi(a_state.m_phi);
  RectMDArray <Real> a_u(a_state.m_u);
  RectMDArray <Real> W_square(domain);
  RectMDArray <Real> WX1(domain);
  RectMDArray <Real> WX2(domain);
  RectMDArray <Real> W2_LOfPhi(domain);
  RectMDArray <Real> LOfU(domain);
  RectMDArray <array <Real, DIM>> GOfPhi(domain);
  RectMDArray <Real> LOfPhi(domain);
  Box innerbox (domain.grow(-1));

  // Step 1 - Calculate different terms of RHS
  for (Point pt = innerbox.getLowCorner(); innerbox.notDone(pt);innerbox.increment(pt))
  {
    a_phi[pt] += a_k.m_phiShift[pt];
    a_u[pt] += a_k.m_uShift[pt];
    g_phi = Operator.getGradient(a_phi,pt,h);
    LOfU[pt]=Operator.getLaplacian(a_u,pt,h); 
    LOfPhi[pt]=Operator.getLaplacian(a_phi,pt,h); 

    int sign=0;
    if (g_phi[0]>0)
      sign=1;
    theta = atan(g_phi[1]/(1E-6+g_phi[0]))+pi*(1-sign);
    W=W0*(1+mu*cos(a0*(theta-theta0)));
    Wp=-W0*mu*a0*sin(a0*(theta-theta0));
    GOfPhi[pt]=g_phi;
    WX1[pt]=W*Wp*g_phi[1];
    WX2[pt]=W*Wp*g_phi[0];
    W_square[pt]=W*W;
    W2_LOfPhi[pt]=W_square[pt]*LOfPhi[pt];
  }

  // Step 2 - Time advancing
  for (Point pt = innerbox.getLowCorner(); innerbox.notDone(pt);innerbox.increment(pt))
  {
    array <Real,DIM> g_W2=Operator.getGradient(W_square,pt,h);
    nu=beta/pi*atan(eta*(um-a_u[pt]));
    a_k.m_phiShift[pt]=((a_phi[pt]*(1-a_phi[pt])*(a_phi[pt]-0.5+nu)
                       -Operator.getGradient(WX1,pt,h)[0]+Operator.getGradient(WX2,pt,h)[1]
		       +dot(GOfPhi[pt],g_W2)+W2_LOfPhi[pt]))*dt/tau;
    a_k.m_uShift[pt]=lh*a_k.m_phiShift[pt]+dt*D*LOfU[pt];
  }
}
