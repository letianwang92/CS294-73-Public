#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Real.H"
#include "RectMDArray.H"
#include "Dendritic.H"
#include "DendriticGrowth.H"
#include "RectMDOperators.H"
#include "WriteRectMDArray.H" 
#include "VisitWriter.H"
#include "RK4.H"
#include "CH_Timer.H"
#include "PowerItoI.H"
using namespace std;

int main(int argc, char* argv[])
{
  unsigned int M;
  unsigned int N;
  int test;
  double timeStop;
  cout << "input log_2(number of grid points)" << endl; 
  cin >> M;
  cout << "input test = 1, 2, 3" << endl;
  cin >> test;
  cout << "enter stopping time" << endl;
  cin >> timeStop;
  N = Power(2,M); // Total number of grid points
  Real h = 0.03; // grid spacing

  /// Initialize Dendritic class
  int low[DIM] = {0,0};
  Point lowCorner(low);
  int high[DIM] = {static_cast<int>(N),static_cast<int>(N)};
  double midpt[DIM] = {0.5*h*N,0.5*h*N};
  Box bx(low,high);
  RectMDArray<Real> phi_int(bx);
  RectMDArray<Real> u_int(bx);
  phi_int.setVal(0.);
  u_int.setVal(0.);
  // Set seed at the center of the domain, r = 0.1
  for (Point pt=lowCorner; bx.notDone(pt); bx.increment(pt))
  {
    double rd = sqrt((pt[0]*h-midpt[0])*(pt[0]*h-midpt[0])+(pt[1]*h-midpt[1])*(pt[1]*h-midpt[1]));
    if(rd < 0.1)
    {
      phi_int[pt] = 1.;
    }
  }
  Dendritic d(bx, phi_int, u_int);
  d.m_h = h;

  /// Initial conditions: Specify parameters.
  if (test == 1)
  {
    d.m_D = 1.; // thermal diffusion constant
    d.m_tau = 0.0003; // relaxiation time
    d.m_beta = 0.9; // material parameter
    d.m_eta = 10.0; // material parameter
    d.m_um = 1.0; // melting temperature
    d.m_W0 = 0.01; // initial interfacial width
    d.m_mu = 0.05; // modulation of interfacial width
    d.m_a0 = 6; // anisotropic mode number
    d.m_theta0 = 0; // orientation angle
    d.m_L = 2; // latent heat
  }
  else if (test == 2) 
  {
    d.m_D = 1.; // thermal diffusion constant
    d.m_tau = 0.0003; // relaxiation time
    d.m_beta = 0.9; // material parameter
    d.m_eta = 10.0; // material parameter
    d.m_um = 1.0; // melting temperature
    d.m_W0 = 0.01; // initial interfacial width
    d.m_mu = 0.05; // modulation of interfacial width
    d.m_a0 = 4; // anisotropic mode number
    d.m_theta0 = 0; // orientation angle
    d.m_L = 2; // latent heat
  }
  else
  {
    d.m_D = 1.; // thermal diffusion constant
    d.m_tau = 0.0003; // relaxiation time
    d.m_beta = 0.9; // material parameter
    d.m_eta = 10.0; // material parameter
    d.m_um = 1.0; // melting temperature
    d.m_W0 = 0.01; // initial interfacial width
    d.m_mu = 0.05; // modulation of interfacial width
    d.m_a0 = 6; // anisotropic mode number
    d.m_theta0 = 0; // orientation angle
    d.m_L = 4; // latent heat
  }

  /// Time advancing
  Real time = 0.;
  Real dt = 0.0001; // time step
  int m = 5000;
  RK4<Dendritic, DendriticGrowth, DendriticShift> integrator;
  MDWrite(d.m_phi);
  for(int i=0; i<m; i++)
  {
    integrator.advance(time, dt, d);
    time = time + dt;
    if ((i % 10) == 0)
    {
      MDWrite(d.m_phi);
      cout << "time = " << time << "  dt " << dt << endl;
    }
    if (time >= timeStop) 
      break;
  }
}
