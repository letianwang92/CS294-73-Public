#include <iostream>
#include <cmath>
#include "RK4.H"
#include "Dendritic.H"
#include "DendriticGrowth.H"
#include "RectMDOperators.H"
#include "PowerItoI.H"
using namespace std;

int main()
{

  unsigned int M;
  unsigned int N;
  //  cout << "input log_2(number of grid points)" << endl; 
  //cin >> M;
  M=4;
  Real timeStop=10;
  N = Power(2,M);
  double h = 1./N;
  int low[DIM] = {0,0};
  int high[DIM] = {static_cast<int>(N),static_cast<int>(N)};

  // Initialize the dendriticfield
  Box bx(low,high);
  RectMDArray<Real> u_int(bx);
  RectMDArray<Real> phi_int(bx);
  Dendritic test_d(bx,phi_int,u_int);

  // Initialize shift class
  DendriticShift kIn,kOut;
  kIn.init(test_d);
  kOut.init(test_d);
  kIn.setToZero();
  
   test_d.m_h=h; // the unit cell distance
   test_d.m_D=0.1;
   test_d.m_tau=0.1;
   test_d.m_beta=0.1;
   test_d.m_eta=0.1;
   test_d.m_um=0.1;
   test_d.m_W0=0.1;
   test_d.m_mu=0.1;
   test_d.m_a0=0.1;
   test_d.m_theta0=0.1;
   test_d.m_L=2;

   assert(test_d.m_isInitialized);
 
  Real time = 0.;
  Real dt = 140*.025/N;
  int m = 5000;
  RK4<Dendritic,DendriticGrowth,DendriticShift> integrator;
  for(int i=0; i<m; i++)
    {
      integrator.advance(time, dt, test_d);
      time = time + dt;
      cout << "time = " << time << "  dt " << dt << endl;

      if (time >= timeStop) 
        {
          break;
        }
    }
}
