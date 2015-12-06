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
  
  DendriticGrowth test_dG;
  test_dG.m_h=h; // the unit cell distance
   test_dG.m_D=1;
   test_dG.m_tau=1;
   test_dG.m_beta=1;
   test_dG.m_eta=1;
   test_dG.m_um=1;
   test_dG.m_W0=1;
   test_dG.m_mu=1;
   test_dG.m_a0=1;
   test_dG.m_theta0=0;

   assert(test_dG.isDefined());
 
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
