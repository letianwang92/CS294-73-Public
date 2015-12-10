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
using namespace std;

int main(int argc, char* argv[])
{
  /// Input from screen, TODO: add more inputs
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
  double h = 1./N;

  /// Initialize Dendritic class
  Dendritic d;
  d.m_h = h;
  int low[DIM] = {0,0};
  int high[DIM] = {static_cast<int>(N),static_cast<int>(N)};
  Box bx(low,high);
  d.m_box=bx;

  /// Initial conditions: Specify Phi and u fields
  if (test == 1)
  {

  }
  else if (test == 2) 
  {

  }
  else
  {

  }

  /// Time advancing
  RK4<Dendritic,DendriticGrowth,DendriticShift> integrator;
  MDWrite(&d.m_phi);
  for(int i=0; i<m; i++)
    {
      integrator.advance(time, dt, d);
      time = time + dt;
      cout << "time = " << time << "  dt " << dt << endl;
      MDWrite(&d.m_phi);
      if (time >= timeStop) 
          break;
    }
}
