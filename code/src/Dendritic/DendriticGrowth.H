#ifndef _DendriticGrowth_H_
#define _DendriticGrowth_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "RectMDArray.H"
#include "RectMDOperators.H"
using namespace std;

class DendriticGrowth
{
public:
  DendriticGrowth();
  // Evolve one step forward
  void stepForward(const Real& a_time, const Real& dt);

private:
  // Calculate RHS
  void calcRHS();

  // Increase using RK4
  void RK4();

  RectMDArray<Real> m_phi;
  RectMDArray<Real> m_u;
  RectMDArray<Real> m_phiShift;
  RectMDArray<Real> m_uShift;
  RectMDArray<array<Real, 2> > m_RHS; //the first number in the array for phi; the second number for u
  
  //Constants
  bool m_isInitialized;
  Real m_D;
  Real m_tau;
  Real m_beta;
  Real m_eta;
  Real m_um;
  Real m_W0;
  Real m_mu;
  Real m_a0;
  Real m_theta0;
};
#endif
