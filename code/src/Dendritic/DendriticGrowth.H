#ifndef _DendriticGrowth_H_
#define _DendriticGrowth_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "RectMDArray.H"
using namespace std;

class DendriticGrowth
{
public:
  DendriticGrowth();
  // Evolve one step forward with RK 4
  void stepForward(const Real& a_time, const Real& dt);

private:
  
  RectMDArray m_phi;
  RectMDArray m_u;
  RectMDArray m_phiShift;
  RectMDArray m_uShift;
  
  //Constants
  Real m_D;
  Real m_tau;
  Real m_beta;
  Real m_eta;
  Real m_um;
  Real W_0;
  Real m_mu;
  Real m_a0;
  Real m_theta0;
};
#endif
