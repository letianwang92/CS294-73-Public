#ifndef _DendriticGrowth_H_
#define _DendriticGrowth_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "RectMDArray.H"
#include "RectMDOperators.H"
#include "Dendritic.H"
using namespace std;



class DendriticGrowth
{
public:
  DendriticGrowth();
 
  // RHS
  
  void operator()(DendriticShift& a_k, 
                     const Real& a_time, const Real& dt, 
                     Dendritic& a_state);
  //Constants
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
