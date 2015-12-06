#ifndef _Dendritic_H_
#define _Dendritic_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "RectMDArray.H"
#include "RectMDOperators.H"
#define P_NUM 10
#define DIM 2
using namespace std;

class Dendritic;
class DendriticShift
{
public:
  Box m_box;
  RectMDArray<Real> m_phiShift;
  RectMDArray<Real> m_uShift;
  void init(const Dendritic& a_dendritic);
  void increment(double a_scale, const DendriticShift& a_rhs );
  void operator*=(double a_scale);
  void setToZero();
};

class Dendritic
{
public:
  Dendritic();
  Dendritic(Box& a_box,RectMDArray<Real>& phi_int,RectMDArray<Real>& u_int);
  bool isDefined();
  void increment(const DendriticShift& a_shift);
  bool m_isDefined;
  Box m_box;
  RectMDArray<Real> m_phi;
  RectMDArray<Real> m_u;

    //Constant
  Real m_h; // the unit cell distance
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
