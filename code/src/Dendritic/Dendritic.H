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
  Dendritic(Box& a_box, RectMDArray<Real>& phi_int, RectMDArray<Real>& u_int);
  bool isDefined();
  void increment(const DendriticShift& a_shift);
  bool m_isDefined;
  Box m_box;
  RectMDArray<Real> m_phi;
  RectMDArray<Real> m_u;

  // Constants
  Real m_h; // unit cell distance
  Real m_D; // thermal diffusion constant
  Real m_tau; // relaxiation time
  Real m_beta; // material parameter
  Real m_eta; // material parameter
  Real m_um; // melting temperature
  Real m_W0; // initial interfacial width
  Real m_mu; // modulation of interfacial width
  Real m_a0; // anisotropic mode number
  Real m_theta0; // orientation angle
};


#endif
