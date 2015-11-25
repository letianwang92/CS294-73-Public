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
  void operator()(RectMDArray& a_phi, RectMDArray& a_T,
                  const Real& a_time, const Real& dt);
};
#endif
