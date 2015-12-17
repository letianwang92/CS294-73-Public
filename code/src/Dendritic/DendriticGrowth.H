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
#include "WriteRectMDArray.H"
#include "VisitWriter.H"
#define P_NUM 10
#define DIM 2
using namespace std;

/// Class that specifies the right hand side of partial differential equations
class DendriticGrowth
{
public:
  DendriticGrowth();
  // RHS
  void operator()(DendriticShift& a_k, 
                     const Real& a_time, const Real& dt, 
                     Dendritic& a_state);
};
#endif
