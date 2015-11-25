#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Real.H"
#include "RectMDArray.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "CutoffKernel.H"
#include "WriteRectMDArray.H" 
#include "VisitWriter.H"
#include "RK4.H"
void outField(ParticleSet& p, int a_coarsenFactor)
{
  int coarsenFactor = a_coarsenFactor;
  Box bx = p.m_box.coarsen(coarsenFactor);
  double h = p.m_dx*coarsenFactor;
  RectMDArray<double> outVort(bx);
  int ipos[DIM];
  double xpos[DIM];
  double weight;
  Point e0 = getUnitv(0);
  Point e1 = getUnitv(1);
  
  outVort.setVal(0.);
  for (int k = 0; k < p.m_particles.size(); k++)
    {
      for (int l = 0; l < DIM; l++)
        {
          double newpos = p.m_particles[k].m_x[l]; 
          ipos[l] = newpos/h;
          xpos[l] = (newpos - ipos[l]*h)/h;
        }
      Point pt(ipos);
      assert(p.m_box.contains(pt));
      for (int l0=0; l0 < DIM;l0++)
        {
          for (int l1=0;l1 < DIM ; l1++)
            {
              outVort[pt+e0*l0 + e1*l1] += 
                (1.-xpos[0] + (2*xpos[0] - 1.)*l0)*
                (1.- xpos[1] + (2*xpos[1] - 1.)*l1)*p.m_particles[k].strength/coarsenFactor/coarsenFactor;
            }
        }
    }
  const char* foo = MDWrite(outVort);
};
int main(int argc, char* argv[])
{
  unsigned int M;
  unsigned int N;
  cout << "input log_2(number of grid points)" << endl; 
  cin >> M;
  cout << "input test = 1,2, other" << endl;
  int test;
  cin >> test;
  cout << "input particle refinement factor" << endl;
  unsigned int cfactor;
  cin >> cfactor;
  cout << "enter stopping time" << endl;
  double timeStop;
  cin >> timeStop;
  ParticleSet p;
  N = Power(2,M);
  double h = 1./N;
  Real hp = h/cfactor; //pow(h,4./3.);
  int Np = 1./hp;
  hp = 1./Np;
  double delta = h;
  int pcfactor = 4/cfactor;
  if (pcfactor < 1 ) pcfactor = 1;
  cout << "number of particles per cell = " << h*h/hp/hp << endl;
  shared_ptr<CutoffKernel> cutkptr = 
    shared_ptr<CutoffKernel>(new CutoffKernel(h,delta));
  shared_ptr<ConvKernel> convkerptr = dynamic_pointer_cast<ConvKernel>(cutkptr);
  p.m_hockney.define(convkerptr,h,M);
  p.m_dx = h;
  Real lowCorner[DIM];
  if (test == 1)
  {
    p.m_particles.resize(1);
      Particle& particle = p.m_particles[0];
      particle.m_x[0] = .5;
      particle.m_x[1] = .25;
      particle.strength = 1./h/h;
  }
  else if (test == 2) 
    {
      p.m_particles.resize(2);
      Particle& particle = p.m_particles[0];
      particle.m_x[0] = .5;
      particle.m_x[1] = .25;
      particle.strength = 1./h/h;
      Particle& particle2 = p.m_particles[1];
      particle2.m_x[0] = .5;
      particle2.m_x[1] = .75;
      particle2.strength = 1./h/h;
    }
  else
    {
      Real xp[DIM];
      for (int i = 0;i < Np;i++)
        {
          xp[0] = i*hp;
          for (int j = 0; j< Np; j++)
            {
              xp[1] = j*hp;
              Real dist1 = sqrt(pow(xp[0] - .375,2) + pow(xp[1] - .5,2));
              Real dist2 = sqrt(pow(xp[0] - .625,2) + pow(xp[1] - .5,2));
              if ((dist1 < .12 ) | (dist2 < .12))
            // if (dist1 < .1125 )
                {
                  Particle part;
                  part.m_x[0] = xp[0];
                  part.m_x[1] = xp[1];
                  part.strength = hp*hp/h/h;
                  p.m_particles.push_back(part);
                }
            }
        }
    }
  Real dx = 1./N;
  cout << "number of particles = " << p.m_particles.size() << endl;
  int low[DIM] = {0,0};
  int high[DIM] = {static_cast<int>(N),static_cast<int>(N)};
  Box bx(low,high);
  p.m_box=bx;
  lowCorner[0] = 0.;
  lowCorner[1] = 0.;
  ParticleShift kIn,kOut;
  kIn.init(p);
  kOut.init(p);
  kIn.setToZero();
  ParticleVelocities pv(); 
  Real time = 0.;
  Real dt = 140*.025/N;
  int m = 5000;
  
  RK4<ParticleSet,ParticleVelocities,ParticleShift> integrator;
  outField(p,pcfactor);
  PWrite(&p);
  for(int i=0; i<m; i++)
    {
      integrator.advance(time, dt, p);
      time = time + dt;
      cout << "time = " << time << "  dt " << dt << endl;
      outField(p,pcfactor);
      PWrite(&p);
      if (time >= timeStop) 
        {
          break;
        }
    }
}
