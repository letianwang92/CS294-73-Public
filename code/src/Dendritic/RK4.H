#ifndef _RK4_H_
#define _RK4_H_

///  generic explicit RK4 algorithm
/**  class X is the type of your state data. class F provides an evaluation
     of the X' for a specified state of X and the time.  If your X' is 
     determined entirely from the current state values then F can ignore.
     F implements the function:
     void operator(X& k, double a_time, double a_dt, const X& a_state, const X& a_shift);

     F interprets the real state as a_state+a_shift

     Note that F includes the customary factor of dT in the RK4 algorithm. (i.e. this algorithm assumes k_i = F(...) NOT k_i = dt*F(...)
*/
template <class X, class F>
class RK4
{
public:
  void advance(double a_time, double a_dt, X& a_state);
protected:
  F m_f;
};

//==============================================



template <class X, class F>
void RK4<X, F>::advance(double a_time, double a_dt, X& a_state)
{
}
  
#endif 
  
