// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.
//
// Imported from licorne

#include <dolfin/evolution/Time.h>

#include <dolfin/log/log.h>
#include <dolfin/math/basic.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Time::Time(real T0, real T1) :
    T_(T0, T1),
    sign_((T0 < T1) - (T1 < T0)),
    t_(T0)
{
}
//-----------------------------------------------------------------------------
Time::Time(Interval I) :
    T_(I),
    sign_((I.first < I.second) - (I.second < I.first)),
    t_(I.first)
{
}
//-----------------------------------------------------------------------------
Time::Time(Time const& other) :
    T_(other.T_),
    sign_((T_.first < T_.second) - (T_.second < T_.first)),
    t_(T_.first)
{
}
//-----------------------------------------------------------------------------
Time::~Time()
{
}
//-----------------------------------------------------------------------------
std::pair<real, real> const& Time::interval() const
{
  return T_;
}
//-----------------------------------------------------------------------------
int Time::sign() const
{
  return sign_;
}
//-----------------------------------------------------------------------------
bool Time::is_valid(real atol) const
{
  return std::fabs(t_ - T_.first)  <= std::fabs(T_.second - T_.first) + atol
      && std::fabs(t_ - T_.second) <= std::fabs(T_.second - T_.first) + atol;
}
//-----------------------------------------------------------------------------
real Time::begin() const
{
  return T_.first;
}
//-----------------------------------------------------------------------------
real Time::end() const
{
  return T_.second;
}
//-----------------------------------------------------------------------------
real Time::measure() const
{
  return std::fabs(T_.second - T_.first);
}
//-----------------------------------------------------------------------------
real const& Time::clock() const
{
  return t_;
}
//-----------------------------------------------------------------------------
real& Time::clock()
{
  return t_;
}
//-----------------------------------------------------------------------------
real Time::elapsed() const
{
  return std::fabs(t_ - T_.first);
}
//-----------------------------------------------------------------------------
real Time::remaining() const
{
  return std::fabs(T_.second - t_);
}
//-----------------------------------------------------------------------------
real Time::elapsed_normalized() const
{
  if(abscmp(T_.second, T_.first))
  {
    return 0.0;
  }
  return std::fabs(t_ - T_.first)/std::fabs(T_.second - T_.first);
}
//-----------------------------------------------------------------------------
real Time::remaining_normalized() const
{
  if(abscmp(T_.second, T_.first))
  {
    return 0.0;
  }
  return std::fabs(T_.second - t_)/std::fabs(T_.second - T_.first);
}
//-----------------------------------------------------------------------------
void Time::show() const
{
  real const p = 100.0 * this->elapsed() / (sign_ == 0 ? 1.0 : this->measure());
  message("----------------------------------------------------------------");
  message("t = %-49.6e [ %6.2f ]", t_, p);
  message("----------------------------------------------------------------");
}
//-----------------------------------------------------------------------------
void Time::disp() const
{
  dolfin::section("Time");
  message("T0        : %f", T_.first);
  message("T1        : %f", T_.second);
  message("t         : %f", t_);
  message("valid     : %d", this->is_valid());
  message("sign      : %d", this->sign());
  message("measure   : %f", this->measure());
  message("elapsed   : %f", this->elapsed());
  message("remaining : %f", this->remaining());
  dolfin::end();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
