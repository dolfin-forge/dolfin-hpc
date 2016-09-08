// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-12-21
// Last changed: 2008-06-20

#include <dolfin/common/timing.h>
#include <dolfin/log/dolfin_log.h>
#include <ctime>

namespace dolfin
{

#define DOLFIN_MAX_TIMER 4096
clock_t _dolfin_timer_[DOLFIN_MAX_TIMER];
uint _dolfin_id_ = 0;

//-----------------------------------------------------------------------------
void tic()
{
  if((_dolfin_id_ + 1) == DOLFIN_MAX_TIMER)
  {
    error("timing : maximum of %u timers has been reached", DOLFIN_MAX_TIMER);
  }
  _dolfin_timer_[_dolfin_id_] = clock();
  ++_dolfin_id_;
}
//-----------------------------------------------------------------------------
real toc()
{
  if(_dolfin_id_ == 0)
  {
    error("timing : returning timer without any prior tic");
  }
  --_dolfin_id_;
  return ((real) (clock() - _dolfin_timer_[_dolfin_id_])) / CLOCKS_PER_SEC;
}
//-----------------------------------------------------------------------------
real tocd(uint level)
{
  real const elapsed_time = toc();
  message(level, "Elapsed time: %8e seconds", elapsed_time);
  return elapsed_time;
}
//-----------------------------------------------------------------------------
real time()
{
  return ((real) (clock())) / CLOCKS_PER_SEC;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
