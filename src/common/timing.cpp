// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/common/timing.h>
#include <dolfin/log/dolfin_log.h>

#include <ctime>
#include <stack>

namespace dolfin
{

std::stack< clock_t > _dolfin_timer_;

//-----------------------------------------------------------------------------
void tic()
{
  _dolfin_timer_.push( clock() );
}
//-----------------------------------------------------------------------------
auto toc() -> real
{
  if ( _dolfin_timer_.empty() )
  {
    error( "timing : returning timer without any prior tic" );
  }

  real t = static_cast<real>( clock() - _dolfin_timer_.top() ) / CLOCKS_PER_SEC;
  _dolfin_timer_.pop();
  return t;
}
//-----------------------------------------------------------------------------
auto tocd( uint level ) -> real
{
  real const elapsed_time = toc();
  message( level, "Elapsed time: %8e seconds", elapsed_time );
  return elapsed_time;
}
//-----------------------------------------------------------------------------
auto time() -> real
{
  return static_cast< real >( clock() ) / CLOCKS_PER_SEC;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
