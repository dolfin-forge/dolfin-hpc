// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TIMER_H
#define __DOLFIN_TIMER_H

#include <dolfin/common/timing.h>

namespace dolfin
{

/// A timer can be used for timing tasks. The basic usage is
///
///   Timer timer("Assembling over cells");
///
/// The timer is started at construction and timing ends
/// when the timer is destroyed (goes out of scope). It is
/// also possible to start and stop a timer explicitly by
///
///   timer.start();
///   timer.stop();

class Timer
{
public:
  /// Create timer
  Timer( std::string const & task_ )
    : task( task_ )
    , t( time() )
    , stopped( false )
  {
  }

  /// Destructor
  ~Timer()
  {
    if ( not stopped )
      stop();
  }

  /// Start timer
  inline void start();

  /// Stop timer
  inline void stop();

private:
  // Name of task
  std::string task;

  // Start time
  real t;

  // True if timer has been stopped
  bool stopped;
};

//-----------------------------------------------------------------------------
inline void Timer::start()
{
  t       = time();
  stopped = false;
}
//-----------------------------------------------------------------------------
inline void Timer::stop()
{
  timing( task.c_str(), time() - t );
  stopped = true;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
