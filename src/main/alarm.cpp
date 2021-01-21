// Copyright (C) 2008-2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/main/alarm.h>

#include <dolfin/main/PE.h>
#include <dolfin/log/log.h>

#include <csignal>
#include <cstdio>
#include <sys/time.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
bool alarm::WALL_CLOCK_LIMIT = false;

//-----------------------------------------------------------------------------
void alarm::action(int sig_code)
{
  switch (sig_code)
    {
    case SIGALRM:
      WALL_CLOCK_LIMIT = true;
      if (PE::rank() == 0)
      {
        warning("alarm : wall clock limit reached");
      }
      itimerval itv;
      itv.it_value.tv_sec     = RESPITE;
      itv.it_value.tv_usec    = 0;
      itv.it_interval.tv_sec  = 0;
      itv.it_interval.tv_usec = 0;
      if (setitimer(ITIMER_REAL, &itv, nullptr) < 0)
      {
        perror("alarm : setitimer failed");
      }
      break;
    default:
      break;
    }
}
//-----------------------------------------------------------------------------
auto alarm::set_limit(long wall_clock_limit) -> bool
{
  if (wall_clock_limit > 0)
  {
    struct sigaction sig_param;
    sig_param.sa_handler = alarm::action;
    sigemptyset(&sig_param.sa_mask);
    sig_param.sa_flags = SA_RESTART;
    if (sigaction(SIGALRM, &sig_param, nullptr) < 0)
    {
      perror("alarm : sigaction failed");
    }

    itimerval itv;
    itv.it_value.tv_sec     = wall_clock_limit;
    itv.it_value.tv_usec    = 0;
    itv.it_interval.tv_sec  = 0;
    itv.it_interval.tv_usec = 0;
    if (setitimer(ITIMER_REAL, &itv, nullptr) < 0)
    {
      perror("alarm : setitimer failed");
    }
    if (PE::rank() == 0)
    {
      message("Wall clock limit set to %ld seconds", wall_clock_limit);
    }
    return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
auto alarm::state() const -> bool
{
  return WALL_CLOCK_LIMIT;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
