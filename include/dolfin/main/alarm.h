// Copyright (C) 2008-2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ALARM_H
#define __DOLFIN_ALARM_H

#include <dolfin/main/init.h>

namespace dolfin
{

struct alarm
{
  static long const RESPITE = 300;

  // Constructor
  alarm() = default;

  // Destructor
  ~alarm() = default;

  // Signal handler
  static void action(int sig_code);

  // Set wall clock limit
  bool set_limit(long wall_clock_limit);

  // Return if limit is being reached
  bool state() const;

private:

  //
  static bool WALL_CLOCK_LIMIT;

};

} /* namespace dolfin */

#endif /* __DOLFIN_ALARM_H */
