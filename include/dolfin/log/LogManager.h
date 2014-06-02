// Copyright (C) 2003-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Thanks to Jim Tilander for many helpful hints.
// Modified by Aurélien Larcher, 2014.
//
// First added:  2003-03-13
// Last changed: 2014-06-01

#ifndef __LOG_MANAGER_H
#define __LOG_MANAGER_H

#include "Logger.h"

namespace dolfin
{

class LogManager
{

public:

  // Meyers singleton
  static Logger& logger()
  {
    static Logger logger_;
    return logger_;
  }

private:

  LogManager();
  ~LogManager();

};

}

#endif
