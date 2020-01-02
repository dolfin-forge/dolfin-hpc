// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PARAMETER_SYSTEM_H
#define __DOLFIN_PARAMETER_SYSTEM_H

#include "ParameterList.h"

namespace dolfin
{

  /// This class holds a global database of parameters for DOLFIN,
  /// implemented as a set of (key, value) pairs. Supported value
  /// types are real, int, bool, and string.

  class ParameterSystem : public ParameterList
  {
  public:

    /// Singleton instance of global parameter database
    static ParameterSystem parameters;

  private:

    // Constructor
    ParameterSystem();
    
  };

}

#endif
