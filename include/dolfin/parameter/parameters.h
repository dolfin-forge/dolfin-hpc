// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PARAMETERS_H
#define __DOLFIN_PARAMETERS_H

#include "Parameter.h"
#include <ostream>

namespace dolfin
{

  /// Get value of parameter with given key
  Parameter dolfin_get(std::string key);
  
  /// Set value of parameter
  void dolfin_set(std::string key, dolfin::Parameter value);

  /// Set special key/value pair
  void dolfin_set(std::string key, std::ostream& ostream);

  /// Add parameter
  void dolfin_add(std::string key, dolfin::Parameter value);
 
}

#endif
