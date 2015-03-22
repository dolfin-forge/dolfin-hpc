// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-04-10
// Last changed: 2014-02-07

#ifndef __PERIODIC_SUB_DOMAIN_H
#define __PERIODIC_SUB_DOMAIN_H

#include <dolfin/common/types.h>

#include "SubDomain.h"

namespace dolfin
{

/// This class defines the interface for definition of periodic sub domains.

class PeriodicSubDomain : public SubDomain
{
public:

  /// Constructor
  PeriodicSubDomain() :
      SubDomain()
  {
  }

  /// Destructor
  virtual ~PeriodicSubDomain()
  {
  }

  /// Return true for points inside the sub domain
  virtual bool inside(real const * x, bool const on_boundary) const = 0;

  /// Map coordinate x in domain H to coordinate y in domain G
  /// (used for periodic boundary conditions)
  virtual void map(real const * xH, real* xG) const = 0;

};

}

#endif
