// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-05-23
// Last changed: 2008-05-23

#ifndef __DOLFIN_DOMAIN_BOUNDARY_H
#define __DOLFIN_DOMAIN_BOUNDARY_H

#include <dolfin/mesh/SubDomain.h>

namespace dolfin
{

/**
 *  @class  DomainBoundary
 *
 *  @brief  This class provides a SubDomain which picks out the boundary of a
 *          mesh, and provides a convenient way to specify boundary conditions
 *          on the entire boundary of a mesh.
 */

class DomainBoundary : public SubDomain
{

public:

  /// Constructor
  DomainBoundary() :
    SubDomain()
  {
  }

  /// Destructor
  ~DomainBoundary()
  {
  }

  /// Return true for points on the boundary
  bool inside(const real* x, bool on_boundary) const
  {
    return on_boundary;
  }

};

} /* namespace dolfin */

#endif /* __DOLFIN_DOMAIN_BOUNDARY_H */
