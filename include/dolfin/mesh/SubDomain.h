// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2013.
//
// First added:  2007-04-10
// Last changed: 2013-08-08

#ifndef __SUB_DOMAIN_H
#define __SUB_DOMAIN_H

#include <dolfin/common/types.h>
#include "MeshFunction.h"

namespace dolfin
{


  class IntersectionDetector;

  /// This class defines the interface for definition of sub domains.
  /// Alternatively, sub domains may be defined by a Mesh and a
  /// MeshFunction<uint> over the mesh.

  class SubDomain
  {
  public:

    /// Constructor
    SubDomain();

    SubDomain(Mesh& bmesh);

    /// Destructor
    virtual ~SubDomain();

    /// Return true for points inside the sub domain
    virtual bool inside(const real* x, bool on_boundary) const;

    /// Map coordinate x in domain H to coordinate y in domain G (used for periodic boundary conditions)
    virtual void map(const real* x, real* y) const;

    /// Set sub domain markers for given sub domain
    void mark(MeshFunction<uint>& sub_domains, uint sub_domain) const;

  private:

    bool intersect(real* x, uint dim, bool on_boundary) const;
    bool intersect(Point p, bool on_boundary) const;

    /// Intersection detector
    mutable IntersectionDetector* intersection_detector;

  };

}

#endif
