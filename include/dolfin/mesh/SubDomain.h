// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2013.
// Modified by Aurélien Larcher 2014.
//
// First added:  2007-04-10
// Last changed: 2014-02-07

#ifndef __DOLFIN_SUB_DOMAIN_H
#define __DOLFIN_SUB_DOMAIN_H

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

  ///
  SubDomain(Mesh& mesh);

  /// Destructor
  virtual ~SubDomain();

  /// Return true for points inside the sub domain
  virtual bool inside(real const * x, bool const on_boundary) const;

  /// Return true if all vertices are inside the subdomain
  bool inside(MeshEntity& entity, bool const on_boundary) const;

  /// Set sub domain markers for given sub domain
  void mark(MeshFunction<uint>& sub_domains, uint sub_domain) const;

protected:

  /// Return if the coordinate is close given internal tolerance
  bool close(real x, real const xref) const;

private:

  ///
  bool intersect(real const * x, uint const dim, bool const on_boundary) const;

  ///
  bool intersect(Point const& p, bool const on_boundary) const;

  /// Intersection detector
  mutable IntersectionDetector* intersection_detector;

protected:

  ///
  real BMARG;

};

}

#endif
