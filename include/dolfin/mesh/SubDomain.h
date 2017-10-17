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
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

class Mesh;
class MeshEntity;

/// This class defines the interface for definition of sub domains.

class SubDomain
{

public:

  /// Constructor
  SubDomain();

  /// Destructor
  virtual ~SubDomain();

  //--- INTERFACE -------------------------------------------------------------

  /// Return true for points inside the sub domain
  virtual bool inside(real const * x, bool on_boundary) const = 0;

  //---------------------------------------------------------------------------

  /// Return true if all vertices of given entity are inside the subdomain
  template <class Entity>
  bool enclosed(Entity& entity, bool on_boundary) const;

  /// Return true if one vertex of given entity is inside the subdomain
  template <class Entity>
  bool overlap(Entity& entity, bool on_boundary) const;

  /// Set sub domain markers for given sub domain
  template <class Entity>
  void mark(MeshValues<uint, Entity>& sub_domains, uint index) const;

  //---------------------------------------------------------------------------

  /// Set geometric absolute tolerance
  inline void set_tolerance(real abstol) { abstol_ = std::fabs(abstol); }

protected:

  /// Return if the coordinate is close given provided tolerance
  bool close(real const x, real const xref, real const abstol) const;

  /// Return if the coordinate is close given internal tolerance
  bool close(real const x, real const xref) const;

private:

  real abstol_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_SUB_DOMAIN_H */
