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

namespace dolfin
{

class Mesh;
class MeshEntity;
template<class T> class MeshFunction;

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
  virtual bool inside(real const * x, bool const on_boundary) const = 0;

  //---------------------------------------------------------------------------

  /// Return true if all vertices of given entity are inside the subdomain
  virtual bool inside(MeshEntity& entity, bool const on_boundary) const;

  /// Return true if one vertex of given entity is inside the subdomain
  virtual bool overlap(MeshEntity& entity, bool const on_boundary) const;

  /// Set sub domain markers for given sub domain
  virtual void mark(MeshFunction<uint>& sub_domains, uint index) const;

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
