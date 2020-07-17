// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SUB_DOMAIN_H
#define __DOLFIN_SUB_DOMAIN_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/VertexIterator.h>

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
  virtual ~SubDomain() = default;

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

  real abstol_{1.0e-6};

};

//-----------------------------------------------------------------------------
template <>
inline bool SubDomain::enclosed( Vertex & entity, bool on_boundary ) const
{
  return inside( entity.x(), on_boundary );
}

//-----------------------------------------------------------------------------
template < class Entity >
inline bool SubDomain::enclosed( Entity & entity, bool on_boundary ) const
{
  for ( VertexIterator v( entity ); !v.end(); ++v )
  {
    if ( !this->inside( v->x(), on_boundary ) )
    {
      return false;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
template <>
inline bool SubDomain::overlap( Vertex & entity, bool on_boundary ) const
{
  return inside( entity.x(), on_boundary );
}

//-----------------------------------------------------------------------------
template < class Entity >
inline bool SubDomain::overlap( Entity & entity, bool on_boundary ) const
{
  for ( VertexIterator v( entity ); !v.end(); ++v )
  {
    if ( this->inside( v->x(), on_boundary ) )
    {
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------------------
inline bool
  SubDomain::close( real const x, real const xref, real const abstol ) const
{
  return ( std::fabs( x - xref ) < abstol );
}

//-----------------------------------------------------------------------------
inline bool SubDomain::close( real const x, real const xref ) const
{
  return ( std::fabs( x - xref ) < abstol_ );
}

} /* namespace dolfin */

#endif /* __DOLFIN_SUB_DOMAIN_H */
