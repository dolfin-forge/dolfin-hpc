// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SUB_DOMAIN_H
#define __DOLFIN_SUB_DOMAIN_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

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
  virtual auto inside( real const * x, bool on_boundary ) const -> bool = 0;

  //---------------------------------------------------------------------------

  /// Return true if all vertices of given entity are inside the subdomain
  template < class Entity >
  auto enclosed( Entity & entity, bool on_boundary ) const -> bool;

  /// Return true if one vertex of given entity is inside the subdomain
  template < class Entity >
  auto overlap( Entity & entity, bool on_boundary ) const -> bool;

  /// Set sub domain markers for given sub domain
  template < class Entity >
  void mark( MeshValues< size_t, Entity > & sub_domains, size_t index ) const;

  //---------------------------------------------------------------------------

  /// Set geometric absolute tolerance
  inline void set_tolerance( real abstol )
  {
    abstol_ = std::fabs( abstol );
  }

protected:
  /// Return if the coordinate is close given provided tolerance
  auto close( real const x, real const xref, real const abstol ) const -> bool;

  /// Return if the coordinate is close given internal tolerance
  auto close( real const x, real const xref ) const -> bool;

private:
  real abstol_ { 1.0e-6 };
};

//-----------------------------------------------------------------------------
template <>
inline auto SubDomain::enclosed( Vertex & entity, bool on_boundary ) const
  -> bool
{
  return inside( entity.x(), on_boundary );
}

//-----------------------------------------------------------------------------
template < class Entity >
inline auto SubDomain::enclosed( Entity & entity, bool on_boundary ) const
  -> bool
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
inline auto SubDomain::overlap( Vertex & entity, bool on_boundary ) const
  -> bool
{
  return inside( entity.x(), on_boundary );
}

//-----------------------------------------------------------------------------
template < class Entity >
inline auto SubDomain::overlap( Entity & entity, bool on_boundary ) const
  -> bool
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
inline auto SubDomain::close( real const x,
                              real const xref,
                              real const abstol ) const -> bool
{
  return ( std::fabs( x - xref ) < abstol );
}

//-----------------------------------------------------------------------------
inline auto SubDomain::close( real const x, real const xref ) const -> bool
{
  return ( std::fabs( x - xref ) < abstol_ );
}

} /* namespace dolfin */

#endif /* __DOLFIN_SUB_DOMAIN_H */
