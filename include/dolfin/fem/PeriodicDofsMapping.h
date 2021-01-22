// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PERIODIC_DOFS_MAPPING
#define __DOLFIN_PERIODIC_DOFS_MAPPING

#include <dolfin/common/types.h>

namespace dolfin
{

class FiniteElementSpace;
class MeshEntity;

class PeriodicDofsMapping
{
public:
  ///
  PeriodicDofsMapping( FiniteElementSpace const & space );

  ///
  ~PeriodicDofsMapping();

  ///
  auto max_local_dimension() const -> size_t;

  ///
  auto num_Gdofs() const -> size_t;
  auto num_Hdofs() const -> size_t;
  auto num_Idofs() const -> size_t;

  ///
  auto is_Gdof( size_t i ) const -> bool;
  auto is_Hdof( size_t i ) const -> bool;
  auto is_Idof( size_t i ) const -> bool;

  ///
  auto get_Gindices() const -> size_t const *;

  ///
  void tabulate_dofs( size_t Gdof, size_t * Hdofs, size_t & count ) const;

  //
  void tabulate_dofs( size_t   i,
                      size_t * Gdof,
                      size_t * Hdofs,
                      size_t & count ) const;

  ///
  void tabulate_coordinates( size_t   Gdof,
                             real *   Gcoords,
                             real **  Hcoords,
                             size_t & count ) const;

  //
  void tabulate_coordinates( size_t   i,
                             size_t * Gdof,
                             real *   Gcoords,
                             real **  Hcoords,
                             size_t & count ) const;

  ///
  void disp() const;

private:
  ///
  void init();

  ///
  void clear();

  FiniteElementSpace const & space_;

  ///
  size_t max_local_dimension_;

  /// Map of G dofs to the offset and count in the arrays
  using OffsetMap = _ordered_map< size_t, size_t >;
  OffsetMap      Goffsets_;
  size_t *       Gindices_;
  real *         Gxcoords_;
  _set< size_t > Hdofs_;
  size_t *       Hcount_;
  size_t *       Hoffsets_;
  size_t *       Hindices_;
  real *         Hxcoords_;
  _set< size_t > Idofs_;
};

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::max_local_dimension() const -> size_t
{
  return max_local_dimension_;
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::num_Gdofs() const -> size_t
{
  return Goffsets_.size();
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::num_Hdofs() const -> size_t
{
  return Hdofs_.size();
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::num_Idofs() const -> size_t
{
  return Idofs_.size();
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::is_Gdof( size_t i ) const -> bool
{
  return ( Goffsets_.find( i ) != Goffsets_.end() );
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::is_Hdof( size_t i ) const -> bool
{
  return ( Hdofs_.count( i ) > 0 );
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::is_Idof( size_t i ) const -> bool
{
  return ( Idofs_.count( i ) > 0 );
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::get_Gindices() const -> size_t const *
{
  return Gindices_;
}

//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_PERIODIC_DOFS_MAPPING */
