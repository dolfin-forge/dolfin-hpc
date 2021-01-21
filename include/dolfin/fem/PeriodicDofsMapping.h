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
  auto max_local_dimension() const -> uint;

  ///
  auto num_Gdofs() const -> uint;
  auto num_Hdofs() const -> uint;
  auto num_Idofs() const -> uint;

  ///
  auto is_Gdof( uint i ) const -> bool;
  auto is_Hdof( uint i ) const -> bool;
  auto is_Idof( uint i ) const -> bool;

  ///
  auto get_Gindices() const -> uint const *;

  ///
  void tabulate_dofs( uint Gdof, uint * Hdofs, uint & count ) const;

  //
  void tabulate_dofs( uint i, uint * Gdof, uint * Hdofs, uint & count ) const;

  ///
  void tabulate_coordinates( uint    Gdof,
                             real *  Gcoords,
                             real ** Hcoords,
                             uint &  count ) const;

  //
  void tabulate_coordinates( uint    i,
                             uint *  Gdof,
                             real *  Gcoords,
                             real ** Hcoords,
                             uint &  count ) const;

  ///
  void disp() const;

private:
  ///
  void init();

  ///
  void clear();

  FiniteElementSpace const & space_;

  ///
  uint max_local_dimension_;

  /// Map of G dofs to the offset and count in the arrays
  using OffsetMap = _ordered_map< uint, uint >;
  OffsetMap    Goffsets_;
  uint *       Gindices_;
  real *       Gxcoords_;
  _set< uint > Hdofs_;
  uint *       Hcount_;
  uint *       Hoffsets_;
  uint *       Hindices_;
  real *       Hxcoords_;
  _set< uint > Idofs_;
};

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::max_local_dimension() const -> uint
{
  return max_local_dimension_;
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::num_Gdofs() const -> uint
{
  return Goffsets_.size();
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::num_Hdofs() const -> uint
{
  return Hdofs_.size();
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::num_Idofs() const -> uint
{
  return Idofs_.size();
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::is_Gdof( uint i ) const -> bool
{
  return ( Goffsets_.find( i ) != Goffsets_.end() );
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::is_Hdof( uint i ) const -> bool
{
  return ( Hdofs_.count( i ) > 0 );
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::is_Idof( uint i ) const -> bool
{
  return ( Idofs_.count( i ) > 0 );
}

//-----------------------------------------------------------------------------
inline auto PeriodicDofsMapping::get_Gindices() const -> uint const *
{
  return Gindices_;
}

//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_PERIODIC_DOFS_MAPPING */
