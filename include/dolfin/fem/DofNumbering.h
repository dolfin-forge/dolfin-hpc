// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_DOF_NUMBERING_H
#define __DOLFIN_DOF_NUMBERING_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>

#include <ufc.h>

namespace dolfin
{

class Cell;
class Mesh;
class UFCCell;

/**
 *  @class  DofNumbering
 *
 *  @brief  Implements an interface for dofs numbering schemes to be called in
 *          the DofMap class.
 *
 */

class DofNumbering
{

public:
  /// Create dof numbering for given ufc dofmap (factory function)
  static auto create( Mesh & mesh, ufc::dofmap & ufc_dofmap ) -> DofNumbering *;

  /// Default constructor given a UFC mesh and dofmap
  DofNumbering( Mesh & mesh, ufc::dofmap & ufc_dofmap );

  // Destructor
  virtual ~DofNumbering();

  //--- INTERFACE -------------------------------------------------------------

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs( size_t *          dofs,
                              ufc::cell const & ufc_cell,
                              Cell const &      cell ) const = 0;

  /// Build the dof numbering
  virtual void build() = 0;

  /// Return is the dof is shared
  /// No checking performed: return false if the index is not known (!)
  virtual auto is_shared( size_t index ) const -> bool = 0;

  /// Return is the dof is ghosted
  /// No checking performed: return false if the index is not known (!)
  virtual auto is_ghost( size_t index ) const -> bool = 0;

  /// Return string description
  virtual auto description() const -> std::string = 0;

  //---------------------------------------------------------------------------

  /// Return local dof numbering offset
  auto offset() const -> size_t;

  /// Return local dof numbering size
  auto size() const -> size_t;

  /// Cached cell tabulated local-to-global mapping
  auto block() const -> size_t const *;

  /// Cached cell tabulated local-to-global mapping
  auto block_size() const -> size_t;

  //---------------------------------------------------------------------------

  /// Compatibility function
  auto tabulate_dofs( size_t * dofs, UFCCell const & ufc_cell ) -> void;

  /// Display information
  auto disp() const -> void;

  /// Initialize UFC mesh and dofmap
  static auto init( Mesh & mesh, ufc::dofmap & ufc_dofmap ) -> void;

protected:
  /// Clear existing data and initialize UFC dofmap
  auto init() -> void;

  ///
  auto set_range( size_t offset, size_t size ) -> void;

  //--- PROTECTED ATTRIBUTES --------------------------------------------------

  //
  Mesh & mesh;

  //
  ufc::dofmap & ufc_dofmap;

  // Caching of cell tabulated dofs array
  mutable size_t   array_size;
  mutable size_t * array;

private:
  /// Copy constructor
  DofNumbering( DofNumbering const & other );

  /// Assignment
  auto operator=( DofNumbering const & other ) -> DofNumbering &;

  ///
  auto clear() -> void;

  /// Pretabulate the local-to-global mapping of dofs for all the cells
  auto pretabulate( size_t *& array, size_t & array_size ) const -> void;

  //
  size_t offset_;

  //
  size_t size_;
};

//-----------------------------------------------------------------------------

inline auto DofNumbering::operator=( DofNumbering const & ) -> DofNumbering &
{
  return *this;
}

//-----------------------------------------------------------------------------

inline auto DofNumbering::offset() const -> size_t
{
  return offset_;
}

//-----------------------------------------------------------------------------

inline auto DofNumbering::size() const -> size_t
{
  return size_;
}

//-----------------------------------------------------------------------------

inline auto DofNumbering::block() const -> size_t const *
{
  if ( array == nullptr )
  {
    pretabulate( array, array_size );
  }
  return array;
}

//-----------------------------------------------------------------------------

inline auto DofNumbering::block_size() const -> size_t
{
  if ( array == nullptr )
  {
    pretabulate( array, array_size );
  }
  return array_size;
}

//-----------------------------------------------------------------------------

inline auto DofNumbering::init() -> void
{
  DofNumbering::clear();
  DofNumbering::init( mesh, ufc_dofmap );
}

//-----------------------------------------------------------------------------

inline auto DofNumbering::clear() -> void
{
  offset_    = 0;
  size_      = 0;
  array_size = 0;
  delete[] array;
  array = nullptr;
}

//-----------------------------------------------------------------------------

inline auto DofNumbering::set_range( size_t offset, size_t size ) -> void
{
  offset_ = offset;
  size_   = size;
}

//-----------------------------------------------------------------------------

inline auto DofNumbering::tabulate_dofs( size_t *        dofs,
                                         UFCCell const & ufc_cell ) -> void
{
  this->tabulate_dofs( dofs, ufc_cell, *ufc_cell );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_DOF_NUMBERING_H */
