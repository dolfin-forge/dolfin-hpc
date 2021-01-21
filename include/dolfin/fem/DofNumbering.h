// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_DOF_NUMBERING_H
#define __DOLFIN_DOF_NUMBERING_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/ufc/ufc.h>

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
  static auto create(Mesh& mesh, ufc::dofmap& ufc_dofmap) -> DofNumbering *;

  /// Default constructor given a UFC mesh and dofmap
  DofNumbering(Mesh& mesh, ufc::dofmap& ufc_dofmap);

  // Destructor
  virtual ~DofNumbering();

  //--- INTERFACE -------------------------------------------------------------

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell, Cell const& cell) const = 0;

  /// Build the dof numbering
  virtual void build() = 0;

  /// Return is the dof is shared
  /// No checking performed: return false if the index is not known (!)
  virtual auto is_shared(uint index) const -> bool = 0;

  /// Return is the dof is ghosted
  /// No checking performed: return false if the index is not known (!)
  virtual auto is_ghost(uint index) const -> bool = 0;

  /// Return string description
  virtual auto description() const -> std::string = 0;

  //---------------------------------------------------------------------------

  /// Return local dof numbering offset
  auto offset() const -> uint;

  /// Return local dof numbering size
  auto size() const -> uint;

  /// Cached cell tabulated local-to-global mapping
  auto block() const -> uint const *;

  /// Cached cell tabulated local-to-global mapping
  auto block_size() const -> uint;

  //---------------------------------------------------------------------------

  /// Compatibility function
  void tabulate_dofs(uint* dofs, UFCCell const& ufc_cell);

  /// Display information
  void disp() const;

  /// Initialize UFC mesh and dofmap
  static void init(Mesh& mesh, ufc::dofmap& ufc_dofmap);

protected:

  /// Clear existing data and initialize UFC dofmap
  void init();

  ///
  void set_range(uint offset, uint size);

  //--- PROTECTED ATTRIBUTES --------------------------------------------------

  //
  Mesh& mesh;

  //
  ufc::dofmap& ufc_dofmap;

  // Caching of cell tabulated dofs array
  mutable uint array_size;
  mutable uint * array;


private:

  /// Copy constructor
  DofNumbering(DofNumbering const& other);

  /// Assignment
  auto operator=(DofNumbering const& other) -> DofNumbering&;

  ///
  void clear();

  /// Pretabulate the local-to-global mapping of dofs for all the cells
  void pretabulate(uint *& array, uint& array_size) const;

  //
  uint offset_;

  //
  uint size_;

};

//-----------------------------------------------------------------------------
inline auto DofNumbering::operator=( DofNumbering const & ) -> DofNumbering &
{
  return *this;
}
//-----------------------------------------------------------------------------
inline auto DofNumbering::offset() const -> uint
{
  return offset_;
}
//-----------------------------------------------------------------------------
inline auto DofNumbering::size() const -> uint
{
  return size_;
}
//-----------------------------------------------------------------------------
inline auto DofNumbering::block() const -> uint const *
{
  if ( array == nullptr )
  {
    pretabulate( array, array_size );
  }
  return array;
}
//-----------------------------------------------------------------------------
inline auto DofNumbering::block_size() const -> uint
{
  if ( array == nullptr )
  {
    pretabulate( array, array_size );
  }
  return array_size;
}
//-----------------------------------------------------------------------------
inline void DofNumbering::init()
{
  DofNumbering::clear();
  DofNumbering::init( mesh, ufc_dofmap );
}
//-----------------------------------------------------------------------------
inline void DofNumbering::clear()
{
  offset_    = 0;
  size_      = 0;
  array_size = 0;
  delete[] array;
  array = nullptr;
}
//-----------------------------------------------------------------------------
inline void DofNumbering::set_range( uint offset, uint size )
{
  offset_ = offset;
  size_   = size;
}
//-----------------------------------------------------------------------------
inline void DofNumbering::tabulate_dofs( uint * dofs, UFCCell const & ufc_cell )
{
  this->tabulate_dofs( dofs, ufc_cell, *ufc_cell );
}

} /* namespace dolfin */

#endif /* __DOLFIN_DOF_NUMBERING_H */
