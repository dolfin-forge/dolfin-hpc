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
  static DofNumbering * create(Mesh& mesh, ufc::dofmap& ufc_dofmap);

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
  virtual bool is_shared(uint index) const = 0;

  /// Return is the dof is ghosted
  /// No checking performed: return false if the index is not known (!)
  virtual bool is_ghost(uint index) const = 0;

  /// Return string description
  virtual std::string description() const = 0;

  //---------------------------------------------------------------------------

  /// Return local dof numbering offset
  uint offset() const;

  /// Return local dof numbering size
  uint size() const;

  /// Cached cell tabulated local-to-global mapping
  uint const * block() const;

  /// Cached cell tabulated local-to-global mapping
  uint block_size() const;

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
  DofNumbering& operator=(DofNumbering const& other);

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
inline DofNumbering & DofNumbering::operator=( DofNumbering const & )
{
  return *this;
}
//-----------------------------------------------------------------------------
inline uint DofNumbering::offset() const
{
  return offset_;
}
//-----------------------------------------------------------------------------
inline uint DofNumbering::size() const
{
  return size_;
}
//-----------------------------------------------------------------------------
inline uint const * DofNumbering::block() const
{
  if ( array == nullptr )
  {
    pretabulate( array, array_size );
  }
  return array;
}
//-----------------------------------------------------------------------------
inline uint DofNumbering::block_size() const
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
