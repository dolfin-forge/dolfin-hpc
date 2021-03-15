// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_DOF_MAP_H
#define __DOLFIN_DOF_MAP_H

#include <dolfin/common/types.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/MeshDependent.h>

#include <ufc.h>

namespace dolfin
{

class Cell;
class Form;
class Mesh;
class SubSytem;
class FiniteElementSpace;

/// This class handles the mapping of degrees of freedom.
/// It wraps a ufc::dofmap on a specific mesh and provides
/// optional precomputation and reordering of dofs.

class DofMap : public MeshDependent
{

public:
  /// Create dof map on mesh for i-th coefficient of given form
  DofMap( Mesh & mesh, ufc::form const & form, size_t const i );

  /// Create dof map on mesh from UFC object
  explicit DofMap( Mesh & mesh, ufc::dofmap const & dofmap );

  /// Create dof map on a subspace
  DofMap( DofMap const & dofmap, size_t i );

  /// Create dof map on a subspace for given subsystem
  DofMap( DofMap const &                dofmap,
          std::vector< size_t > const & sub_system,
          size_t &                      offset );

  /// Copy constructor
  DofMap( DofMap const & copy );

  /// Move constructor
  DofMap( DofMap && move ) = delete;

  /// Destructor
  ~DofMap() override;

  /// assignement operators
  auto operator=( DofMap const & copy ) = delete;
  auto operator=( DofMap && move ) = delete;

  /// Check if the element definitions are identical
  auto operator==( DofMap const & other ) const -> bool;
  auto operator!=( DofMap const & other ) const -> bool;

  /// access ufc interface
  auto ufc() const -> ufc::dofmap const & { return *ufc_dofmap_; }

  //--- Instantiation using the dofmap cache

  /// Acquire dofmap from cache for i-th function of the form.
  static auto acquire( Mesh & mesh, Form const & form, size_t const i )
    -> DofMap &;

  /// Acquire dofmap from cache for the given UFC dofmap.
  static auto acquire( Mesh & mesh, ufc::dofmap const & dofmap, bool owner )
    -> DofMap &;

  /// Release a token for the given dofmap
  static void release( DofMap & dofmap );

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// Tabulate the local-to-global mapping of dofs on a cell
  /// FIXME
  void tabulate_dofs( size_t *          dofs,
                      ufc::cell const & ufc_cell,
                      Cell const &      cell ) const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  /// FIXME
  void tabulate_dofs( size_t * dofs, UFCCell const & ufc_cell ) const;

  /// Extract sub dof map
  auto create_sub_dofmap( std::vector< size_t > const & sub_system ) const
    -> ufc::dofmap *;

  /// Extract sub dof map and compute the offset local to the reference
  /// element
  auto create_sub_dofmap( std::vector< size_t > const & sub_system,
                          size_t & local_offset ) const -> ufc::dofmap *;

  static auto create_sub_dofmap( ufc::dofmap const &           dofmap,
                                 std::vector< size_t > const & sub_system,
                                 size_t & local_offset ) -> ufc::dofmap *;

  /// Get sub dof maps offset (for a mixed element)
  auto sub_dofmaps_dimensions() const -> std::vector< size_t > const &;

  /// Get sub dof maps offset (for a mixed element)
  auto sub_dofmaps_offsets() const -> std::vector< size_t > const &;

  /// Get list of scalar dofmaps ordered by entries
  auto flatten() const -> std::vector< ufc::dofmap const * > const &;

  /// Create flatten representation of given dofmap (append sub dofmaps)
  static void flatten( ufc::dofmap const *                  dofmap,
                       std::vector< ufc::dofmap const * > & stack );

  /// Return if the dofmap can be seen as a vector element dofmap
  auto is_vectorizable() const -> bool;

  /// Return if the list of dofmap can be seen as a vector element
  static auto can_vectorize( std::vector< ufc::dofmap const * > flattened )
    -> bool;

  /// Unique identifier
  auto hash() const -> std::string const &;

  /// Return the dofmap local size i.e the process range (only owned dofs)
  auto local_size() const -> size_t;

  //--- Management of local-to-global mapping

  /// Return local to global mapping
  auto dofsmapping() const -> size_t const *;

  /// Return the size of the local to global mapping
  auto dofsmapping_size() const -> size_t;

  //--- Management of periodic dofs mapping

  /// Return local to global mapping
  auto periodic_mapping( FiniteElementSpace const & space ) const
    -> PeriodicDofsMapping const &;

  //--- Dof ownership

  /// Return is the dof is shared, return false if the index is not known (!)
  auto is_shared( size_t index ) const -> bool;

  /// Return is the dof is ghosted, return false if the index is not known (!)
  auto is_ghost( size_t index ) const -> bool;

  /// Display mapping
  void disp() const;

  //---

  /// Returns the dofmap signature corresponding to a given finite element
  static auto make_signature( std::string const & finite_element )
    -> std::string const;

  /// Create unique string identifiers for dofmap from UFC dofmap
  static auto make_hash( Mesh & mesh, ufc::dofmap const & ufc_dofmap )
    -> std::string const;

  /// Check consistency of ghosted entities
  auto check( bool throw_error = false ) -> bool;

private:
  // UFC dof map
  ufc::dofmap * const ufc_dofmap_;

public:
  std::string const signature;
  size_t const tdim;
  size_t const global_dim;
  size_t const num_sub_dofmaps;
  size_t const num_global_support_dofs;
  size_t const num_element_support_dofs;
  size_t const num_element_dofs;
  size_t const num_facet_dofs;
  std::vector< size_t > const num_entity_dofs;
  std::vector< size_t > const num_entity_closure_dofs;

private:
  /// Build dof numbering
  void build();

  // Forward declaration of offset to be update at creation of ufc::dofmap
  size_t offset_;

  // Use type to allow optimization of dof map ordering
  DofNumbering * const numbering_;

  // Dofmap hash
  std::string const hash_;

  // Sub dof maps offsets
  std::vector< size_t > sub_dofmaps_dims_;

  // Sub dof maps offsets
  std::vector< size_t > sub_dofmaps_offs_;

  //
  std::vector< ufc::dofmap const * > flattened_;

  // Periodic dofs mapping
  mutable PeriodicDofsMapping * periodic_dofmap_;
};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------

inline auto DofMap::operator==( DofMap const & other ) const -> bool
{
  return ( this->hash() == other.hash() );
}

//-----------------------------------------------------------------------------

inline auto DofMap::operator!=( DofMap const & other ) const -> bool
{
  return !( *this == other );
}

//-----------------------------------------------------------------------------

inline auto DofMap::acquire( Mesh & mesh, Form const & form, size_t const i )
  -> DofMap &
{
  return DofMapCache::instance().acquire( mesh, form, i );
}

//-----------------------------------------------------------------------------

inline auto DofMap::acquire( Mesh & mesh, ufc::dofmap const & dofmap, bool owner )
  -> DofMap &
{
  return DofMapCache::instance().acquire( mesh, dofmap, owner );
}

//-----------------------------------------------------------------------------

inline void DofMap::release( DofMap & dofmap )
{
  return DofMapCache::instance().release( dofmap );
}

//-----------------------------------------------------------------------------

inline void DofMap::tabulate_dofs( size_t *          dofs,
                                   ufc::cell const & ufc_cell,
                                   Cell const &      cell ) const
{
  dolfin_assert( dofs != nullptr );
  numbering_->tabulate_dofs( dofs, ufc_cell, cell );
}

//-----------------------------------------------------------------------------

inline void DofMap::tabulate_dofs( size_t *        dofs,
                                   UFCCell const & ufc_cell ) const
{
  dolfin_assert( dofs != nullptr );
  numbering_->tabulate_dofs( dofs, ufc_cell, *ufc_cell );
}

//-----------------------------------------------------------------------------

inline auto DofMap::is_vectorizable() const -> bool
{
  return DofMap::can_vectorize( this->flatten() );
}

//-----------------------------------------------------------------------------

inline auto DofMap::sub_dofmaps_dimensions() const
  -> std::vector< size_t > const &
{
  return sub_dofmaps_dims_;
}

//-----------------------------------------------------------------------------

inline auto DofMap::sub_dofmaps_offsets() const -> std::vector< size_t > const &
{
  return sub_dofmaps_offs_;
}

//-----------------------------------------------------------------------------

inline auto DofMap::flatten() const -> std::vector< ufc::dofmap const * > const &
{
  return flattened_;
}

//-----------------------------------------------------------------------------

inline auto DofMap::local_size() const -> size_t
{
  return numbering_->size();
}

//-----------------------------------------------------------------------------

inline auto DofMap::dofsmapping() const -> size_t const *
{
  return numbering_->block();
}

//-----------------------------------------------------------------------------

inline auto DofMap::dofsmapping_size() const -> size_t
{
  return numbering_->block_size();
}

//-----------------------------------------------------------------------------

inline auto DofMap::is_shared( size_t index ) const -> bool
{
  return numbering_->is_shared( index );
}

//-----------------------------------------------------------------------------

inline auto DofMap::is_ghost( size_t index ) const -> bool
{
  return numbering_->is_ghost( index );
}

//-----------------------------------------------------------------------------

inline auto DofMap::hash() const -> std::string const &
{
  return hash_;
}

//-----------------------------------------------------------------------------

inline auto DofMap::make_signature( std::string const & finite_element )
  -> std::string const
{
  return "FFC dofmap for " + finite_element;
}

//-----------------------------------------------------------------------------

inline auto DofMap::make_hash( Mesh & mesh, ufc::dofmap const & ufc_dofmap )
  -> std::string const
{
  return std::string( ufc_dofmap.signature() ) + mesh.hash();
}

//-----------------------------------------------------------------------------

}

#endif
