// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_DOF_MAP_H
#define __DOLFIN_DOF_MAP_H

#include <dolfin/common/Array.h>
#include <dolfin/common/types.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/MeshDependent.h>
#include <dolfin/ufc/ufc.h>

namespace dolfin
{

class Cell;
class Form;
class Mesh;
class SubSytem;

/// This class handles the mapping of degrees of freedom.
/// It wraps a ufc::dofmap on a specific mesh and provides
/// optional precomputation and reordering of dofs.

class DofMap : public ufc::dofmap, public MeshDependent
{

public:
  /// Create dof map on mesh for i-th coefficient of given form
  DofMap( Mesh & mesh, ufc::form const & form, size_t const i );

  /// Create dof map on mesh from UFC object
  /// Ownership of the UFC object is transfered to the instance if the boolean
  /// is set to true, otherwise a clone of the dofmap is created.
  /// In any case the instance the member attribute will be destroyed.
  DofMap( Mesh & mesh, ufc::dofmap & dofmap, bool const owner );

  /// Create dof map on a subspace
  DofMap( DofMap const & dofmap, size_t i );

  /// Create dof map on a subspace for given subsystem
  DofMap( DofMap const &          dofmap,
          Array< size_t > const & sub_system,
          size_t &                offset );

  /// Destructor
  ~DofMap() override;

  /// Check if the element definitions are identical
  bool operator==( DofMap const & other ) const;
  bool operator!=( DofMap const & other ) const;

  //--- Instantiation using the dofmap cache

  /// Acquire dofmap from cache for i-th function of the form.
  static DofMap & acquire( Mesh & mesh, Form const & form, size_t const i );

  /// Acquire dofmap from cache for the given UFC dofmap.
  static DofMap & acquire( Mesh & mesh, ufc::dofmap & dofmap, bool owner );

  /// Release a token for the given dofmap
  static void release( DofMap & dofmap );

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Return a string identifying the dof map
  char const * signature() const override;

  /// Return true iff mesh entities of topological dimension d are needed
  bool needs_mesh_entities( size_t d ) const override;

  /// Return the topological dimension of the associated cell shape
  size_t topological_dimension() const override;

  /// Return the dimension of the global finite element function space
  size_t global_dimension(
    std::vector< size_t > const & num_global_mesh_entities ) const override;

  /// Return the dimension of the local finite element function space
  /// Return the number of dofs with global support (i.e. global constants)
  size_t num_global_support_dofs() const override;

  /// Return the dimension of the local finite element function space
  /// for a cell (not including global support dofs)
  /// was local_dimension() before
  size_t num_element_support_dofs() const override;

  /// Return the dimension of the local finite element function space
  /// for a cell (old version including global support dofs)
  size_t num_element_dofs() const override;

  /// Return number of facet dofs
  size_t num_facet_dofs() const override;

  /// Return the number of dofs associated with each cell entity of dimension d
  size_t num_entity_dofs( size_t d ) const override;

  /// Return the number of dofs associated with the closure
  /// of each cell entity dimension d
  size_t num_entity_closure_dofs( size_t d ) const override;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs( size_t *                      dofs,
                      std::vector< size_t > const & num_global_entities,
                      std::vector< std::vector< size_t > > const &
                        entity_indices ) const override;

  /// Tabulate local-local facet dofs
  void tabulate_facet_dofs( size_t * dofs, size_t local_facet ) const override;

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  void tabulate_entity_dofs( size_t * dofs, size_t d, size_t i ) const override;

  // /// Tabulate the coordinates of all dofs on a cell
  // void tabulate_coordinates(real** coordinates,
  //                           const ufc::cell& ufc_cell) const override;

  /// Tabulate the local-to-local mapping of dofs on the closure of entity (d, i)
  void tabulate_entity_closure_dofs( size_t * dofs,
                                     size_t   d,
                                     size_t   i ) const override;

  //// Return the number of sub dof maps (for a mixed element)
  size_t num_sub_dofmaps() const override;

  /// Create a new dofmap for sub dof map i (for a mixed element)
  ufc::dofmap * create_sub_dofmap( size_t i ) const override;

  /// Create a new instance
  ufc::dofmap * create() const override;

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// FIXME
  // /// Return the dimension of the local finite element function space
  // size_t macro_local_dimension() const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  /// FIXME
  void tabulate_dofs( size_t *          dofs,
                      ufc::cell const & ufc_cell,
                      Cell const &      cell ) const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  /// FIXME
  void tabulate_dofs( size_t * dofs, UFCCell const & ufc_cell ) const;

  /// Extract sub dof map
  ufc::dofmap * create_sub_dofmap( Array< size_t > const & sub_system ) const;

  /// Extract sub dof map and compute the offset local to the reference element
  ufc::dofmap * create_sub_dofmap( Array< size_t > const & sub_system,
                                   size_t &                local_offset ) const;

  static ufc::dofmap * create_sub_dofmap( ufc::dofmap const &     dofmap,
                                          Array< size_t > const & sub_system,
                                          size_t & local_offset );

  /// Get sub dof maps offset (for a mixed element)
  Array< size_t > const & sub_dofmaps_dimensions() const;

  /// Get sub dof maps offset (for a mixed element)
  Array< size_t > const & sub_dofmaps_offsets() const;

  /// Get list of scalar dofmaps ordered by entries
  Array< ufc::dofmap const * > const & flatten() const;

  /// Create flatten representation of given dofmap (append sub dofmaps)
  static void flatten( ufc::dofmap const *            dofmap,
                       Array< ufc::dofmap const * > & stack,
                       size_t                         maxlevel );

  /// Create flatten representation of given dofmap (append sub dofmaps)
  static void flatten( ufc::dofmap const *            dofmap,
                       Array< ufc::dofmap const * > & stack );

  /// Return if the dofmap can be seen as a vector element dofmap
  bool is_vectorizable() const;

  /// Return if the list of dofmap can be seen as a vector element
  static bool can_vectorize( Array< ufc::dofmap const * > flattened );

  /// Unique identifier
  std::string const & hash() const;

  /// Return the dofmap local size i.e the process range (only owned dofs)
  size_t local_size() const;

  //--- Management of local-to-global mapping

  /// Return local to global mapping
  size_t const * dofsmapping() const;

  /// Return the size of the local to global mapping
  size_t dofsmapping_size() const;

  //--- Management of periodic dofs mapping

  /// Return local to global mapping
  PeriodicDofsMapping const & periodic_mapping() const;

  //--- Dof ownership

  /// Return is the dof is shared, return false if the index is not known (!)
  bool is_shared( size_t index ) const;

  /// Return is the dof is ghosted, return false if the index is not known (!)
  bool is_ghost( size_t index ) const;

  //--- Debugging

  /// Return renumbering (used for testing)
  _ordered_map< size_t, size_t > get_map() const;

  /// Display mapping
  void disp() const;

  /// Return if the dof map has been renumbered
  bool renumbered() const;

  //---

  /// Returns the dofmap signature corresponding to a given finite element
  static std::string const make_signature( std::string const & finite_element );

  /// Create unique string identifiers for dofmap from UFC dofmap
  static std::string const make_hash( Mesh &              mesh,
                                      ufc::dofmap const & ufc_dofmap );

  //--- Debugging

  /// Check consistency of ghosted entities
  bool check( bool throw_error = false );

private:
  /// Initialize
  void init();

  /// Build dof numbering
  void build();

  // Forward declaration of offset to be update at creation of ufc::dofmap
  size_t offset_;

  // UFC dof map
  ufc::dofmap * const ufc_dofmap_;

  // Use type to allow optimization of dof map ordering
  DofNumbering * const numbering_;

  // Dofmap hash
  std::string const hash_;

  // Sub dof maps offsets
  Array< size_t > sub_dofmaps_dims_;

  // Sub dof maps offsets
  Array< size_t > sub_dofmaps_offs_;

  //
  mutable Array< ufc::dofmap const * > flattened_;

  //
  size_t num_leaf_spaces_;

  // Periodic dofs mapping
  mutable PeriodicDofsMapping * periodic_dofmap_;

  // Provide easy access to map for testing
  _ordered_map< size_t, size_t > map_;
};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------

inline bool DofMap::operator==( DofMap const & other ) const
{
  return ( this->hash() == other.hash() );
}

//-----------------------------------------------------------------------------

inline bool DofMap::operator!=( DofMap const & other ) const
{
  return !( *this == other );
}

//-----------------------------------------------------------------------------

inline DofMap & DofMap::acquire( Mesh & mesh, Form const & form, size_t const i )
{
  return DofMapCache::instance().acquire( mesh, form, i );
}

//-----------------------------------------------------------------------------

inline DofMap & DofMap::acquire( Mesh & mesh, ufc::dofmap & dofmap, bool owner )
{
  return DofMapCache::instance().acquire( mesh, dofmap, owner );
}

//-----------------------------------------------------------------------------

inline void DofMap::release( DofMap & dofmap )
{
  return DofMapCache::instance().release( dofmap );
}

//-----------------------------------------------------------------------------

inline char const * DofMap::signature() const
{
  return ufc_dofmap_->signature();
}

//-----------------------------------------------------------------------------

inline bool DofMap::needs_mesh_entities( size_t d ) const
{
  return ufc_dofmap_->needs_mesh_entities( d );
}

//-----------------------------------------------------------------------------

inline size_t DofMap::topological_dimension() const
{
  return ufc_dofmap_->topological_dimension();
}

//-----------------------------------------------------------------------------

inline size_t DofMap::global_dimension(
  std::vector< size_t > const & num_global_mesh_entities ) const
{
  return ufc_dofmap_->global_dimension( num_global_mesh_entities );
}

//-----------------------------------------------------------------------------

inline size_t DofMap::num_global_support_dofs() const
{
  return ufc_dofmap_->num_global_support_dofs();
}

//-----------------------------------------------------------------------------

inline size_t DofMap::num_element_support_dofs() const
{
  return ufc_dofmap_->num_element_support_dofs();
}

//-----------------------------------------------------------------------------

inline size_t DofMap::num_element_dofs() const
{
  return ufc_dofmap_->num_element_dofs();
}

//-----------------------------------------------------------------------------

inline size_t DofMap::num_facet_dofs() const
{
  return ufc_dofmap_->num_facet_dofs();
}

//-----------------------------------------------------------------------------

inline size_t DofMap::num_entity_dofs( size_t d ) const
{
  return ufc_dofmap_->num_entity_dofs( d );
}

//-----------------------------------------------------------------------------

inline size_t DofMap::num_entity_closure_dofs( size_t d ) const
{
  return ufc_dofmap_->num_entity_closure_dofs( d );
}

//-----------------------------------------------------------------------------

inline void DofMap::tabulate_dofs(
  size_t *                                     dofs,
  std::vector< size_t > const &                num_global_entities,
  std::vector< std::vector< size_t > > const & entity_indices ) const
{
  ufc_dofmap_->tabulate_dofs( dofs, num_global_entities, entity_indices );
}

//-----------------------------------------------------------------------------

inline void DofMap::tabulate_facet_dofs( size_t * dofs,
                                         size_t   local_facet ) const
{
  ufc_dofmap_->tabulate_facet_dofs( dofs, local_facet );
}

//-----------------------------------------------------------------------------

inline void
  DofMap::tabulate_entity_dofs( size_t * dofs, size_t d, size_t i ) const
{
  ufc_dofmap_->tabulate_entity_dofs( dofs, d, i );
}

//-----------------------------------------------------------------------------

// inline void DofMap::tabulate_coordinates( real **           coordinates,
//                                           const ufc::cell & ufc_cell ) const
// {
//   ufc_dofmap_->tabulate_coordinates( coordinates, ufc_cell );
// }

inline void DofMap::tabulate_entity_closure_dofs( size_t * dofs,
                                                  size_t   d,
                                                  size_t   i ) const
{
  ufc_dofmap_->tabulate_entity_closure_dofs( dofs, d, i );
}

//-----------------------------------------------------------------------------

inline size_t DofMap::num_sub_dofmaps() const
{
  return ufc_dofmap_->num_sub_dofmaps();
}

//-----------------------------------------------------------------------------

inline ufc::dofmap * DofMap::create_sub_dofmap( size_t i ) const
{
  return ufc_dofmap_->create_sub_dofmap( i );
}

//-----------------------------------------------------------------------------

inline ufc::dofmap * DofMap::create() const
{
  return ufc_dofmap_->create();
}

//-----------------------------------------------------------------------------

// inline size_t DofMap::macro_local_dimension() const
// {
//   return 2 * ufc_dofmap_->local_dimension();
// }

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

inline bool DofMap::is_vectorizable() const
{
  return DofMap::can_vectorize( this->flatten() );
}

//-----------------------------------------------------------------------------

inline Array< size_t > const & DofMap::sub_dofmaps_dimensions() const
{
  return sub_dofmaps_dims_;
}

//-----------------------------------------------------------------------------

inline Array< size_t > const & DofMap::sub_dofmaps_offsets() const
{
  return sub_dofmaps_offs_;
}
//-----------------------------------------------------------------------------

inline size_t DofMap::local_size() const
{
  return numbering_->size();
}

//-----------------------------------------------------------------------------

inline size_t const * DofMap::dofsmapping() const
{
  return numbering_->block();
}

//-----------------------------------------------------------------------------

inline size_t DofMap::dofsmapping_size() const
{
  return numbering_->block_size();
}

//-----------------------------------------------------------------------------

inline PeriodicDofsMapping const & DofMap::periodic_mapping() const
{
  if ( periodic_dofmap_ == nullptr )
  {
    periodic_dofmap_ = new PeriodicDofsMapping( *this );
  }
  return *periodic_dofmap_;
}

//-----------------------------------------------------------------------------

inline _ordered_map< size_t, size_t > DofMap::get_map() const
{
  return map_;
}

//-----------------------------------------------------------------------------

inline bool DofMap::is_shared( size_t index ) const
{
  return numbering_->is_shared( index );
}

//-----------------------------------------------------------------------------

inline bool DofMap::is_ghost( size_t index ) const
{
  return numbering_->is_ghost( index );
}

//-----------------------------------------------------------------------------

inline std::string const & DofMap::hash() const
{
  return hash_;
}

//-----------------------------------------------------------------------------

inline std::string const
  DofMap::make_signature( std::string const & finite_element )
{
  return "FFC dofmap for " + finite_element;
}

//-----------------------------------------------------------------------------

inline std::string const DofMap::make_hash( Mesh &              mesh,
                                            ufc::dofmap const & ufc_dofmap )
{
  return std::string( ufc_dofmap.signature() ) + mesh.hash();
}

//-----------------------------------------------------------------------------

}

#endif
