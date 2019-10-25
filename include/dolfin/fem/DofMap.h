// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_DOF_MAP_H
#define __DOLFIN_DOF_MAP_H

#include <dolfin/mesh/MeshDependent.h>

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

#include <ufc.h>

#include <map>

namespace dolfin
{

class Cell;
class DofNumbering;
class Form;
class PeriodicDofsMapping;
class Mesh;
class SubSytem;
class UFCCell;
class UFCMesh;

/// This class handles the mapping of degrees of freedom.
/// It wraps a ufc::dofmap on a specific mesh and provides
/// optional precomputation and reordering of dofs.

class DofMap : public ufc::dofmap, public MeshDependent
{

public:

  /// Create dof map on mesh for i-th coefficient of given form
  DofMap(Mesh& mesh, ufc::form const& form, uint const i);

  /// Create dof map on mesh from UFC object
  /// Ownership of the UFC object is transfered to the instance if the boolean
  /// is set to true, otherwise a clone of the dofmap is created.
  /// In any case the instance the member attribute will be destroyed.
  DofMap(Mesh& mesh, ufc::dofmap& dofmap, bool const owner);

  /// Create dof map on a subspace
  DofMap(DofMap const& dofmap, uint i);

  /// Create dof map on a subspace for given subsystem
  DofMap(DofMap const& dofmap, Array<uint> const& sub_system, uint& offset);

  /// Destructor
  ~DofMap();

  /// Check if the element definitions are identical
  bool operator ==(DofMap const& other) const;
  bool operator !=(DofMap const& other) const;

  //--- Instantiation using the dofmap cache

  /// Acquire dofmap from cache for i-th function of the form.
  static DofMap& acquire(Mesh& mesh, Form const& form, uint const i);

  /// Acquire dofmap from cache for the given UFC dofmap.
  static DofMap& acquire(Mesh& mesh, ufc::dofmap& dofmap, bool owner);

  /// Release a token for the given dofmap
  static void release(DofMap& dofmap);

  //--- INTERFACE -------------------------------------------------------------
  /// @remark Exposes a subset of UFC v2.0

  /// Return a string identifying the dof map
  /// UFC @since 1.1
  char const * signature() const;

  /// Return true iff mesh entities of topological dimension d are needed
  /// UFC @since 1.1
  bool needs_mesh_entities(uint d) const;

  /// Initialize dofmap for mesh (return true iff init_cell() is needed)
  /// UFC @since 1.1
  bool init_mesh(const ufc::mesh& mesh);

  /// Initialize dofmap for given cell
  /// UFC @since 1.1
  void init_cell(const ufc::mesh& m, const ufc::cell& c);

  /// Finish initialization of dofmap for cells
  /// UFC @since 1.1
  void init_cell_finalize();

  /// Return the topological dimension of the associated cell shape
  /// UFC @since 2.0
  uint topological_dimension() const;

  /// Return the geometric dimension of the coordinates this dof map provides
  /// UFC @since 1.1 but not implemented
  uint geometric_dimension() const;

  /// Return the dimension of the global finite element function space
  /// UFC @since 1.1
  uint global_dimension() const;

  /// Return the dimension of the local finite element function space
  /// UFC @since 1.1
  uint local_dimension() const;

  /// Return number of facet dofs
  /// UFC @since 1.1
  uint num_facet_dofs() const;

  /// Return the number of dofs associated with each cell entity of dimension d
  /// UFC @since 1.1 but not implemented
  uint num_entity_dofs(uint d) const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  /// UFC @since 1.1, specific cell local versions are called in priority
  void tabulate_dofs(uint* dofs, const ufc::mesh& m, const ufc::cell& c) const;

  /// Tabulate local-local facet dofs
  /// UFC @since 1.1
  void tabulate_facet_dofs(uint* dofs, uint local_facet) const;

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  /// UFC @since 1.1 but not implemented
  void tabulate_entity_dofs(uint* dofs, uint d, uint i) const;

  /// Tabulate the coordinates of all dofs on a cell
  /// UFC @since 1.1
  void tabulate_coordinates(real** coordinates,
                            const ufc::cell& ufc_cell) const;

  //// Return the number of sub dof maps (for a mixed element)
  /// UFC @since 1.1
  uint num_sub_dofmaps() const;

  /// Create a new dofmap for sub dof map i (for a mixed element)
  /// UFC @since 1.1
  ufc::dofmap* create_sub_dofmap(uint i) const;

  /// Create a new instance
  /// UFC @since 2.0
  ufc::dofmap* create() const;

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// Return the dimension of the local finite element function space
  uint macro_local_dimension() const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell, Cell const& cell) const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(uint* dofs, UFCCell const& ufc_cell) const;

  /// Extract sub dof map
  ufc::dofmap* create_sub_dofmap(Array<uint> const& sub_system) const;

  /// Extract sub dof map and compute the offset local to the reference element
  ufc::dofmap* create_sub_dofmap(Array<uint> const& sub_system,
                                 uint& local_offset) const;

  static ufc::dofmap* create_sub_dofmap(ufc::dofmap const& dofmap,
                                        Array<uint> const& sub_system,
                                        uint& local_offset);

  /// Get sub dof maps offset (for a mixed element)
  Array<uint> const& sub_dofmaps_dimensions() const;

  /// Get sub dof maps offset (for a mixed element)
  Array<uint> const& sub_dofmaps_offsets() const;

  /// Get list of scalar dofmaps ordered by entries
  Array<ufc::dofmap const *> const& flatten() const;

  /// Create flatten representation of given dofmap (append sub dofmaps)
  static void flatten(ufc::dofmap const * dofmap,
                      Array<ufc::dofmap const *>& stack, uint maxlevel);

  /// Create flatten representation of given dofmap (append sub dofmaps)
  static void flatten(ufc::dofmap const * dofmap,
                      Array<ufc::dofmap const *>& stack);

  /// Return if the dofmap can be seen as a vector element dofmap
  bool is_vectorizable() const;

  /// Return if the list of dofmap can be seen as a vector element
  static bool can_vectorize(Array<ufc::dofmap const *> flattened);

  /// Unique identifier
  std::string const& hash() const;

  /// Return the dofmap local size i.e the process range (only owned dofs)
  uint local_size() const;

  //--- Management of local-to-global mapping

  /// Return local to global mapping
  uint const * dofsmapping() const;

  /// Return the size of the local to global mapping
  uint dofsmapping_size() const;

  //--- Management of periodic dofs mapping

  /// Return local to global mapping
  PeriodicDofsMapping const& periodic_mapping() const;

  //--- Dof ownership

  /// Return is the dof is shared, return false if the index is not known (!)
  bool is_shared(uint index) const;

  /// Return is the dof is ghosted, return false if the index is not known (!)
  bool is_ghost(uint index) const;

  //--- Debugging

  /// Return renumbering (used for testing)
  std::map<uint, uint> get_map() const;

  /// Display mapping
  void disp() const;

  /// Return if the dof map has been renumbered
  bool renumbered() const;

  //---

  /// Returns the dofmap signature corresponding to a given finite element
  static std::string const make_signature(std::string const& finite_element);

  /// Create unique string identifiers for dofmap from UFC dofmap
  static std::string const make_hash(Mesh& mesh, ufc::dofmap const& ufc_dofmap);

  //--- Debugging

  /// Check consistency of ghosted entities
  bool check(bool throw_error = false);

private:

  /// Initialize
  void init();

  /// Build dof numbering
  void build();

  // Forward declaration of offset to be update at creation of ufc::dofmap
  uint offset_;

  // UFC dof map
  ufc::dofmap * const ufc_dofmap_;

  // Use type to allow optimization of dof map ordering
  DofNumbering * const numbering_;

  // Dofmap hash
  std::string const hash_;

  // Sub dof maps offsets
  Array<uint> sub_dofmaps_dims_;

  // Sub dof maps offsets
  Array<uint> sub_dofmaps_offs_;

  //
  mutable Array<ufc::dofmap const *> flattened_;

  //
  uint num_leaf_spaces_;

  // Periodic dofs mapping
  mutable PeriodicDofsMapping * periodic_dofmap_;

  // Provide easy access to map for testing
  std::map<uint, uint> map_;

};

//--- INLINES -----------------------------------------------------------------

inline char const * DofMap::signature() const
{
  return ufc_dofmap_->signature();
}

//-----------------------------------------------------------------------------
inline bool DofMap::needs_mesh_entities(uint d) const
{
  return ufc_dofmap_->needs_mesh_entities(d);
}

//-----------------------------------------------------------------------------
inline bool DofMap::init_mesh(const ufc::mesh& mesh)
{
  return ufc_dofmap_->init_mesh(mesh);
}

//-----------------------------------------------------------------------------
inline void DofMap::init_cell(const ufc::mesh& m, const ufc::cell& c)
{
  ufc_dofmap_->init_cell(m, c);
}

//-----------------------------------------------------------------------------
inline void DofMap::init_cell_finalize()
{
  ufc_dofmap_->init_cell_finalize();
}

//-----------------------------------------------------------------------------
inline uint DofMap::topological_dimension() const
{
  return ufc_dofmap_->topological_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::geometric_dimension() const
{
  return ufc_dofmap_->geometric_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::global_dimension() const
{
  return ufc_dofmap_->global_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::local_dimension() const
{
  return ufc_dofmap_->local_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::num_facet_dofs() const
{
  return ufc_dofmap_->num_facet_dofs();
}

//-----------------------------------------------------------------------------
inline uint DofMap::num_entity_dofs(uint d) const
{
  return ufc_dofmap_->num_entity_dofs(d);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_dofs(uint* dofs, const ufc::mesh& m,
                                  const ufc::cell& c) const
{
  ufc_dofmap_->tabulate_dofs(dofs, m, c);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_facet_dofs(uint* dofs, uint local_facet) const
{
  ufc_dofmap_->tabulate_facet_dofs(dofs, local_facet);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_entity_dofs(uint* dofs, uint d, uint i) const
{
  ufc_dofmap_->tabulate_entity_dofs(dofs, d, i);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_coordinates(real** coordinates,
                                         const ufc::cell& ufc_cell) const
{
  ufc_dofmap_->tabulate_coordinates(coordinates, ufc_cell);
}

//-----------------------------------------------------------------------------
inline uint DofMap::num_sub_dofmaps() const
{
  return ufc_dofmap_->num_sub_dofmaps();
}

//-----------------------------------------------------------------------------
inline ufc::dofmap* DofMap::create_sub_dofmap(uint i) const
{
  return ufc_dofmap_->create_sub_dofmap(i);
}

//-----------------------------------------------------------------------------
inline ufc::dofmap* DofMap::create() const
{
  return ufc_dofmap_->create();
}

//-----------------------------------------------------------------------------
inline uint DofMap::macro_local_dimension() const
{
  return 2 * ufc_dofmap_->local_dimension();
}

//-----------------------------------------------------------------------------

}

#endif
