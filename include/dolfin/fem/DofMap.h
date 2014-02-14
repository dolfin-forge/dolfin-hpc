// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Martin Alnes, 2008
// Modified by Aurélien Larcher, 2013 (extension and partial rewrite)
//
// First added:  2007-03-01
// Last changed: 2008-04-10

#ifndef __DOF_MAP_H
#define __DOF_MAP_H

#include "UFCCell.h"
#include "UFCMesh.h"

#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>

#include <ufc.h>

#include <map>

namespace dolfin
{

class SubSytem;
class UFC;

/// This class handles the mapping of degrees of freedom.
/// It wraps a ufc::dof_map on a specific mesh and provides
/// optional precomputation and reordering of dofs.

class DofMap : protected ufc::dof_map
{

  static std::string const SIGN_PREFIX;

public:

  /// Returns the dofmap signature corresponding to a given finite element
  static std::string const dofmap_signature(std::string const& fe_signature);

  /// Returns the finite element signature corresponding to a given dofmap
  static std::string const finite_element_signature(
      std::string const& dofmap_signature);

  /// Create unique string identifiers for dofmap from signature
  static std::string const make_hash(std::string const& dofmap_signature,
                                     Mesh& mesh);

  /// Create unique string identifiers for dofmap from UFC dofmap
  static std::string const make_hash(ufc::dof_map& ufc_dofmap, Mesh& mesh);

  /// Create dof map on mesh
  DofMap(ufc::dof_map& dof_map, Mesh& mesh, bool const dof_map_local = false);

  /// Create dof map on mesh (parallel)
  DofMap(ufc::dof_map& dof_map, Mesh& mesh, MeshFunction<uint>& partitions,
         bool const dof_map_local = false);

  /// Create dof map on mesh
  DofMap(ufc::form const& form, uint const& i, Mesh& mesh);

  /// Create dof map on mesh (parallel)
  DofMap(ufc::form const& form, uint const& i, Mesh& mesh,
         MeshFunction<uint>& partitions);

  /// Create dof map on mesh
  DofMap(std::string const& signature, Mesh& mesh);

  /// Create dof map on mesh (parallel)
  DofMap(std::string const& signature, Mesh& mesh,
         MeshFunction<uint>& partitions);

  /// Create dof map on a sub space
  DofMap(DofMap const& dofmap, Array<uint> const& sub_system, uint& offset);

  /// Destructor
  ~DofMap();

  //--- INTERFACE -------------------------------------------------------------
  /// @remark Exposes a subset of UFC v1.1

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

  /// Return the dimension of the global finite element function space
  /// UFC @since 1.1
  uint global_dimension() const;

  /// Return the dimension of the local finite element function space
  /// UFC @since 1.1
  uint local_dimension() const;

  /// Return the geometric dimension of the coordinates this dof map provides
  /// UFC @since 1.1 but not implemented
  uint geometric_dimension() const;

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
  uint num_sub_dof_maps() const;

  /// Create a new dof_map for sub dof map i (for a mixed element)
  /// UFC @since 1.1
  ufc::dof_map* create_sub_dof_map(uint i) const;

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// Return the dimension of the local finite element function space
  uint macro_local_dimension() const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(uint* dofs, ufc::cell& ufc_cell, uint cell_index) const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(uint* dofs, const ufc::cell& ufc_cell,
                     uint cell_index) const;

  // FIXME: Can this function eventually be removed?
  /// Tabulate the local-to-global mapping of dofs on a ufc cell
  void tabulate_dofs(uint* dofs, const ufc::cell& cell) const;

  /// Get sub dof maps offset (for a mixed element)
  Array<uint> const& sub_dof_maps_dimensions() const;

  /// Get sub dof maps offset (for a mixed element)
  Array<uint> const& sub_dof_maps_offsets() const;

  /// Extract sub dof map
  ufc::dof_map* create_sub_dof_map(Array<uint> const& sub_system,
                                   uint& offset) const;

  /// Extract sub DofMap
  ufc::dof_map* create_sub_dof_map(ufc::dof_map const& dof_map,
                                   Array<uint> const& sub_system,
                                   uint& offset) const;

  /// Unique identifier
  std::string const& hash() const;

  /// Return the dofmap size on the mesh partition
  uint local_size() const;

  //--- Management of local-to-global mapping

  /// Return local to global mapping
  uint const * dofsmapping() const;

  /// Return the size of the local to global mapping
  uint dofsmapping_size() const;

  /// Return local cell mapping from component grouped numbering to node grouped
  uint const * cellmapping() const;

  /// Return the size of the local cell mapping
  uint cellmapping_size() const;

  //--- Mesh information

  /// Return mesh associated with map
  Mesh& mesh() const;

  ///
  std::string const& mesh_hash() const;

  //--- Debugging

  /// Return renumbering (used for testing)
  std::map<uint, uint> getMap() const;

  /// Display mapping
  void disp() const;

  /// Return if the dof map has been renumbered
  bool renumbered() const;

private:

  /// Initialise DofMap
  void init();

  /// Build parallel dof map
  void build();

  ///
  void pretabulate_all_dofs() const;

  ///
  void pretabulate_cell_dofs() const;

  // Forward declaration of internals for parallelism to ensure initialization
  int _type_;
  uint _offset_;
  uint _local_size_;
  uint * _v_map_;
  uint * _dof_map_; // Parallel dof map

  // Is UFC dof map lifetime managed by the DofMap instance
  bool const ufc_dof_map_local_;

  // UFC dof map
  ufc::dof_map * const ufc_dof_map_;

  // DOLFIN mesh
  Mesh& dolfin_mesh;

  // Partitions
  MeshFunction<uint>* partitions;

  // Number of cells in the mesh
  uint num_cells_;

  // UFC mesh
  UFCMesh ufc_mesh;

  // Hash
  std::string const mesh_hash_;

  // Dofmap hash
  std::string const hash_;

  // Sub dof maps offsets
  Array<uint> sub_dof_maps_dims_;

  // Sub dof maps offsets
  Array<uint> sub_dof_maps_offs_;

  //
  mutable uint* local_to_global_;
  uint local_to_global_size_;

  //
  mutable uint* local_to_cell_;
  uint local_to_cell_size_;

  // Provide easy access to map for testing
  std::map<uint, uint> map;
};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline std::string const DofMap::dofmap_signature(
    std::string const& fe_signature)
{
  return SIGN_PREFIX + fe_signature;
}

//-----------------------------------------------------------------------------
inline std::string const DofMap::finite_element_signature(
    std::string const& dofmap_signature)
{
  std::string s(dofmap_signature);
  s.erase(0, SIGN_PREFIX.size());
  return s;
}

//-----------------------------------------------------------------------------
inline std::string const DofMap::make_hash(std::string const& dofmap_signature,
                                           Mesh& mesh)
{
  std::stringstream ss;
  ss << dofmap_signature << "+" << mesh.hash();
  return ss.str();
}

//-----------------------------------------------------------------------------
inline std::string const DofMap::make_hash(ufc::dof_map& ufc_dofmap, Mesh& mesh)
{
  return make_hash(ufc_dofmap.signature(), mesh);
}

//-----------------------------------------------------------------------------
inline char const * DofMap::signature() const
{
  return ufc_dof_map_->signature();
}

//-----------------------------------------------------------------------------
inline bool DofMap::needs_mesh_entities(uint d) const
{
  return ufc_dof_map_->needs_mesh_entities(d);
}

//-----------------------------------------------------------------------------
inline bool DofMap::init_mesh(const ufc::mesh& mesh)
{
  return ufc_dof_map_->init_mesh(mesh);
}

//-----------------------------------------------------------------------------
inline void DofMap::init_cell(const ufc::mesh& m, const ufc::cell& c)
{
  ufc_dof_map_->init_cell(m,c);
}

//-----------------------------------------------------------------------------
inline void DofMap::init_cell_finalize()
{
  ufc_dof_map_->init_cell_finalize();
}

//-----------------------------------------------------------------------------
inline uint DofMap::global_dimension() const
{
  return ufc_dof_map_->global_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::local_dimension() const
{
  return ufc_dof_map_->local_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::geometric_dimension() const
{
  return ufc_dof_map_->geometric_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::num_facet_dofs() const
{
  return ufc_dof_map_->num_facet_dofs();
}

//-----------------------------------------------------------------------------
inline uint DofMap::num_entity_dofs(uint d) const
{
  return ufc_dof_map_->num_entity_dofs(d);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_dofs(uint* dofs, const ufc::mesh& m,
                                  const ufc::cell& c) const
{
  ufc_dof_map_->tabulate_dofs(dofs, m, c);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_dofs(uint* dofs, const ufc::cell& cell) const
{
  ufc_dof_map_->tabulate_dofs(dofs, ufc_mesh, cell);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_facet_dofs(uint* dofs, uint local_facet) const
{
  ufc_dof_map_->tabulate_facet_dofs(dofs, local_facet);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_entity_dofs(uint* dofs, uint d, uint i) const
{
  ufc_dof_map_->tabulate_entity_dofs(dofs, d, i);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_coordinates(real** coordinates,
                                         const ufc::cell& ufc_cell) const
{
  ufc_dof_map_->tabulate_coordinates(coordinates, ufc_cell);
}

//-----------------------------------------------------------------------------
inline uint DofMap::num_sub_dof_maps() const
{
  return ufc_dof_map_->num_sub_dof_maps();
}

//-----------------------------------------------------------------------------
inline ufc::dof_map* DofMap::create_sub_dof_map(uint i) const
{
  return ufc_dof_map_->create_sub_dof_map(i);
}

//-----------------------------------------------------------------------------
inline uint DofMap::macro_local_dimension() const
{
  return ufc_dof_map_->local_dimension();
}

//-----------------------------------------------------------------------------
inline Array<uint> const& DofMap::sub_dof_maps_dimensions() const
{
  return sub_dof_maps_dims_;
}

//-----------------------------------------------------------------------------
inline Array<uint> const& DofMap::sub_dof_maps_offsets() const
{
  return sub_dof_maps_offs_;
}

//-----------------------------------------------------------------------------
inline Mesh& DofMap::mesh() const
{
  return dolfin_mesh;
}

//-----------------------------------------------------------------------------
inline bool DofMap::renumbered() const
{
  return (_dof_map_ || _type_ > -1 || _v_map_);
}

//-----------------------------------------------------------------------------
inline uint DofMap::local_size() const
{
  return _local_size_;
}

//-----------------------------------------------------------------------------
inline uint const * DofMap::dofsmapping() const
{
  if (local_to_global_ == NULL)
    pretabulate_all_dofs();
  return local_to_global_;
}

//-----------------------------------------------------------------------------
inline uint DofMap::dofsmapping_size() const
{
  return local_to_global_size_;
}

//-----------------------------------------------------------------------------
inline uint const * DofMap::cellmapping() const
{
  if (local_to_cell_ == NULL)
    pretabulate_cell_dofs();
  return local_to_cell_;
}

//-----------------------------------------------------------------------------
inline uint DofMap::cellmapping_size() const
{
  return local_to_cell_size_;
}

}

#endif
