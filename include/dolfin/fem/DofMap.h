// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

// Modified by Martin Alnes, 2008

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

class DofMap
{

  static std::string const SIGN_PREFIX;

public:

  static std::string const dofmap_signature(std::string const& fe_signature)
  {
    return SIGN_PREFIX + fe_signature;
  }

  static std::string const finite_element_signature(
      std::string const& dofmap_signature)
  {
    std::string s(dofmap_signature);
    s.erase(0, SIGN_PREFIX.size());
    return s;
  }

  static std::string const make_hash(std::string const& dofmap_signature,
                                     Mesh& mesh)
  {
    std::stringstream ss;
    ss << dofmap_signature << "+" << mesh.hash();
    return ss.str();
  }

  static std::string const make_hash(ufc::dof_map& ufcdofmap, Mesh& mesh)
  {
    return make_hash(ufcdofmap.signature(), mesh);
  }

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

  /// Destructor
  ~DofMap();

  /// Return a string identifying the dof map
  char const * signature() const;

  /// Return the dimension of the global finite element function space
  uint global_dimension() const;

  /// Return the dimension of the local finite element function space
  uint local_dimension() const;

  /// Return the dimension of the local finite element function space
  uint macro_local_dimension() const;

  /// Return number of facet dofs
  uint num_facet_dofs() const;

  /// Return local to global mapping
  uint const * dofsmapping() const;

  /// Return the size of the local to global mapping
  uint dofsmapping_size() const;

  /// Return local cell mapping from component grouped numbering to node grouped
  uint const * cellmapping() const;

  /// Return the size of the local cell mapping
  uint cellmapping_size() const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(uint* dofs, ufc::cell& ufc_cell, uint cell_index);

  /// Tabulate the local-to-global mapping of dofs on a cell
  void
  tabulate_dofs(uint* dofs, const ufc::cell& ufc_cell, uint cell_index) const;

  /// Tabulate local-local facet dofs
  void tabulate_facet_dofs(uint* dofs, uint local_facet) const;

  // FIXME: Can this function eventually be removed?
  /// Tabulate the local-to-global mapping of dofs on a ufc cell
  void tabulate_dofs(uint* dofs, const ufc::cell& cell) const;

  void
  tabulate_coordinates(real** coordinates, const ufc::cell& ufc_cell) const;

  /// Extract sub dof map
  DofMap* extractDofMap(const Array<uint>& sub_system, uint& offset) const;

  /// Return mesh associated with map
  Mesh& mesh() const;

  /// Return renumbering (used for testing)
  std::map<uint, uint> getMap();  // const;

  ///
  std::string const& mesh_hash() const;

  ///
  std::string const& hash() const;

  /// Display mapping
  void disp() const;

  bool renumbered();

  uint local_size();

private:

  /// Initialise DofMap
  void init();

  /// Build parallel dof map
  void build();

  ///
  void pretabulate_all_dofs() const;

  ///
  void pretabulate_cell_dofs() const;

  /// Extract sub DofMap
  ufc::dof_map* extractDofMap(const ufc::dof_map& dof_map, uint& offset,
                              const Array<uint>& sub_system) const;

  // Local UFC dof map
  bool const ufc_dof_map_local;

  // UFC dof map
  ufc::dof_map * const ufc_dof_map;

  // Hash
  std::string const mesh_hash_;
  std::string const hash_;

  // Parallel dof map
  uint* dof_map;

  //
  mutable uint* local_to_global_;
  uint local_to_global_size_;

  //
  mutable uint* local_to_cell_;
  uint local_to_cell_size_;

  // UFC mesh
  UFCMesh ufc_mesh;

  // DOLFIN mesh
  Mesh& dolfin_mesh;

  // Number of cells in the mesh
  uint num_cells;

  // Partitions
  MeshFunction<uint>* partitions;

  // Provide easy access to map for testing
  std::map<uint, uint> map;

  int _type_;
  uint _offset_;
  uint _local_size;

  uint *v_map;
};

//--- Inlined -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline char const * DofMap::signature() const
{
  return ufc_dof_map->signature();
}

//-----------------------------------------------------------------------------
inline uint DofMap::global_dimension() const
{
  return ufc_dof_map->global_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::local_dimension() const
{
  return ufc_dof_map->local_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::macro_local_dimension() const
{
  return ufc_dof_map->local_dimension();
}

//-----------------------------------------------------------------------------
inline uint DofMap::num_facet_dofs() const
{
  return ufc_dof_map->num_facet_dofs();
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

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_facet_dofs(uint* dofs, uint local_facet) const
{
  ufc_dof_map->tabulate_facet_dofs(dofs, local_facet);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_dofs(uint* dofs, const ufc::cell& cell) const
{
  ufc_dof_map->tabulate_dofs(dofs, ufc_mesh, cell);
}

//-----------------------------------------------------------------------------
inline void DofMap::tabulate_coordinates(real** coordinates,
                                         const ufc::cell& ufc_cell) const
{
  ufc_dof_map->tabulate_coordinates(coordinates, ufc_cell);
}

//-----------------------------------------------------------------------------
inline Mesh& DofMap::mesh() const
{
  return dolfin_mesh;
}

//-----------------------------------------------------------------------------
inline bool DofMap::renumbered()
{
  return (dof_map > 0 || _type_ > -1 || v_map > 0);
}

//-----------------------------------------------------------------------------
inline uint DofMap::local_size()
{
  return _local_size;
}

}

#endif
