// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

// Modified by Martin Alnes, 2008

// First added:  2007-03-01
// Last changed: 2008-04-10

#ifndef __DOF_MAP_H
#define __DOF_MAP_H

#include <map>
#include <ufc.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include "UFCCell.h"
#include "UFCMesh.h"

namespace dolfin
{
class SubSytem;
class UFC;

/// This class handles the mapping of degrees of freedom.
/// It wraps a ufc::dof_map on a specific mesh and provides
/// optional precomputation and reordering of dofs.

class DofMap
{
public:

	static std::string const dofmap_signature(std::string const& fe_signature)
	{
	#if UFC_VERSION_MAJOR == 1
		return "FFC dof map for "+fe_signature;
	#elif UFC_VERSION_MAJOR == 2
		return "FFC dofmap for " + fe_signature;
	#else
		error("Invalid UFC version, supported major are 1.1.x and 2.x");
		return "";
	#endif
	}

	static std::string const make_hash(std::string const& dofmap_signature, Mesh& mesh);

	static std::string const make_hash(ufc::dof_map& ufcdofmap, Mesh& mesh);

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
	unsigned int global_dimension() const;

	/// Return the dimension of the local finite element function space
	unsigned int local_dimension() const;

	/// Return the dimension of the local finite element function space
	unsigned int macro_local_dimension() const;

	/// Return number of facet dofs
	unsigned int num_facet_dofs() const;

	/// Tabulate the local-to-global mapping of dofs on a cell
	void tabulate_dofs(uint* dofs, ufc::cell& ufc_cell, uint cell_index);

	/// Tabulate the local-to-global mapping of dofs on a cell
	void tabulate_dofs(uint* dofs, const ufc::cell& ufc_cell,
						uint cell_index) const;

	/// Tabulate local-local facet dofs
	void tabulate_facet_dofs(uint* dofs, uint local_facet) const;

	// FIXME: Can this function eventually be removed?
	/// Tabulate the local-to-global mapping of dofs on a ufc cell
	void tabulate_dofs(uint* dofs, const ufc::cell& cell) const;

	void tabulate_coordinates(real** coordinates,
								const ufc::cell& ufc_cell) const;

	/// Extract sub dof map
	DofMap* extractDofMap(const Array<uint>& sub_system, uint& offset) const;

	/// Return mesh associated with map
	Mesh& mesh() const;

	/// Build parallel dof map
	void build();

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
inline unsigned int DofMap::global_dimension() const
{
	return ufc_dof_map->global_dimension();
}

//-----------------------------------------------------------------------------
inline unsigned int DofMap::local_dimension() const
{
	return ufc_dof_map->local_dimension();
}

//-----------------------------------------------------------------------------
inline unsigned int DofMap::macro_local_dimension() const
{
	return ufc_dof_map->local_dimension();
}

//-----------------------------------------------------------------------------
inline unsigned int DofMap::num_facet_dofs() const
{
	return ufc_dof_map->num_facet_dofs();
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

//-----------------------------------------------------------------------------
inline std::string const DofMap::make_hash(std::string const& ufcdofmap_sign,
											Mesh& mesh)
{
	std::stringstream ss;
	ss << ufcdofmap_sign << "+" << mesh.hash();
	return ss.str();
}

//-----------------------------------------------------------------------------
inline std::string const DofMap::make_hash(ufc::dof_map& ufcdofmap,
											Mesh& mesh)
{
	return make_hash(ufcdofmap.signature(), mesh);
}

}

#endif
