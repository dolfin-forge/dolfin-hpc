// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman 2007.
// Modified by Magnus Vikstrøm 2007.
// Modified by Garth N. Wells 2007.
// Modified by Balthasar Reuter, 2013.
// Modified by Aurélien Larcher, 2014.
//
// First added:  2006-05-08
// Last changed: 2013-03-22

#ifndef __DOLFIN_MESH_H
#define __DOLFIN_MESH_H

#include <dolfin/common/Variable.h>

#include <dolfin/common/types.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshTopology.h>

#include <string>

namespace dolfin
{

class BoundaryMesh;
class Cell;
class IntersectionDetector;
class MappedManifold;
template<class T, class E, uint N = 1> class MeshValues;
class PeriodicSubDomain;
class Space;
class Vertex;

/// A Mesh consists of a set of connected and numbered mesh entities.
///
/// Both the representation and the interface are dimension-independent,
/// but a concrete interface is also provided for standard named mesh
/// entities:
///
///     Entity  Dimension  Codimension
///
///     Vertex      0           -
///     Edge        1           -
///     Face        2           -
///
///     Facet       -           1
///     Cell        -           0
///
/// When working with mesh iterators, all entities and connectivity
/// are precomputed automatically the first time an iterator is
/// created over any given topological dimension or connectivity.
///
/// Note that for efficiency, only entities of dimension zero
/// (vertices) and entities of the maximal dimension (cells) exist
/// when creating a Mesh. Other entities must be explicitly created
/// by calling init(). For example, all edges in a mesh may be created
/// by a call to mesh.init(1). Similarly, connectivities such as
/// all edges connected to a given vertex must also be explicitly
/// created (in this case by a call to mesh.init(0, 1)).

class Mesh : public Variable
{

  friend class MeshEditor;

public:

  /// Create empty mesh
  Mesh();

  /// Constructor from cell type and space
  Mesh(CellType const& type, Space const& space);

  /// Copy constructor
  Mesh(Mesh const& mesh);

  /// Create mesh from data file
  Mesh(std::string const& filename);

  /// Destructor
  virtual ~Mesh();

  /// Assignment
  Mesh const& operator=(Mesh const& mesh);

  /// Swap instances
  void swap(Mesh& other);

  /// Identity
  bool operator ==(Mesh const& other) const;
  bool operator !=(Mesh const& other) const;

  /// Return if mesh is empty
  bool empty() const;

  /// Return mesh cell type
  CellType const& type() const;

  /// Return mesh space
  Space const& space() const;

  //--- TOPOLOGY --------------------------------------------------------------

  /// Return mesh topology (non-const version)
  MeshTopology& topology();

  /// Return mesh topology (const)
  MeshTopology const& topology() const;

  /// Return number of entities of given topological dimension
  uint size(uint dim) const;

  /// Return number of vertices of mesh partition
  uint num_vertices() const;

  /// Return number of cells of mesh partition
  uint num_cells() const;

  /// Return connectivity for all cells
  uint * cells();

  /// Return connectivity for all cells
  uint const * cells() const;

  /// Compute entities of given topological dimension
  void init(uint dim) const;

  /// Compute connectivity between given pair of dimensions
  void init(uint d0, uint d1) const;

  /// Compute all entities and connectivity
  void init() const;

  /// Return if the topology should be reordered
  virtual bool reordering() const { return true; }

  //---

  /// Return exterior boundary of the mesh
  BoundaryMesh& exterior_boundary();

  /// Return interior boundary of the mesh
  BoundaryMesh& interior_boundary();

  //---

  //// Return if the mesh should be read/written in serial
  bool serial_io() const;

  //// Return if the mesh should be read/written in parallel
  bool parallel_io() const;

  /// Return whether the mesh is distributed i.e iff the topology is distributed
  bool is_distributed() const;

  /// Return mesh distribution data (non-const version)
  MeshDistributedData& distdata();

  /// Return mesh distribution data (const)
  MeshDistributedData const& distdata() const;

  /// Return global number of entities of given topological dimension
  uint global_size(uint dim) const;

  /// Return number of vertices of global mesh
  uint num_global_vertices() const;

  /// Return number of cells of global mesh
  uint num_global_cells() const;

  //--- GEOMETRY --------------------------------------------------------------

  /// Return mesh geometry (non-const version)
  MeshGeometry& geometry();

  /// Return mesh geometry (const)
  MeshGeometry const& geometry() const;

  //---

  /// Return intersection detector for the mesh
  IntersectionDetector& intersector();

  //--- PERIODICITY -----------------------------------------------------------

  /// Return whether the mesh has periodic subdomain
  bool has_periodic_constraint() const;

  /// Add a periodic subdomain
  void add_periodic_constraint(PeriodicSubDomain const& periodic);

  /// Return the list of periodic mappings
  Array<MappedManifold *> const& periodic_mappings() const;

  //---------------------------------------------------------------------------

  /*
   *  Mesh partitioning routines
   */

  /// Partition mesh into num_processes partitions
  void partition(MeshValues<uint, Cell>& partitions);

  /// Partition mesh into num_partitions = numProc with weights on vertices
  /// ^H^H^H on the *fucking* *cells*
  void partition(MeshValues<uint, Cell>& partitions, MeshValues<uint, Cell>& weight);

  /// Partition mesh into num_partitions = numProc
  void partition_geom(MeshValues<uint, Vertex>& partitions);

  /*
   *  Mesh distribution routines
   */

  /// Distribute a mesh according to a mesh function
  void distribute();

  /// Distribute a mesh according to a mesh function
  void distribute(MeshValues<uint, Cell>& distribution);
  void distribute(MeshValues<uint, Vertex>& distribution);

  /*
   *  Mesh refinement routines
   */

  /// Refine mesh uniformly
  void refine();

  /// Refine mesh according to cells marked for refinement
  void refine(MeshValues<bool, Cell>& cell_markers, bool refine_boundary = true,
              bool load_balance = true);

  //---

  /// Return hash to identify the state of the mesh
  std::string const hash() const;

  /// Display mesh data
  void disp() const;

private:

  /// Intialize mesh given cell type and space
  void init(CellType const& type, Space const& space);

  /// Clear all mesh data
  void clear();

  // Cell type
  CellType * cell_type_;

  // Mesh topology
  MeshTopology topology_;

  // Space
  Space * space_;

  // Mesh geometry
  MeshGeometry geometry_;

  /// Exterior boundary mesh
  mutable BoundaryMesh * exterior_boundary_;

  /// Interior boundary mesh
  mutable BoundaryMesh * interior_boundary_;

  /// Intersection detector
  mutable IntersectionDetector * intersection_detector_;

  /// Periodic constraints
  mutable Array<MappedManifold *> periodic_mappings_;

  int timestamp_;

};

//--- TEMPLATE SPECIALIZATIONS ------------------------------------------------

template<class E>
inline uint entity_dimension(Mesh& m) { return dimension<E>(m.type()); }

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_H */
