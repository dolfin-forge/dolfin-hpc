// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_H
#define __DOLFIN_MESH_H

#include <dolfin/common/Variable.h>

#include <dolfin/common/types.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/parameter/parameters.h>

#include <string>

namespace dolfin
{

class BoundaryMesh;
class Cell;
class IntersectionDetector;
class MappedManifold;
class MeshData;
template<class T, class E, uint N = 1> struct MeshValues;
class PeriodicSubDomain;
struct Space;
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

  /// Constructor from cell type, space, and communicator
  Mesh(CellType const& ctype, Space const& space, Comm& comm = DOLFIN_COMM );

  /// Copy constructor
  Mesh(Mesh const& mesh);

  /// Create mesh from data file
  Mesh(std::string const& filename);

  /// Destructor
  virtual ~Mesh();

  /// Swap instances
  friend void swap( Mesh& a, Mesh& b );

  /// Assignment
  auto operator=(Mesh const& other) -> Mesh const&;

  /// Identity
  auto operator ==(Mesh const& other) const -> bool;
  auto operator !=(Mesh const& other) const -> bool;

  /// Return mesh cell type
  auto type() const -> CellType const&;

  /// Return mesh space
  auto space() const -> Space const&;

  /// Return a Point describing the minimum extent of the mesh
  /// this point might not be part of the mesh itself
  auto extent_min() const -> Point const &;

  /// Return a Point describing the maximum extent of the mesh
  /// this point might not be part of the mesh itself
  auto extent_max() const -> Point const &;

  /// update the mesh extent
  void extent_update();

  //--- TOPOLOGY --------------------------------------------------------------

  /// Return mesh topology (non-const version)
  auto topology() -> MeshTopology&;

  /// Return mesh topology (const)
  auto topology() const -> MeshTopology const&;

  /// Return topological dimension
  auto topology_dimension() const -> uint;

  /// Return number of entities of given topological dimension
  auto size(uint dim) const -> uint;

  /// Return number of vertices of mesh partition
  auto num_vertices() const -> uint;

  /// Return number of cells of mesh partition
  auto num_cells() const -> uint;

  /// Return connectivity for all cells
  auto cells() -> Array< Array< uint > > &; //!< @tod remove this function

  /// Return connectivity for all cells
  auto cells() const -> Array< Array< uint > > const &; //!< @tod remove this function

  /// Compute entities of given topological dimension
  void init(uint dim) const;

  /// Compute connectivity between given pair of dimensions
  void init(uint d0, uint d1) const;

  /// Compute all entities and connectivity
  void init() const;

  /// Return if the topology should be reordered
  virtual auto reordering() const -> bool { return true; }

  //---

  /// Return exterior boundary of the mesh
  auto exterior_boundary() -> BoundaryMesh&;

  /// Return interior boundary of the mesh
  auto interior_boundary() -> BoundaryMesh&;

  //---

  //// Return if the mesh should be read/written in serial
  auto serial_io() const -> bool;

  //// Return if the mesh should be read/written in parallel
  auto parallel_io() const -> bool;

  /// Return whether the mesh is distributed i.e iff the topology is distributed
  auto is_distributed() const -> bool; //!< @todo remove this function

  /// Return mesh distribution data (non-const version)
  auto distdata() -> MeshDistributedData&; //!< @todo remove this function

  /// Return mesh distribution data (const)
  auto distdata() const -> MeshDistributedData const&; //!< @tod remove this function

  /// Return global number of entities of given topological dimension
  auto global_size(uint dim) const -> uint;

  /// Return number of vertices of global mesh
  auto num_global_vertices() const -> uint;

  /// Return number of cells of global mesh
  auto num_global_cells() const -> uint;

  //--- GEOMETRY --------------------------------------------------------------

  /// Return mesh geometry (non-const version)
  auto geometry() -> MeshGeometry&;

  /// Return mesh geometry (const)
  auto geometry() const -> MeshGeometry const&;

  /// Return geometric dimension
  auto geometry_dimension() const -> uint;

  //---

  /// Return intersection detector for the mesh
  auto intersector() -> IntersectionDetector&;

  //--- PERIODICITY -----------------------------------------------------------

  /// Return whether the mesh has periodic subdomain
  auto has_periodic_constraint() const -> bool;

  /// Add a periodic subdomain
  void add_periodic_constraint(PeriodicSubDomain const& periodic);

  /// Return the list of periodic mappings
  auto periodic_mappings() const -> Array<MappedManifold *> const&;

  //---------------------------------------------------------------------------

  /*
   *  Mesh partitioning routines
   */

  /// Partition mesh into num_processes partitions
  void partition(MeshValues<uint, Cell>& partitions);

  /// Partition mesh into num_partitions = numProc with weights
  /// on the vertices of the dual graph
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

  /// Distribute a mesh according to a cell function and transfer data
  void distribute(MeshValues<uint, Cell>& distribution, MeshData& data);

  /*
   *  Mesh refinement routines
   */

  /// Refine mesh uniformly
  void refine();

  //---------------------------------------------------------------------------

  /// Return hash to identify the state of the mesh
  auto hash() const -> std::string const;

  /// Display mesh data
  void disp() const;

  /// Check mesh consistency
  void check() const;

private:

  // Mesh topology
  MeshTopology topology_;

  // Mesh geometry
  MeshGeometry geometry_;

  /// Exterior boundary mesh
  mutable BoundaryMesh * exterior_boundary_{nullptr};

  /// Interior boundary mesh
  mutable BoundaryMesh * interior_boundary_{nullptr};

  /// Intersection detector
  mutable IntersectionDetector * intersection_detector_{nullptr};

  /// Periodic constraints
  mutable Array<MappedManifold *> periodic_mappings_;

  Point extent_min_;
  Point extent_max_;

  int timestamp_;
};

//-----------------------------------------------------------------------------
inline auto Mesh::operator!=( Mesh const & other ) const -> bool
{
  return !( *this == other );
}

//-----------------------------------------------------------------------------
inline auto Mesh::type() const -> CellType const &
{
  return topology_.type();
}

//-----------------------------------------------------------------------------
inline auto Mesh::topology() -> MeshTopology &
{
  return topology_;
}

//-----------------------------------------------------------------------------
inline auto Mesh::topology() const -> MeshTopology const &
{
  return topology_;
}

//-----------------------------------------------------------------------------
inline auto Mesh::topology_dimension() const -> uint
{
  return topology_.dim();
}

//-----------------------------------------------------------------------------
inline auto Mesh::size( uint dim ) const -> uint
{
  return topology_.size( dim );
}

//-----------------------------------------------------------------------------
inline auto Mesh::num_vertices() const -> uint
{
  return topology_.size( 0 );
}

//-----------------------------------------------------------------------------
inline auto Mesh::num_cells() const -> uint
{
  return topology_.size( topology_.dim() );
}

//-----------------------------------------------------------------------------
inline auto Mesh::cells() -> Array< Array< uint > > &
{
  return topology_( topology_.dim(), 0 )();
}

//-----------------------------------------------------------------------------
inline auto Mesh::cells() const -> Array< Array< uint > > const &
{
  return topology_( topology_.dim(), 0 )();
}

//-----------------------------------------------------------------------------
inline void Mesh::init( uint dim ) const
{
  topology_.size( dim );
}

//-----------------------------------------------------------------------------
inline void Mesh::init( uint d0, uint d1 ) const
{
  topology_( d0, d1 ).order();
}

//-----------------------------------------------------------------------------
inline auto Mesh::serial_io() const -> bool
{
  return ( PE::size() == 1 ) || dolfin_get< bool >( "Mesh read in serial" );
}

//-----------------------------------------------------------------------------
inline auto Mesh::parallel_io() const -> bool
{
  return !this->serial_io();
}

//-----------------------------------------------------------------------------
inline auto Mesh::is_distributed() const -> bool
{
  return topology().distributed();
}

//-----------------------------------------------------------------------------
inline auto Mesh::distdata() -> MeshDistributedData &
{
  return topology().distdata();
}

//-----------------------------------------------------------------------------
inline auto Mesh::distdata() const -> MeshDistributedData const &
{
  return topology().distdata();
}

//-----------------------------------------------------------------------------
inline auto Mesh::global_size( uint dim ) const -> uint
{
  return topology_.global_size( dim );
}

//-----------------------------------------------------------------------------
inline auto Mesh::num_global_vertices() const -> uint
{
  return topology_.global_size( 0 );
}

//-----------------------------------------------------------------------------
inline auto Mesh::num_global_cells() const -> uint
{
  return topology_.global_size( topology_.dim() );
}

//-----------------------------------------------------------------------------
inline auto Mesh::space() const -> Space const &
{
  return geometry_.space();
}

//-----------------------------------------------------------------------------
inline auto Mesh::geometry() -> MeshGeometry &
{
  return geometry_;
}

//-----------------------------------------------------------------------------
inline auto Mesh::geometry() const -> MeshGeometry const &
{
  return geometry_;
}

//-----------------------------------------------------------------------------
inline auto Mesh::geometry_dimension() const -> uint
{
  return geometry_.dim();
}

//-----------------------------------------------------------------------------
inline auto Mesh::extent_min() const -> Point const &
{
  return extent_min_;
}

//-----------------------------------------------------------------------------
inline auto Mesh::extent_max() const -> Point const &
{
  return extent_max_;
}

//--- TEMPLATE SPECIALIZATIONS ------------------------------------------------

template<class E>
inline auto entity_dimension(Mesh& m) -> uint { return dimension<E>(m.type()); }

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_H */
