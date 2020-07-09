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

  /// Constructor from cell type and space
  Mesh(CellType const& ctype, Space const& space);

  /// Constructor from cell type, space, and communicator
  Mesh(CellType const& ctype, Space const& space, Comm& comm);

  /// Copy constructor
  Mesh(Mesh const& mesh);

  /// Create mesh from data file
  Mesh(std::string const& filename);

  /// Destructor
  virtual ~Mesh();

  /// Swap instances
  friend void swap( Mesh& a, Mesh& b );

  /// Assignment
  Mesh const& operator=(Mesh const& other);

  /// Identity
  bool operator ==(Mesh const& other) const;
  bool operator !=(Mesh const& other) const;

  /// Return mesh cell type
  CellType const& type() const;

  /// Return mesh space
  Space const& space() const;

  //--- TOPOLOGY --------------------------------------------------------------

  /// Return mesh topology (non-const version)
  MeshTopology& topology();

  /// Return mesh topology (const)
  MeshTopology const& topology() const;

  /// Return topological dimension
  uint topology_dimension() const;

  /// Return number of entities of given topological dimension
  uint size(uint dim) const;

  /// Return number of vertices of mesh partition
  uint num_vertices() const;

  /// Return number of cells of mesh partition
  uint num_cells() const;

  /// Return connectivity for all cells
  Array< Array< uint > > & cells(); //!< @tod remove this function

  /// Return connectivity for all cells
  Array< Array< uint > > const & cells() const; //!< @tod remove this function

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
  bool is_distributed() const; //!< @todo remove this function

  /// Return mesh distribution data (non-const version)
  MeshDistributedData& distdata(); //!< @todo remove this function

  /// Return mesh distribution data (const)
  MeshDistributedData const& distdata() const; //!< @tod remove this function

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

  /// Return geometric dimension
  uint geometry_dimension() const;

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
  std::string const hash() const;

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
  mutable BoundaryMesh * exterior_boundary_;

  /// Interior boundary mesh
  mutable BoundaryMesh * interior_boundary_;

  /// Intersection detector
  mutable IntersectionDetector * intersection_detector_;

  /// Periodic constraints
  mutable Array<MappedManifold *> periodic_mappings_;

  int timestamp_;
};

//-----------------------------------------------------------------------------
inline bool Mesh::operator!=( Mesh const & other ) const
{
  return !( *this == other );
}

//-----------------------------------------------------------------------------
inline CellType const & Mesh::type() const
{
  return topology_.type();
}

//-----------------------------------------------------------------------------
inline MeshTopology & Mesh::topology()
{
  return topology_;
}

//-----------------------------------------------------------------------------
inline MeshTopology const & Mesh::topology() const
{
  return topology_;
}

//-----------------------------------------------------------------------------
inline uint Mesh::topology_dimension() const
{
  return topology_.dim();
}

//-----------------------------------------------------------------------------
inline uint Mesh::size( uint dim ) const
{
  return topology_.size( dim );
}

//-----------------------------------------------------------------------------
inline uint Mesh::num_vertices() const
{
  return topology_.size( 0 );
}

//-----------------------------------------------------------------------------
inline uint Mesh::num_cells() const
{
  return topology_.size( topology_.dim() );
}

//-----------------------------------------------------------------------------
inline Array< Array< uint > > & Mesh::cells()
{
  return topology_( topology_.dim(), 0 )();
}

//-----------------------------------------------------------------------------
inline Array< Array< uint > > const & Mesh::cells() const
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
inline bool Mesh::serial_io() const
{
  return ( PE::size() == 1 ) || dolfin_get< bool >( "Mesh read in serial" );
}

//-----------------------------------------------------------------------------
inline bool Mesh::parallel_io() const
{
  return !this->serial_io();
}

//-----------------------------------------------------------------------------
inline bool Mesh::is_distributed() const
{
  return topology().distributed();
}

//-----------------------------------------------------------------------------
inline MeshDistributedData & Mesh::distdata()
{
  return topology().distdata();
}

//-----------------------------------------------------------------------------
inline MeshDistributedData const & Mesh::distdata() const
{
  return topology().distdata();
}

//-----------------------------------------------------------------------------
inline uint Mesh::global_size( uint dim ) const
{
  return topology_.global_size( dim );
}

//-----------------------------------------------------------------------------
inline uint Mesh::num_global_vertices() const
{
  return topology_.global_size( 0 );
}

//-----------------------------------------------------------------------------
inline uint Mesh::num_global_cells() const
{
  return topology_.global_size( topology_.dim() );
}

//-----------------------------------------------------------------------------
inline Space const & Mesh::space() const
{
  return geometry_.space();
}

//-----------------------------------------------------------------------------
inline MeshGeometry & Mesh::geometry()
{
  return geometry_;
}

//-----------------------------------------------------------------------------
inline MeshGeometry const & Mesh::geometry() const
{
  return geometry_;
}

//-----------------------------------------------------------------------------
inline uint Mesh::geometry_dimension() const
{
  return geometry_.dim();
}

//--- TEMPLATE SPECIALIZATIONS ------------------------------------------------

template<class E>
inline uint entity_dimension(Mesh& m) { return dimension<E>(m.type()); }

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_H */
