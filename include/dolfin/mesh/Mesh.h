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

#ifndef __MESH_H
#define __MESH_H

#include <string>

#include <dolfin/common/types.h>
#include <dolfin/common/Variable.h>
#include "ALEType.h"
#include "MeshDistributedData.h"
#include "MeshGeometry.h"
#include "MeshTopology.h"
#include "CellType.h"

#ifdef HAVE_LIBGEOM
namespace libgeom
{
  class Geometry;
}
#endif

namespace dolfin
{

class BoundaryMesh;
class IntersectionDetector;
template<class T> class MeshFunction;
class MeshData;

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

class Mesh: public Variable
{

  // Friends
  friend class MeshEditor;
  friend class MeshTopology;
  friend class MeshGeometry;
  friend class MPIMeshCommunicator;

public:

  /// Create empty mesh
  Mesh();

  /// Copy constructor
  Mesh(Mesh const& mesh);

  /// Create mesh from data file
  Mesh(std::string filename);

  /// Destructor
  ~Mesh();

  /// Assignment
  const Mesh& operator=(const Mesh& mesh);

  /// Return number of vertices of mesh partition
  uint numVertices() const;

  /// Return number of ghost vertices
  uint numGhostVertices() const;

  /// Return number of vertices of global mesh
  uint global_numVertices() const;

  /// Return number of edges of mesh partition
  uint numEdges() const;

  /// Return number of edges of global mesh
  uint global_numEdges() const;

  /// Return number of faces of mesh partition
  uint numFaces() const;

  /// Return number of faces of global mesh
  uint global_numFaces() const;

  /// Return number of facets of mesh partition
  uint numFacets() const;

  /// Return number of facets of global mesh
  uint global_numFacets() const;

  /// Return number of cells of mesh partition
  uint numCells() const;

  /// Return number of cells of global mesh
  uint global_numCells() const;

  /// Return coordinates of all vertices
  real* coordinates();

  /// Return coordinates of all vertices
  const real* coordinates() const;

  /// Return connectivity for all cells
  uint* cells();

  /// Return connectivity for all cells
  const uint* cells() const;

  /// Return number of entities of given topological dimension
  uint size(uint dim) const;

  /// Return mesh topology (non-const version)
  MeshTopology& topology();

  /// Return mesh topology (const version)
  const MeshTopology& topology() const;

  /// Return mesh geometry (non-const version)
  MeshGeometry& geometry();

  /// Return mesh geometry (const version)
  const MeshGeometry& geometry() const;

  /// Return mesh distribution data
  MeshDistributedData& distdata();

  /// Return mesh distribution data (const version)
  const MeshDistributedData& distdata() const;

  /// Return mesh data
  MeshData& data();

  /// Return mesh cell type
  CellType& type();

  /// Return mesh cell type
  const CellType& type() const;

  /// Return exterior boundary of the mesh
  BoundaryMesh& exterior_boundary();

  /// Return interior boundary of the mesh
  BoundaryMesh& interior_boundary();

  /// Return intersection detector for the mesh
  IntersectionDetector& intersector();

  /// Compute entities of given topological dimension and return number of entities
  uint init(uint dim) const;

  /// Compute connectivity between given pair of dimensions
  void init(uint d0, uint d1) const;

  /// Compute all entities and connectivity
  void init() const;

  /// Clear all mesh data
  void clear();

  /// Order all mesh entities (not needed if "mesh order entities" is set)
  void order();

  /// Return true iff topology is ordered according to the UFC numbering
  bool ordered() const;

  /// Refine mesh uniformly
  void refine();

#ifdef HAVE_LIBGEOM
  /// Refine mesh uniformly including geometry informations -surfaces
  void refine(libgeom::Geometry& geom, MeshFunction<int>& patch_id_list,
              MeshFunction<float>& bnd_u, MeshFunction<float>& bnd_v);

  /// Refine mesh uniformly including geometry informations -curves
  void refine(libgeom::Geometry& geom, MeshFunction<int>& patch_id_list,
              MeshFunction<float>& bnd_u);
#endif

  /// Refine mesh according to cells marked for refinement
  void refine(MeshFunction<bool>& cell_markers, bool refine_boundary = true,
              bool load_balance = true);

  /// Coarsen mesh uniformly
  void coarsen();

  /// Coarsen mesh according to cells marked for coarsening
  void coarsen(MeshFunction<bool>& cell_markers, bool coarsen_boundary = false);

  /// Move coordinates of mesh according to new boundary coordinates
  void move(Mesh& boundary, ALEType method = lagrange);

  /// Smooth mesh using Lagrangian mesh smoothing
  void smooth();

  /// Partition mesh into num_processes partitions
  void partition(MeshFunction<uint>& partitions);

  /// Partition mesh into num_partitions partitions
  void partition(MeshFunction<uint>& partitions, uint num_partitions);

  /// Partition mesh into num_partitions = numProc with weights on vertices
  void partition(MeshFunction<uint>& partitions, MeshFunction<uint>& weight);

  /// Partition mesh into num_partitions = numProc
  void partition_geom(MeshFunction<uint>& partitions);

  /// Distribute a mesh according to a mesh function
  void distribute(MeshFunction<uint>& distribution);

  /// Distribute a mesh according to a mesh function and transfer marked cells
  void distribute(MeshFunction<uint>& distribution,
                  MeshFunction<bool>& cell_markers,
                  MeshFunction<bool>& new_cell_markers);

  /// Distribute a mesh according to a mesh function and transfer cell functions
  ///
  /// cell_functions contains pairs as <old_function,new_function>
  void distribute(MeshFunction<uint>& distribution,
                  Array<std::pair<MeshFunction<uint> *,
                  MeshFunction<uint> *> >& cell_functions);

  /// Distribute a mesh according to a mesh function and transfer vertex functions
  ///
  /// vertex_functions contains pairs as <old_function,new_function>
  void distribute(MeshFunction<uint>& distribution,
                  Array<std::pair<MeshFunction<double> *,
                  MeshFunction<double> *> >& vertex_functions);

  /// Distribute a mesh according to a mesh function and transfer cell and
  /// vertex functions
  ///
  /// cell_functions contains pairs as <old_function,new_function>
  ///
  /// vertex_functions contains pairs as <old_function,new_function>
  void distribute(MeshFunction<uint>& distribution,
                  Array<std::pair<MeshFunction<uint> *,
                  MeshFunction<uint> *> >& cell_functions,
                  Array<std::pair<MeshFunction<double> *,
                  MeshFunction<double> *> >& vertex_functions);

  /// Return whether the mesh is distributed
  bool distributed() const;

  /// Renumber mesh global numbering
  void renumber();

  /// Return hash to identify the state of the mesh
  std::string const hash() const;

  /// Display mesh data
  void disp() const;

  /// Return a short desriptive string
  std::string str() const;

  /// Output
  friend LogStream& operator<<(LogStream& stream, const Mesh& mesh);

private:

  bool _is_distributed;

  // Mesh topology
  MeshTopology _topology;

  // Mesh geometry
  MeshGeometry _geometry;

  // Auxiliary mesh data
  MeshData* _data;

  // Cell type
  CellType* _cell_type;

  /// Exterior boundary mesh
  mutable BoundaryMesh * _exterior_boundary;

  /// Interior boundary mesh
  mutable BoundaryMesh * _interior_boundary;

  /// Intersection detector
  mutable IntersectionDetector * _intersection_detector;

  int _timestamp;

};
//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline uint Mesh::numVertices() const
{
  return _topology.size(0);
}

//-----------------------------------------------------------------------------
inline uint Mesh::numGhostVertices() const
{
  return ((dolfin::MPI::numProcesses() == 1) ? 0 : distdata().num_ghost(0));
}

//-----------------------------------------------------------------------------
inline uint Mesh::global_numVertices() const
{
  return (
      (dolfin::MPI::numProcesses() == 1) ?
          numVertices() : distdata().global_numVertices());
}

//-----------------------------------------------------------------------------
inline uint Mesh::numEdges() const
{
  return _topology.size(1);
}

//-----------------------------------------------------------------------------
inline uint Mesh::global_numEdges() const
{
  return (
      (dolfin::MPI::numProcesses() == 1) ?
          numEdges() : distdata().global_numEdges());
}

//-----------------------------------------------------------------------------
inline uint Mesh::numFaces() const
{
  return _topology.size(2);
}

//-----------------------------------------------------------------------------
inline uint Mesh::global_numFaces() const
{
  return (
      (dolfin::MPI::numProcesses() == 1) ?
          numFaces() : distdata().global_numFaces());
}

//-----------------------------------------------------------------------------
inline uint Mesh::numFacets() const
{
  return _topology.size(_topology.dim() - 1);
}

//-----------------------------------------------------------------------------
inline uint Mesh::global_numFacets() const
{
  uint ret = numFacets();

  if (dolfin::MPI::numProcesses() > 1)
  {
    switch(_topology.dim())
    {
      case 3:
        ret= distdata().global_numFaces();
        break;
      case 2:
        ret= distdata().global_numEdges();
        break;
      case 1:
        ret= distdata().global_numVertices();
        break;
    }
  }

  return ret;
}

//-----------------------------------------------------------------------------
inline uint Mesh::numCells() const
{
  return _topology.size(_topology.dim());
}

//-----------------------------------------------------------------------------
inline uint Mesh::global_numCells() const
{
  return (
      (dolfin::MPI::numProcesses() == 1) ?
          numCells() : distdata().global_numCells());
}

//-----------------------------------------------------------------------------
inline real* Mesh::coordinates()
{
  return _geometry.x();
}

//-----------------------------------------------------------------------------
inline const real* Mesh::coordinates() const
{
  return _geometry.x();
}

//-----------------------------------------------------------------------------
inline uint* Mesh::cells()
{
  return _topology(_topology.dim(), 0)();
}

//-----------------------------------------------------------------------------
inline const uint* Mesh::cells() const
{
  return _topology(_topology.dim(), 0)();
}

//-----------------------------------------------------------------------------
inline uint Mesh::size(uint dim) const
{
  return _topology.size(dim);
}

//-----------------------------------------------------------------------------
inline MeshTopology& Mesh::topology()
{
  return _topology;
}

//-----------------------------------------------------------------------------
inline const MeshTopology& Mesh::topology() const
{
  return _topology;
}

//-----------------------------------------------------------------------------
inline MeshGeometry& Mesh::geometry()
{
  return _geometry;
}

//-----------------------------------------------------------------------------
inline const MeshGeometry& Mesh::geometry() const
{
  return _geometry;
}

//-----------------------------------------------------------------------------
inline MeshDistributedData& Mesh::distdata()
{
  return topology().distdata();
}

//-----------------------------------------------------------------------------
const inline MeshDistributedData& Mesh::distdata() const
{
  return topology().distdata();
}

//-----------------------------------------------------------------------------
inline CellType& Mesh::type()
{
  dolfin_assert(_cell_type);
  return *_cell_type;
}

//-----------------------------------------------------------------------------
inline const CellType& Mesh::type() const
{
  dolfin_assert(_cell_type);
  return *_cell_type;
}

//-----------------------------------------------------------------------------
inline bool Mesh::distributed() const
{
  return _is_distributed;
}

}

#endif
