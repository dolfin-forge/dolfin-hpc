// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2013.
// Modified by Balthasar Reuter, 2013
// Modified by Aurelien Larcher, 2015
//

#ifndef __DOLFIN_D_MESH_H
#define __DOLFIN_D_MESH_H

#include <dolfin/common/types.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/EdgeKey.h>

#include <vector>
#include <list>
#include <set>
#include <map>

namespace dolfin
{

template<class T> class Array;
class CellType;
class DCell;
class DVertex;
class Mesh;
template<class T, class E, uint N> class MeshValues;

/// Dynamic mesh class for on-the-fly changes to the mesh. It is used by the
/// recursive RivaraRefinement and the EdgeCollapse-MeshCoarsening.
///
/// It provides import- and export-routines from and to the regular Mesh class
/// of DOLFIN but cannot be used directly for MeshFunction or assembly
/// routines.
///
/// For mesh refinement it has methods for cell-based bisection.
///
class DMesh
{
public:

  /// Constructor on a mesh (preferred)
  DMesh(Mesh& mesh);

  /// Destructor
  ~DMesh();

  /// Export to a regular mesh
  void exp(Mesh& mesh);

  /// Bisect marked cells
  void bisectMarked(MeshValues<bool, Cell> const& marked_ids);

private:

  /// Edge data structure for propagation
  typedef struct __edge__
  {
    uint mv;    //< global index of midpoint vertex
    uint v1;    //< global index of endpoint
    uint v2;    //< global index of endpoint
    uint owner;  //< rank of owner
  } prop_edge;

  /// Pair datatype for propagation
  typedef std::pair<uint, prop_edge> Propagation;

  /// Erase removed entities from datastructures
  ///
  /// removeVertex() and removeCell() only mark entitites for deletion but are
  /// not actually erased
  void eraseRemovedEntities();

  /// Find Vertex by its local id
  DVertex* getVertex(int local_id);

  /// Find Cell by its local id
  DCell* getCell(int local_id);

  /// Add a new vertex
  void add_vertex(DVertex* v);

  /// Add a new cell with vertices vs and inside existing cell parent_id
  void add_cell(DCell* c, std::vector<DVertex*> vs, int parent_id);

  /// Remove a vertex
  ///
  /// Entity is just marked as deleted, but not yet erased
  void removeVertex(DVertex* v);

  /// Remove a cell
  ///
  /// Entity is just marked as deleted, but not yet erased
  void removeCell(DCell* c);

  /// Export to a regular mesh but keep numbering in the DMesh
  ///
  /// An optional mapping between old and new indices is generated. The Arrays
  /// have to have at least the size of the old numbering (as it was at the
  /// time of the import).
  void expKeepNumbering(Mesh& mesh, Array<int> * old2new_cells = 0,
                        Array<int> * old2new_vertices = 0);

  /// Renumber mesh entities locally
  ///
  /// An optional mapping between old and new indices is generated. The Arrays
  /// have to have the size of the old numbering
  void number(Array<int> *old2new_cells = 0, Array<int> *old2new_vertices = 0);

  /// Bisect cell dcell
  ///
  /// The edge for the bisection is given by hv0 and hv1 and hangv is the
  /// hanging node of the opposite cell
  void bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1);

  /// Get opposite cell with respect to vertices v1 and v2
  DCell* opposite(DCell* dcell, DVertex* v1, DVertex* v2);

  /// Propagate refinement
  ///
  /// TODO: what are the arguments???
  inline void propagate_refinement(std::vector<Propagation>& propagated,
                                   bool& empty)
  {
    if (PE::size() == 1) return;
    if (PE::size() & (PE::size() - 1))
    {
      propagate_naive(propagated, empty);
    }
    else
    {
      propagate_hypercube(propagated, empty);
    }
  }

  /// Naive refinement propagation with pairwise communication
  void propagate_naive(std::vector<Propagation>& propagated, bool& empty);

  /// Refinement propagation within hypercube
  void propagate_hypercube(std::vector<Propagation>& propagated, bool& empty);

  /// Vertices contained in the mesh
  typedef std::set<DVertex *> VertexSet;
  VertexSet vertices;

  /// Cells contained in the mesh
  typedef std::list<DCell *> CellList;
  std::list<DCell *> cells;

  /// Propagation buffer
  std::vector<Propagation> propagate;

  /// Map between global number of boundary vertex to vertex
  typedef _map<uint, DVertex*> BoundaryVertices;
  BoundaryVertices bc_dvs;

  /// Refined edges
  typedef _map<EdgeKey, DVertex*> RefinedEdges;
  RefinedEdges ref_edge;

  /// Comparison operator for index/value pairs
  struct less_pair : public std::binary_function<std::pair<uint, prop_edge>,
  std::pair<uint, prop_edge>, bool>
  {
    bool operator()(std::pair<uint, prop_edge> x, std::pair<uint, prop_edge> y)
    {
      return x.first < y.first;
    }
  };

  /// Mesh
  Mesh& mesh_;

  /// CellType of mesh
  CellType const * const ctype_;
  Space    const * const space_;

  /// Shared edges
  typedef std::map<EdgeKey, uint>  SharedEdges;
  typedef std::pair<EdgeKey, uint> SharedEdgeItem;
  SharedEdges * shared_edges_;

  /// Maximum global index of vertices
  /// Implemented as number of vertices in the *global* mesh
  uint glb_max_;

  /// enumeration salt for bisect
  uint salt_;

};

} /* namespace dolfin */

#endif
