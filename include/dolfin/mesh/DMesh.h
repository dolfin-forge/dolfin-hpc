// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2019.
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
struct DCell;
struct DVertex;
class Mesh;
template<class T, class E, uint N> struct MeshValues;

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
    long mv;    //< global index of midpoint vertex
    long v1;    //< global index of endpoint
    long v2;    //< global index of endpoint
    uint owner;	//< rank of owner
  } prop_edge;

  /// Pair datatype for propagation
  typedef std::pair<uint, prop_edge> Propagation;

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

  /// Renumber mesh entities locally
  ///
  /// An optional mapping between old and new indices is generated. The Arrays
  /// have to have the size of the old numbering
  void number(Array<int> *old2new_cells = NULL,
              Array<int> *old2new_vertices = NULL);

  /// Renumber global indicies to fit within the range of an unsigned int
  void renumber_glb(_map<long, uint>& new_global);

  /// Bisect cell dcell
  ///
  /// The edge for the bisection is given by hv0 and hv1 and hangv is the
  /// hanging node of the opposite cell
  void bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1);

  /// Get opposite cell with respect to vertices v1 and v2
  DCell* opposite(DCell* dcell, DVertex* v1, DVertex* v2);

  /// Propagate refinement
  void propagate_refinement(Mesh& mesh,
                            Array<Propagation>& propagation, bool& empty);

  /// Naive refinement propagation with pairwise communication
  void propagate_naive(Mesh& mesh,
                       Array<Propagation>& propagation, bool& empty);

  /// Refinement propagation within hypercube
  void propagate_hypercube(Mesh& mesh,
                           Array<Propagation>& propagation, bool& empty);

  /// Vertices contained in the mesh
  typedef std::set<DVertex *> VertexSet;
  VertexSet vertices;

  /// Cells contained in the mesh
  typedef std::list<DCell *> CellList;
  CellList cells;

  /// Propagation buffer
  Array<Propagation> propagate;

  /// Map between global number of boundary vertex to vertex
  typedef _map<long, DVertex*> BoundaryVertices;
  BoundaryVertices bc_dvs;

  /// Refined edges
  typedef _map<EdgeKey<long>, DVertex*> RefinedEdges;
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

  /// Maximum global index of vertices
  /// Implemented as number of vertices in the *global* mesh
  long glb_max_;

  /// enumeration salt for bisect
  long salt_;

  /// Count of deleted
  uint cdeleted_;
  uint vdeleted_;

};

} /* namespace dolfin */

#endif
