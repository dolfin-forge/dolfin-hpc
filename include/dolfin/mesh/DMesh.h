// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_D_MESH_H
#define __DOLFIN_D_MESH_H

#include <dolfin/common/types.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/EdgeKey.h>
#include <dolfin/mesh/entities/Cell.h>

#include <list>
#include <vector>

namespace dolfin
{

class CellType;
struct DCell;
struct DVertex;
class Mesh;
template < class T, class E, size_t N >
struct MeshValues;

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
  DMesh( Mesh & mesh );

  /// Destructor
  ~DMesh();

  /// Export to a regular mesh
  void exp( Mesh & mesh );

  /// Bisect marked cells
  void bisectMarked( MeshValues< bool, Cell > const & marked_ids );

private:
  /// Edge data structure for propagation
  typedef struct __edge__
  {
    size_t mv;    //< global index of midpoint vertex
    size_t v1;    //< global index of endpoint
    size_t v2;    //< global index of endpoint
    size_t owner; //< rank of owner
  } prop_edge;

  static auto less_pair_comp( std::pair< size_t, prop_edge > const & x,
                              std::pair< size_t, prop_edge > const & y )
    -> bool;

  /// Pair datatype for propagation
  using Propagation = std::pair< size_t, prop_edge >;

  /// Add a new vertex
  void add_vertex( DVertex * v );

  /// Add a new cell with vertices vs and inside existing cell parent_id
  void add_cell( DCell * c, std::vector< DVertex * > vs, int parent_id );

  /// Remove a vertex
  ///
  /// Entity is just marked as deleted, but not yet erased
  void removeVertex( DVertex * v );

  /// Remove a cell
  ///
  /// Entity is just marked as deleted, but not yet erased
  void removeCell( DCell * c );

  /// Renumber mesh entities locally
  ///
  /// An optional mapping between old and new indices is generated. The
  /// std::vectors have to have the size of the old numbering
  void number( std::vector< int > * old2new_cells    = nullptr,
               std::vector< int > * old2new_vertices = nullptr );

  /// Renumber global indicies to fit within the range of an unsigned int
  void renumber_glb( _map< size_t, size_t > & new_global );

  /// Bisect cell dcell
  ///
  /// The edge for the bisection is given by hv0 and hv1 and hangv is the
  /// hanging node of the opposite cell
  void bisect( DCell * dcell, DVertex * hangv, DVertex * hv0, DVertex * hv1 );

  /// Get opposite cell with respect to vertices v1 and v2
  auto opposite( DCell * dcell, DVertex * v1, DVertex * v2 ) -> DCell *;

  /// Propagate refinement
  void propagate_refinement( Mesh &                       mesh,
                             std::vector< Propagation > & propagation,
                             bool &                       empty );

  /// Naive refinement propagation with pairwise communication
  void propagate_naive( Mesh &                       mesh,
                        std::vector< Propagation > & propagation,
                        bool &                       empty );

  /// Refinement propagation within hypercube
  void propagate_hypercube( Mesh &                       mesh,
                            std::vector< Propagation > & propagation,
                            bool &                       empty );

  /// Vertices contained in the mesh
  using VertexSet = _ordered_set< DVertex * >;
  VertexSet vertices;

  /// Cells contained in the mesh
  using CellList = std::list< DCell * >;
  CellList cells;

  /// Propagation buffer
  std::vector< Propagation > propagate;

  /// Map between global number of boundary vertex to vertex
  using BoundaryVertices = _map< size_t, DVertex * >;
  BoundaryVertices bc_dvs;

  /// Refined edges
  using RefinedEdges = _map< EdgeKey< size_t >, DVertex * >;
  RefinedEdges ref_edge;

  /// Mesh
  Mesh & mesh_;

  /// CellType of mesh
  CellType const * const ctype_;
  Space const * const    space_;

  /// Maximum global index of vertices
  /// Implemented as number of vertices in the *global* mesh
  size_t glb_max_;

  /// enumeration salt for bisect
  size_t salt_;

  /// Count of deleted
  size_t cdeleted_;
  size_t vdeleted_;
};
//-----------------------------------------------------------------------------
/// Comparison operator for index/value pairs
inline auto DMesh::less_pair_comp( std::pair< size_t, prop_edge > const & x,
                                   std::pair< size_t, prop_edge > const & y )
  -> bool
{
  return x.first < y.first;
}

} /* namespace dolfin */

#endif
