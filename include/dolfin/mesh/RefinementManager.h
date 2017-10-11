// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2015.
//
// First added:  2008-01-21
// Last changed: 2008-07-23

#ifndef __DOLFIN_REFINEMENT_MANAGER_H
#define __DOLFIN_REFINEMENT_MANAGER_H

#include <dolfin/common/Array.h>
#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/EdgeKey.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/mesh/FaceKey.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/MeshValues.h>

#include <map>
#include <set>

namespace dolfin
{

class Mesh;
class RefinementPattern;

/**
 *  @class  RefinementManager
 *
 *  @brief
 *
 */

class RefinementManager
{

public:

  /// Refine mesh with default refinement pattern
  RefinementManager(Mesh& mesh, Mesh& refined_mesh);

  /// Refine mesh with provided refinement pattern
  RefinementManager(Mesh& mesh, Mesh& refined_mesh,
                    RefinementPattern const& pattern);

  /// Destructor
  ~RefinementManager();

  /// Return cell refinement pattern
  RefinementPattern const& pattern() const;

  /// Add original vertex with given index
  void add(Vertex& v, uint index);

  /// Add edge-based vertex with given index
  void add(Edge& e, uint index);

  /// Add face-based vertex with given index
  void add(Face& f, uint index);

  /// Check if the edge lies on the interprocess boundary
  bool on_boundary(Edge& e);

  /// Check if the face lies on the interprocess boundary
  bool on_boundary(Face& f);

  /// Apply refinement and finalize numbering
  void apply();

  ///--- EDGE BISECTION ONLY  -------------------------------------------------

  void mark_localboundary(MeshValues<bool, Cell>& cell_marker,
                          uint& num_new_vertices, uint& num_new_cells);

  /// Check if the cell has received a propagation
  bool forbidden_cell(Cell& cell);

  /// Check if the edge are forbidden from propagated refinement
  bool forbidden_edge(Edge& edge);

  ///
  uint edge_refined(Cell& cell);

  ///--- EDGE BISECTION ONLY  -------------------------------------------------

private:

  // Initialize mesh connectivities and internal data structures
  void init();

  // Map global numbers to unassigned shared vertices
  void map_new_vertices(Array<uint> shared_edge);

  //
  Mesh& mesh_;
  Mesh& refined_mesh_;
  bool const is_distributed_;
  RefinementPattern const * const pattern_;
  uint start_offset_;

  // Edge vertices
  Array<uint> shared_edge_;
  std::map< EdgeKey, uint> new_edge_global_;
  std::map< EdgeKey, uint> new_edge_vertex_;

  // Face vertices
  Array<uint> shared_face_;
  std::map< FaceKey, uint> new_face_global_;
  std::map< FaceKey, uint> new_face_vertex_;

  //--- ONLY EDGE BISECTION ---------------------------------------------------

  _set<uint> boundary_cells_;
  _map<uint, uint> cell_refedge_;
  MeshValues<bool, Cell> * cell_forbidden_;
  MeshValues<bool, Edge> * edge_forbidden_;
  std::map<EdgeKey, uint> edge_keymap_;
  std::map<EdgeKey, bool> refined_edge_;

  //--- ONLY EDGE BISECTION ---------------------------------------------------

};

//--- INLINES -----------------------------------------------------------------

inline RefinementPattern const& RefinementManager::pattern() const
{
  return (*pattern_);
}

//-----------------------------------------------------------------------------
inline void RefinementManager::add(Vertex& v, uint index)
{
  if (is_distributed_)
  {
    refined_mesh_.distdata()[0].set_map(index, v.global_index());
    if(v.is_ghost())
    {
      refined_mesh_.distdata()[0].set_ghost(index, v.owner());
    }
    else if (v.is_shared())
    {
      refined_mesh_.distdata()[0].set_shared(index);
    }
  }
}
//-----------------------------------------------------------------------------
inline void RefinementManager::add(Edge& e, uint index)
{
  if (is_distributed_)
  {
    if(on_boundary(e))
    {
      // Store edge key in shared list
      EdgeKey key(e.entities(0)[0], e.entities(0)[1]);
      new_edge_global_[key] = start_offset_;
      new_edge_vertex_[key] = index;

      // Buffer edge information for mapping phase
      shared_edge_.push_back(key.first);
      shared_edge_.push_back(key.second);
      shared_edge_.push_back(index);
    }
    refined_mesh_.distdata()[0].set_map(index, start_offset_++);
  }
}
//-----------------------------------------------------------------------------
inline void RefinementManager::add(Face& f, uint index)
{
  if (is_distributed_)
  {
    if(on_boundary(f))
    {
      // Store edge key in shared list
      FaceKey key(f);
      new_face_global_[key] = start_offset_;
      new_face_vertex_[key] = index;

      // Buffer face information for mapping phase
      for (uint i = 0; i < f.num_entities(0); ++i)
      {
        shared_face_.push_back(f.entities(0)[i]);
      }
      shared_face_.push_back(index);
    }
    refined_mesh_.distdata()[0].set_map(index, start_offset_++);
  }
}
//-----------------------------------------------------------------------------
inline bool RefinementManager::on_boundary(Edge& e)
{
  return e.is_shared();
}
//-----------------------------------------------------------------------------
inline bool RefinementManager::on_boundary(Face& f)
{
  return f.is_shared();
}
//-----------------------------------------------------------------------------
inline bool RefinementManager::forbidden_cell(Cell& cell)
{
  return (cell_forbidden_ ? (*cell_forbidden_)(cell) : false);
}

//-----------------------------------------------------------------------------
inline bool RefinementManager::forbidden_edge(Edge& edge)
{
  return (edge_forbidden_ ? (*edge_forbidden_)(edge) : false);
}

//-----------------------------------------------------------------------------
inline uint RefinementManager::edge_refined(Cell& cell)
{
  dolfin_assert(cell_refedge_.count(cell.index()));
  return cell_refedge_[cell.index()];
}

}

#endif
