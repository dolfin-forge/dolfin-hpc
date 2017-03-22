// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2015.
//
// First added:  2008-01-21
// Last changed: 2008-08-12

#include <dolfin/mesh/RefinementManager.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/common/timing.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/RefinementPattern.h>

#include <cstdlib>
#include <map>
#include <set>

#include <time.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
RefinementManager::RefinementManager(Mesh& mesh, Mesh& refined_mesh) :
    mesh_(mesh),
    refined_mesh_(refined_mesh),
    is_distributed_(mesh_.is_distributed()),
    pattern_(&mesh.type()),
    start_offset_(0)
{
  // Initialize internal data structures
  init();
}
//-----------------------------------------------------------------------------
RefinementManager::RefinementManager(Mesh& mesh, Mesh& refined_mesh,
                                     RefinementPattern const& pattern) :
    mesh_(mesh),
    refined_mesh_(refined_mesh),
    is_distributed_(mesh_.is_distributed()),
    pattern_(&pattern),
    start_offset_(0)
{
  // Initialize internal data structures
  init();
}
//-----------------------------------------------------------------------------
RefinementManager::~RefinementManager()
{

}
//-----------------------------------------------------------------------------
void RefinementManager::apply()
{
  map_new_vertices(shared_edge_);
}
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
//-----------------------------------------------------------------------------
void RefinementManager::init()
{
  uint const tdim = mesh_.topology().dim();

  // Generate entity - vertex connectivity if not generated
  for (uint d = 1; d < tdim; ++d)
  {
    if (pattern_->refinement_needs_entities(d))
    {
      mesh_.init(d);
    }
  }

  // Generate facet - cell connectivity if not generated
  if (tdim > 0) mesh_.init(mesh_.type().facet_dim(), tdim);

  // No further step is required in serial
  if (!is_distributed_)
  {
    return;
  }

  // Get number of new vertices from uniform refinement pattern
  uint num_new_vertices = 0;
  for (uint i = 1; i <= tdim; ++i)
  {
    num_new_vertices += mesh_.topology().size(i)
                          * pattern_->num_refined_vertices(i);
  }

  // Assign a safe range for each processor
  uint glb_max = mesh_.topology().global_size(0);
  start_offset_ = 0;
  MPI::offset(num_new_vertices, start_offset_);
  start_offset_ += glb_max;

  // Initialize data structures for interprocess boundary

  cell_forbidden_.init(mesh_, tdim);

  edge_forbidden_.init(mesh_, 1);

  DistributedData& distdata = mesh_.distdata()[0];
  uint const facet_dim = mesh_.type().facet_dim();
  for (SharedIterator it(mesh_.distdata()[facet_dim]); !it.end(); ++it)
  {
    Facet f(mesh_, it.index());
    boundary_cells_.insert(f.entities(tdim)[0]);
  }

  //--- ONLY EDGE BISECTION ---------------------------------------------------

  if (tdim > 1 && pattern_->refinement_needs_entities(1))
  {
    for (SharedIterator it(mesh_.distdata()[1]); !it.end(); ++it)
    {
      Edge e(mesh_, it.index());
      uint const * edge_v = e.entities(0);
      EdgeKey key(distdata.get_global(edge_v[0]),
                  distdata.get_global(edge_v[1]));
      refined_edge_[key] = false;
      edge_keymap_[key] = e.index();
    }
  }

  //--- ONLY EDGE BISECTION ---------------------------------------------------

}
//-----------------------------------------------------------------------------
void RefinementManager::map_new_vertices(Array<uint> shared_edge)
{
  dolfin_assert(refined_mesh_.size(0) == refined_mesh_.geometry().size());

  if(!is_distributed_)
  {
    // Remap local vertices contiguously instead of by supporting entity type.
    // In parallel this is already taken care of by the global renumbering.
    uint const num_local_vertices = refined_mesh_.topology().size(0);
    message(1, "RefinementManager : map %d new vertices in serial",
            num_local_vertices);
    Array<uint> vertex_map(num_local_vertices);
    vertex_map = num_local_vertices;
    uint vertex_count = 0;
    for (CellIterator c(refined_mesh_); !c.end(); ++c)
    {
      for (VertexIterator v(*c); !v.end(); ++v)
      {
        dolfin_assert(v->index() < num_local_vertices);
        if (vertex_map[v->index()] == num_local_vertices)
        {
          vertex_map[v->index()] = vertex_count;
          ++vertex_count;
        }
      }
    }

    if(vertex_count != refined_mesh_.size(0))
    {
      error("RefinementManager : invalid count of vertices '%u' instead of %u",
            vertex_count, refined_mesh_.size(0));
    }

    // Reorder geometry
    dolfin_assert(vertex_count == refined_mesh_.geometry().size());
    refined_mesh_.geometry().remap(vertex_map);

    // Reorder connectivities
    dolfin_assert(vertex_count == refined_mesh_.topology().size(0));
    refined_mesh_.topology().remap(0, vertex_map);

    return;
  }

  message(1, "Map new vertices %u");
  tic();

  DistributedData& olddistdata = mesh_.distdata()[0];
  DistributedData& newdistdata = refined_mesh_.distdata()[0];

  int rank = MPI::rank();
  int pe_size = MPI::size();

  uint num_unass = 0;
  srand((uint) ::time(0) + rank);
  Array<uint> send_buff, send_buff_id;
  std::map<EdgeKey, uint> edge_id;
  std::map<EdgeKey, bool> owns_edge;
  EdgeKey key;
  for (uint i = 0; i < shared_edge.size(); i += 3)
  {

    EdgeKey key(shared_edge[i], shared_edge[i + 1]);
    dolfin_assert(edge_id.count(key) == 0);
    edge_id[key] = (uint) std::rand() + (uint) std::rand() + (uint) rank;
    owns_edge[key] = true;

    send_buff.push_back(olddistdata.get_global(shared_edge[i]));
    send_buff.push_back(olddistdata.get_global(shared_edge[i + 1]));
    send_buff_id.push_back(edge_id[key]);
  }

  // Assign ownership of shared edges
  MPI_Status status;
  uint src, dest;
  int max_un, num_un, max_id, num_id, recv_count;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  num_id = send_buff_id.size();
  MPI_Allreduce(&num_id, &max_id, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  uint *recv_buff = new uint[max_un];
  uint *recv_buff_id = new uint[max_id];
  for (int j = 1; j < pe_size; j++)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff_id[0], num_id, MPI_UNSIGNED, dest, 1, recv_buff_id,
                 max_id, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, recv_buff,
                 max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (uint i = 0; i < (uint) recv_count; i += 2)
    {
      // Check if I have the vertices
      if (olddistdata.has_global(recv_buff[i])
          && olddistdata.has_global(recv_buff[i + 1]))
      {

        // Generate edge key
        EdgeKey key(olddistdata.get_local(recv_buff[i]),
                    olddistdata.get_local(recv_buff[i + 1]));

        // Check if I have the corresponding edge
        if (edge_id.count(key))
        {
          newdistdata.set_shared_adj(new_edge_vertex_[key], src);
          if (recv_buff_id[i >> 1] < edge_id[key]
              || (recv_buff_id[i >> 1] == edge_id[key]
                  && status.MPI_SOURCE < rank))
          {
            owns_edge[key] = false;
            new_edge_global_.erase(key);
            edge_id.erase(key);
            num_unass++;
          }
        }
      }
    }
  }

  //Exchange assigned global numbers
  Array<uint> global_buff;
  uint index;
  for (int j = 1; j < pe_size; j++)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, recv_buff,
                 max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (uint i = 0; i < (uint) recv_count; i += 2)
    {
      if (olddistdata.has_global(recv_buff[i])
          && olddistdata.has_global(recv_buff[i + 1]))
      {

        EdgeKey key(olddistdata.get_local(recv_buff[i]),
                    olddistdata.get_local(recv_buff[i + 1]));

        if (owns_edge[key])
        {
          global_buff.push_back(i);
          global_buff.push_back(new_edge_global_[key]);
          //newdistdata.set_shared_adj(new_edge_vertex_[key], src);
        }
      }
    }

    MPI_Sendrecv(&global_buff[0], global_buff.size(), MPI_UNSIGNED, src, 2,
                 recv_buff, max_un, MPI_UNSIGNED, dest, 2, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (uint i = 0; i < (uint) recv_count; i += 2)
    {
      index = shared_edge[(recv_buff[i] >> 1) * 3 + 2];
      newdistdata.set_map(index, recv_buff[i + 1], true);
      newdistdata.set_ghost(index, status.MPI_SOURCE);
      num_unass--;
    }
    global_buff.clear();
  }

  dolfin_assert(num_unass == 0);

  delete[] recv_buff;
  delete[] recv_buff_id;

  // Finalize and renumber globally so that vertices are indexed correctly after
  // apply the refinement
  newdistdata.finalize();

  uint num_shared_edges = 0;
  MPI::allReduceSum(mesh_.topology().num_shared(1), num_shared_edges);

  dolfin_assert(newdistdata.global_size() == mesh_.global_size(0) + num_shared_edges);

  newdistdata.renumber_global();

  tocd(1);
}
//-----------------------------------------------------------------------------
void RefinementManager::mark_localboundary(MeshFunction<bool>& cell_marker,
                                           uint& num_new_vertices,
                                           uint& num_new_cells)
{
  DistributedData& olddistdata = mesh_.distdata()[0];

  uint rank = MPI::rank();
  uint pe_size = MPI::size();
  srand((uint) ::time(0) + rank);

  Array<uint> send_buff;
  uint edge[2];
  EdgeKey key;
  // Set of cell indices with longest edge
  std::map<EdgeKey, uint> num_ref, removed, edge_id;
  _set<uint> cell_forbidden_edges;
  std::set<EdgeKey> forbidden_propagation;
  _map<uint, uint> edge_vote;

  num_new_cells = 0;
  num_new_vertices = 0;
  real max, l;
  uint index = 0;

  // Process cells between processors
  _set<uint>::iterator bc;

  // Reset forbidden edges and cells
  edge_forbidden_ = false;
  cell_forbidden_ = false;

  for (bc = boundary_cells_.begin(); bc != boundary_cells_.end(); bc++)
  {
    Cell c(mesh_, *bc);
    if (cell_marker.get(c) && !cell_forbidden_.get(c))
    {
      max = 0.0;
      for (EdgeIterator e(c); !e.end(); ++e)
      {
        if (edge_forbidden_.get(*e))
          continue;
        edge_vote[e->index()] = 0;
        l = e->length();
        if (max < l)
        {
          max = l;
          index = e->index();
        }
      }

      if (max > 0.0)
      {
        Edge longest_edge(mesh_, index);
        if (on_boundary(longest_edge))
        {
          const uint *edge_v = longest_edge.entities(0);
          edge_vote[longest_edge.index()] = (uint) std::rand();
          send_buff.push_back(olddistdata.get_global(edge_v[0]));
          send_buff.push_back(olddistdata.get_global(edge_v[1]));
          send_buff.push_back(edge_vote[longest_edge.index()]);
          for (CellIterator nc(longest_edge); !nc.end(); ++nc)
          {
            cell_forbidden_.set(*nc, true);
            for (EdgeIterator e(*nc); !e.end(); ++e)
            {
              edge_forbidden_.set(*e, true);
            }
          }
        }
      }
    }
  }

  // Decide ownership of refinement
  MPI_Status status;
  uint src, dest;
  int max_un, num_un, recv_count;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  uint *recv_buff = new uint[max_un];
  Array<uint> forbidden;

  for (uint j = 1; j < pe_size; j++)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, dest, 1,
                 recv_buff, max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i += 3)
    {
      EdgeKey key(recv_buff[i], recv_buff[i + 1]);

      // If rank has edge
      if (edge_keymap_.count(key))
      {
        Edge e(mesh_, edge_keymap_[key]);
        edge_vote[e.index()] += recv_buff[i + 2];
      }
    }
  }

  send_buff.clear();
  cell_forbidden_ = false;
  edge_forbidden_ = false;

  for (bc = boundary_cells_.begin(); bc != boundary_cells_.end(); bc++)
  {
    Cell c(mesh_, *bc);
    if (cell_marker.get(c) && !cell_forbidden_.get(c))
    {
      max = 0.0;
      for (EdgeIterator e(c); !e.end(); ++e)
      {
        if (edge_forbidden_.get(*e))
          continue;
        l = (real) edge_vote[e->index()];
        if (max < l)
        {
          max = l;
          index = e->index();
        }
      }

      if (max > 0.0)
      {
        Edge longest_edge(mesh_, index);
        if (on_boundary(longest_edge))
        {
          const uint *edge_v = longest_edge.entities(0);
          edge[0] = olddistdata.get_global(edge_v[0]);
          edge[1] = olddistdata.get_global(edge_v[1]);
          EdgeKey key(edge[0], edge[1]);
          refined_edge_[key] = true;
          num_new_vertices++;
          send_buff.push_back(edge[0]);
          send_buff.push_back(edge[1]);

          for (CellIterator nc(longest_edge); !nc.end(); ++nc)
          {
            cell_forbidden_.set(*nc, true);
            num_ref[key]++;
            num_new_cells++;
            cell_refedge_[nc->index()] = longest_edge.index();
            for (EdgeIterator e(*nc); !e.end(); ++e)
            {
              edge_forbidden_.set(*e, true);
              cell_forbidden_edges.insert(e->index());
            }
          }
        }
      }
    }
  }

  Array<uint> terminated;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  delete[] recv_buff;
  recv_buff = new uint[max_un];

  for (uint j = 1; j < pe_size; j++)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, dest, 1,
                 recv_buff, max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i += 2)
    {
      EdgeKey key(recv_buff[i], recv_buff[i + 1]);
      // If rank has edge
      if (edge_keymap_.count(key))
      {
        if (!refined_edge_[key])
        {
          Edge e(mesh_, edge_keymap_[key]);
          if (cell_forbidden_edges.count(e.index()) == 0)
          {
            num_new_vertices++;
            refined_edge_[key] = true;
            for (CellIterator c(e); !c.end(); ++c)
            {
              num_new_cells++;
              num_ref[key]++;
              cell_forbidden_.set(*c, true);
              cell_refedge_[c->index()] = e.index();
              for (EdgeIterator ce(*c); !ce.end(); ++ce)
              {
                edge_forbidden_.set(*ce, true);
                cell_forbidden_edges.insert(ce->index());
              }
            }
          }
          else
          {
            forbidden.push_back(recv_buff[i]);
            forbidden.push_back(recv_buff[i + 1]);
          }
        }
      }
    }
    MPI_Sendrecv(&forbidden[0], forbidden.size(), MPI_UNSIGNED, src, 2,
                 recv_buff, max_un, MPI_UNSIGNED, dest, 2, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i += 2)
    {
      EdgeKey key(recv_buff[i], recv_buff[i + 1]);
      if (removed.count(key) == 0)
      {
        num_new_vertices--;
        removed[key] = 1;
        Edge e(mesh_, edge_keymap_[key]);
        for (CellIterator c(e); !c.end(); ++c)
        {
          if (cell_forbidden_.get(*c))
            num_new_cells--;
          cell_forbidden_.set(*c, false);
          cell_refedge_.erase(c->index());
        }
        terminated.push_back(recv_buff[i]);
        terminated.push_back(recv_buff[i + 1]);
      }
    }
    forbidden.clear();
  }

  num_un = terminated.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  delete[] recv_buff;
  recv_buff = new uint[max_un];

  for (uint j = 1; j < pe_size; j++)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&terminated[0], terminated.size(), MPI_UNSIGNED, dest, 1,
                 recv_buff, max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i += 2)
    {
      EdgeKey key(recv_buff[i], recv_buff[i + 1]);
      // If rank has edge
      if (edge_keymap_.count(key))
      {
        // Remove refinement if edge is refined with terminated refinement
        if (refined_edge_[key])
        {
          if (removed.count(key) == 0)
          {
            num_new_vertices--;
            removed[key] = 1;
            Edge e(mesh_, edge_keymap_[key]);
            for (CellIterator c(e); !c.end(); ++c)
            {
              if (cell_forbidden_.get(*c))
                num_new_cells--;
              cell_forbidden_.set(*c, false);
              cell_refedge_.erase(c->index());
            }
            refined_edge_.erase(key);
          }
        }
      }
    }
  }

  // Mark unrefined cells shared edges as forbidden
  for (bc = boundary_cells_.begin(); bc != boundary_cells_.end(); bc++)
  {
    Cell c(mesh_, *bc);
    if (!cell_forbidden_.get(c))
    {
      cell_refedge_.erase(c.index());
      for (EdgeIterator e(c); !e.end(); ++e)
      {
        if (on_boundary(*e))
          edge_forbidden_.set(*e, true);
        else
        {
          bool ok_to_remove = true;
          for (CellIterator ec(*e); !ec.end(); ++ec)
          {
            if (cell_forbidden_.get(*ec))
            {
              ok_to_remove = false;
              break;
            }
          }
          if (ok_to_remove)
            edge_forbidden_.set(*e, false);
        }
      }
    }
  }

  delete[] recv_buff;
}
//-----------------------------------------------------------------------------
#else
void RefinementManager::init()
{
}
//-----------------------------------------------------------------------------
void RefinementManager::mark_localboundary( MeshFunction<bool>& cell_marker,
                                            uint& num_new_vertices,
                                            uint& num_new_cells)
{
}
//-----------------------------------------------------------------------------
void RefinementManager::map_new_vertices(Array<uint> shared_edge)
{
}
//-----------------------------------------------------------------------------
#endif

//-----------------------------------------------------------------------------

}
