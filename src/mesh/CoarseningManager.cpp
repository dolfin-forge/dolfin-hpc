// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-03-25
// Last changed: 2013-06-01

#include <dolfin/mesh/CoarseningManager.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/LoadBalancer.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/mesh/MeshFunctionConverter.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/DCell.h>
#include <dolfin/mesh/DVertex.h>
#include <dolfin/parameter/ParameterSystem.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
CoarseningManager::CoarseningManager() :
    dmesh_(NULL)
{
  // do nothing
}
//-----------------------------------------------------------------------------
CoarseningManager::CoarseningManager(Mesh& mesh,
                                     MeshFunction<bool>& cell_marker,
                                     bool coarsen_boundary) :
    dmesh_(NULL)
{
  init(mesh, cell_marker, coarsen_boundary);
}
//-----------------------------------------------------------------------------
CoarseningManager::~CoarseningManager()
{
  delete dmesh_;
}
//-----------------------------------------------------------------------------
void CoarseningManager::init(Mesh& mesh, MeshFunction<bool>& cell_marker,
                             bool coarsen_boundary)
{
  dolfin_assert(&mesh == &(cell_marker.mesh()));

  num_migrated_cells_ = 1; // choosing initally > 0 s.t. coarsening won't exit
                           // in first migration if all marked cells are requested

  if (ParameterSystem::parameters.defined("coarsening quality threshold")) quality_threshold_ =
      dolfin_get("coarsening quality threshold");
  else quality_threshold_ = 1e-3;

  if (MPI::processNumber() == 0) message("Quality threshold: %d",
                                         quality_threshold_);

  // attempt count is for sake of simplicity in the migration phase always one
  // larger than the number of actually performed attempts
  MeshFunction<uint> attempt_count;
  MeshFunctionConverter::cast(cell_marker, attempt_count);

  initCommon(mesh, &attempt_count);

  // find independent set
  findIndependentSet(mesh, coarsen_boundary);

  dolfin_assert(forbidden_vertices_.size() == cell_marker.mesh().size(0));
}
//-----------------------------------------------------------------------------
void CoarseningManager::initCommon(Mesh& mesh,
                                   MeshFunction<uint> *attempt_count)
{
  dolfin_assert(&mesh == &(attempt_count->mesh()));

  /// Delete old dmesh and import new mesh into dmesh
  delete dmesh_;
  dmesh_ = new DMesh;
  dmesh_->imp(mesh);

  orig_num_cells_ = mesh.num_cells();
  orig_num_vertices_ = mesh.size(0);

  // build list of cells to coarsen
  findCellsToCoarsen(attempt_count);

  // reset request lists
  cells_to_request_.clear();
  vertices_to_request_.clear();
}
//-----------------------------------------------------------------------------
void CoarseningManager::findIndependentSet(Mesh& mesh, bool coarsen_boundary)
{
  // set to false on whole domain
  forbidden_vertices_.resize(mesh.size(0));
  forbidden_vertices_ = false;

  uint num_independent_vertices(0);
  // if boundary coarsening is forbidden: put boundary vertices into set first
  if (!coarsen_boundary)
  {
    BoundaryMesh boundary(mesh, BoundaryMesh::exterior);

    if (boundary.size(0) > 0)
    {
      // add boundary vertices to set
      for (VertexIterator v_it(boundary); !v_it.end(); ++v_it)
      {
        forbidden_vertices_.at(boundary.vertex_index(v_it->index())) = true;
        if (!mesh.distdata()[0].is_ghost(boundary.vertex_index(v_it->index())))
        {
          ++num_independent_vertices;
        }
      }
    }
  }

  // iterate over remaining vertices
  for (VertexIterator v_it(mesh); !v_it.end(); ++v_it)
  {
    if (!forbidden_vertices_[v_it->index()])
    {
      bool independent = isIndependentVertex(*v_it);
      forbidden_vertices_.at(v_it->index()) = independent;
      if (independent && !mesh.distdata()[0].is_ghost(v_it->index())) ++num_independent_vertices;
    }
  }
#ifdef HAVE_MPI
  uint global_num_independent_vertices;
  MPI_Reduce(&num_independent_vertices, &global_num_independent_vertices, 1,
             MPI_UNSIGNED, MPI_SUM, 0, MPI::DOLFIN_COMM);
#endif
}
//-----------------------------------------------------------------------------
bool CoarseningManager::isIndependentVertex(Vertex& v)
{
  // Check all vertices connected to that vertex
  for (VertexIterator v_it(v); !v_it.end(); ++v_it)
  {
    if (v_it->index() == v.index()) continue;

    if (forbidden_vertices_.at(v_it->index())) return false;
  }

  // no neighboring point in set: independent vertex
  return true;
}
//-----------------------------------------------------------------------------
void CoarseningManager::findCellsToCoarsen(MeshFunction<uint> * attempt_count)
{
  cells_to_coarsen_.clear();

  for (std::list<DCell *>::iterator c_it(dmesh_->cells.begin());
      c_it != dmesh_->cells.end(); ++c_it)
  {
    DCell * dc = *c_it;
    if (attempt_count->get(dc->id) > 0) cells_to_coarsen_.push_back(
        std::make_pair(dc, attempt_count->get(dc->id)));
  }
}
//-----------------------------------------------------------------------------
void CoarseningManager::removeErasedCellsFromCoarseningList()
{
  // remove deleted entities from coarsening list
  for (List<std::pair<DCell *, uint> >::iterator it(cells_to_coarsen_.begin());
      it != cells_to_coarsen_.end();)
  {
    DCell * dc = it->first;
    if (dc->deleted) it = cells_to_coarsen_.erase(it);
    else ++it;
  }
}
//-----------------------------------------------------------------------------
void CoarseningManager::buildMFArrays(
    Mesh& mesh,
    Array<int>& old2new_cells,
    Array<int>& old2new_vertices,
    Array<std::pair<MeshFunction<uint>*, MeshFunction<uint>*> >& cell_functions,
    Array<std::pair<MeshFunction<real>*, MeshFunction<real>*> >& vertex_functions)
{
  MeshFunction<real> * forbidden_vertices = new MeshFunction<real>();
  MeshFunction<real> * forbidden_vertices_new = new MeshFunction<real>();
  MeshFunction<uint> * attempts = new MeshFunction<uint>();
  MeshFunction<uint> * attempts_new = new MeshFunction<uint>();

  if (mesh.num_cells() != 0)
  {
    // Prepare forbidden vertices for exchange
    forbidden_vertices->init(mesh, 0);
    *forbidden_vertices = 0.0;
    for (uint old_idx(0); old_idx < forbidden_vertices_.size(); ++old_idx)
    {
      int new_idx = old2new_vertices[old_idx];
      if (new_idx < 0) continue;
      forbidden_vertices->set(new_idx, double(isForbiddenVertex(old_idx)));
    }

    // Prepare cell_marker and attempt count for exchange
    attempts->init(mesh, mesh.topology().dim());
    *attempts = 0u;
    for (List<std::pair<DCell *, uint> >::iterator it(
        cells_to_coarsen_.begin()); it != cells_to_coarsen_.end(); ++it)
    {
      DCell * dc = it->first;
      int new_idx = old2new_cells[dc->id];
      if (new_idx < 0) continue;
      attempts->set(new_idx, it->second);
      dolfin_assert(it->second > 0);
    }
  }

  vertex_functions.push_back(
      std::make_pair(forbidden_vertices, forbidden_vertices_new));
  cell_functions.push_back(std::make_pair(attempts, attempts_new));
}
//-----------------------------------------------------------------------------
void CoarseningManager::cleanupMFArrays(
    Array<std::pair<MeshFunction<uint>*, MeshFunction<uint>*> >& cell_functions,
    Array<std::pair<MeshFunction<real>*, MeshFunction<real>*> >& vertex_functions)
{
  for (Array<std::pair<MeshFunction<uint>*, MeshFunction<uint>*> >::iterator it(
      cell_functions.begin()); it != cell_functions.end(); ++it)
  {
    delete it->first;
    delete it->second;
  }

  for (Array<std::pair<MeshFunction<real>*, MeshFunction<real>*> >::iterator it(
      vertex_functions.begin()); it != vertex_functions.end(); ++it)
  {
    delete it->first;
    delete it->second;
  }
}
//-----------------------------------------------------------------------------
void CoarseningManager::updateIndependentSet(
    Mesh& mesh, MeshFunction<real>& forbidden_vertices_new)
{
  forbidden_vertices_.resize(mesh.size(0));
  forbidden_vertices_ = false;
  for (VertexIterator v_it(mesh); !v_it.end(); ++v_it)
    if (forbidden_vertices_new.get(v_it->index()) > 0.5) forbidden_vertices_[v_it->index()] =
        true;
}
//-----------------------------------------------------------------------------
void CoarseningManager::exchangeRequests(Mesh& mesh, Array<int>& old2new_cells,
                                         Array<int>& old2new_vertices,
                                         uint max_num_requested_vertices,
                                         MeshFunction<uint> *& partitions)
{
#ifdef HAVE_MPI
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();

  // List of vertices to request from owner
  Array<uint> *send_list_requests = new Array<uint> [pe_size];

  // Received requests
  uint * recv_buff_requests = new uint[2 * max_num_requested_vertices];

  // Map of requested vertices (that this process owns) and requesting processes
  std::map<uint, uint> requested_vertices;

  // Build send lists of requests
  for (List<uint>::iterator it(vertices_to_request_.begin());
      it != vertices_to_request_.end(); ++it)
  {
    dolfin_assert(old2new_vertices[*it] >= 0);
    uint global_index = mesh.distdata()[0].get_global(old2new_vertices[*it]);

    // vertex belongs to other process: request has to be distributed by owner
    if (mesh.distdata()[0].is_ghost(old2new_vertices[*it]))
    {
      uint owner = mesh.distdata()[0].get_owner(old2new_vertices[*it]);
      send_list_requests[owner].push_back(global_index);
    }
    // vertex belongs to this process: put request into map
    else
    {
      requested_vertices[global_index] = rank;
    }
  }

  // pairwise communication to exchange requests
  MPI_Status status;
  int recv_size;
  for (uint i(1); i < pe_size; ++i)
  {
    uint src = (rank - i + pe_size) % pe_size;
    uint dest = (rank + i) % pe_size;

    MPI_Sendrecv(&send_list_requests[dest][0], send_list_requests[dest].size(),
                 MPI_UNSIGNED, dest, 0, recv_buff_requests,
                 max_num_requested_vertices, MPI_UNSIGNED, src, 0,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_size);

    // process received requests and puts them into the map
    std::map<uint, uint>::iterator m_it;
    for (uint i(0); i < recv_size; ++i)
    {
      // search for this index in the map
      uint requested_vertex = recv_buff_requests[i];
      m_it = requested_vertices.find(requested_vertex);

      // Another process also requested this vertex: lower rank wins
      if (m_it != requested_vertices.end()) m_it->second = std::min(
          src, m_it->second);
      // no conflict: simply insert into map
      else requested_vertices[requested_vertex] = src;
    }
  }

  // clear buffers
  for (uint i(0); i < pe_size; ++i)
  {
    send_list_requests[i].clear();
  }

  // New partitions according to requests. Initialized with current partitions
  partitions = new MeshFunction<uint>(mesh, mesh.topology().dim());
  *partitions = rank;
  uint num_send_cells(0);
  MeshFunction<bool> requested_cells(mesh, mesh.topology().dim());
  requested_cells = false;

  // Build send list of vertices with requesting processes
  for (std::map<uint, uint>::iterator m_it(requested_vertices.begin());
      m_it != requested_vertices.end(); ++m_it)
  {
    uint local_index = mesh.distdata()[0].get_local(m_it->first);
    uint pe = m_it->second;

    // set of processes that share this vertex
    _set<uint> shared_adj = mesh.distdata()[0].get_shared_adj(local_index);
    for (_set<uint>::iterator s_it(shared_adj.begin());
    s_it != shared_adj.end(); ++s_it )
    {
      if ( *s_it == rank ) continue;

      send_list_requests[*s_it].push_back(m_it->first);
      send_list_requests[*s_it].push_back(m_it->second);
    }

    // select lowest process index for all cells around requested vertex
    Vertex v(mesh, local_index);
    uint target_proc = m_it->second;

    for (CellIterator c_it(v); !c_it.end(); ++c_it)
    {
      if (partitions->get(*c_it) == rank && target_proc != rank) ++num_send_cells;
      if (requested_cells(*c_it)) partitions->set(
          *c_it, std::min(target_proc, partitions->get(*c_it)));
      else partitions->set(*c_it, target_proc);
      requested_cells.set(*c_it, true);
    }
  }

  // pairwise communication to exchange requests
  for (uint i(1); i < pe_size; ++i)
  {
    uint src = (rank - i + pe_size) % pe_size;
    uint dest = (rank + i) % pe_size;

    MPI_Sendrecv(&send_list_requests[dest][0], send_list_requests[dest].size(),
                 MPI_UNSIGNED, dest, 0, recv_buff_requests,
                 2 * max_num_requested_vertices, MPI_UNSIGNED, src, 0,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_size);

    // process received requests and marks cells accordingly
    for (uint i(0); i < recv_size; i += 2)
    {
      uint local_index = mesh.distdata()[0].get_local(recv_buff_requests[i]);
      Vertex v(mesh, local_index);

      // select lowest process index for all cells around requested vertex
      uint target_proc = recv_buff_requests[i + 1];
      for (CellIterator c_it(v); !c_it.end(); ++c_it)
        target_proc = std::min(target_proc, partitions->get(*c_it));

      // set partitions
      for (CellIterator c_it(v); !c_it.end(); ++c_it)
      {
        if (partitions->get(*c_it) == rank && target_proc != rank) ++num_send_cells;
        partitions->set(*c_it, target_proc);
      }
    }
  }

  // Check that at least one cell resides here
  if (num_send_cells >= mesh.num_cells())
  {
    partitions->set(0, rank);
    --num_send_cells;
  }

  num_migrated_cells_ = num_send_cells;

  // clear buffers
  for (uint i(0); i < pe_size; ++i)
  {
    send_list_requests[i].clear();
  }
  delete[] send_list_requests;
  delete[] recv_buff_requests;
#endif
}
//-----------------------------------------------------------------------------
bool CoarseningManager::migrate(uint num_cells_coarsened)
{
#ifdef HAVE_MPI
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();

  if (rank == 0) message("Starting migration.");

  // Cleanup Coarsening List
  removeErasedCellsFromCoarseningList();

  if (pe_size == 1)
  {
    dolfin_assert(cells_to_request_.size() == 0);
    dolfin_assert(vertices_to_request_.size() == 0);
    return (num_cells_coarsened > 0) && vertices_to_request_.size() == 0;
  }

  // determine global numbers of coarsened cells and remaining cells
  uint global_nums[4];
  uint local_nums[4];
  local_nums[0] = num_cells_coarsened;
  local_nums[1] = cells_to_coarsen_.size();
  local_nums[2] = num_migrated_cells_;
  local_nums[3] = vertices_to_request_.size();
  MPI_Allreduce(local_nums, global_nums, 4, MPI_UNSIGNED, MPI_MAX,
                MPI::DOLFIN_COMM);

  // no cells remaining: we're done!
  if (global_nums[1] == 0) return false;
  // No migration happened the last time and nothing coarsened: we're done!
  else if (global_nums[0] == 0 && global_nums[2] == 0) return false;

  uint max_num_requested_vertices = global_nums[3];

  // no vertices requested: no migration necessary, but obviously still work left
  if (max_num_requested_vertices == 0) return true;

  // Export DMesh to Mesh
  Array<int> old2new_cells(orig_num_cells_);
  Array<int> old2new_vertices(orig_num_vertices_);
  Mesh omesh;
  dmesh_->expKeepNumbering(omesh, &old2new_cells, &old2new_vertices);

  // Compute cell-vertex connectivity
  omesh.init(0, omesh.topology().dim());

  // Pointer to MeshFunction with new partitions
  MeshFunction<uint> *partitions = 0;

  // Lists of MeshFunctions for exchange
  Array<std::pair<MeshFunction<uint> *, MeshFunction<uint> *> > cell_functions;
  Array<std::pair<MeshFunction<real> *, MeshFunction<real> *> > vertex_functions;

  buildMFArrays(omesh, old2new_cells, old2new_vertices, cell_functions,
                vertex_functions);

  // Exchange requests
  exchangeRequests(omesh, old2new_cells, old2new_vertices,
                   max_num_requested_vertices, partitions);

  // distribute partitioning and MeshFunctions
  omesh.distribute(*partitions, cell_functions, vertex_functions);
  omesh.renumber();

  // Re-initialize
  initCommon(omesh, cell_functions[0].second);

  // Update independent set
  updateIndependentSet(omesh, *(vertex_functions.front().second));

  // Cleanup MeshFunction Arrays
  cleanupMFArrays(cell_functions, vertex_functions);

  // free memory
  delete partitions;
#endif
  return true;
}
//-----------------------------------------------------------------------------

