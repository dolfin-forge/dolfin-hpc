// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-03-25
// Last changed: 2013-06-01

#include <dolfin/mesh/CoarseningManager.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/LoadBalancer.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/mesh/MeshFunctionConverter.h>
#include <dolfin/main/MPI.h>

#ifdef ____USE_D_MESH____
#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/DCell.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;
//-----------------------------------------------------------------------------
CoarseningManager::CoarseningManager() 
#ifdef ____USE_D_MESH____
: _dmesh(0)
#endif
{
  // do nothing
}
//-----------------------------------------------------------------------------
CoarseningManager::CoarseningManager(Mesh& mesh, MeshFunction<bool>& cell_marker, 
                                     bool coarsen_boundary)
#ifdef ____USE_D_MESH____
: _dmesh(0)
#endif
{
  init(mesh, cell_marker, coarsen_boundary);
}
//-----------------------------------------------------------------------------
CoarseningManager::~CoarseningManager()
{
#ifdef ____USE_D_MESH____
  if ( _dmesh )
    delete _dmesh;
#endif
}
//-----------------------------------------------------------------------------
void CoarseningManager::init(Mesh& mesh, MeshFunction<bool>& cell_marker, 
                             bool coarsen_boundary)
{
  dolfin_assert( &mesh == &(cell_marker.mesh()) );

  _global_max_coarsened_cells = 0;
  _global_max_remaining_cells = -1;
  _migrations = 0;
  _load_balances = 0;

  MeshFunction<uint> attempt_count(mesh, mesh.topology().dim());
  attempt_count = 0;
  initCommon(mesh, cell_marker, &attempt_count);

  //dolfin_assert( _int_bnd_vertices.size() == cell_marker.mesh().numVertices() );
  //dolfin_assert( _int_bnd_cells.size() == cell_marker.mesh().numCells() );
  dolfin_assert( _cells_to_coarsen.size() >= 0 );

  // find independent set
  findIndependentSet(mesh, coarsen_boundary);

  dolfin_assert( _forbidden_vertices.size() == cell_marker.mesh().numVertices() );
}
//-----------------------------------------------------------------------------
template<typename T>
void CoarseningManager::initCommon(Mesh& mesh, MeshFunction<T>& cell_marker,
                                   MeshFunction<uint> * attempt_count)
{
  dolfin_assert( &mesh == &(cell_marker.mesh()) );

#ifdef ____USE_D_MESH____
  /// Delete old dmesh and import new mesh into dmesh
  if ( _dmesh )
    delete _dmesh;
  _dmesh = new DMesh;
  _dmesh->imp(mesh);
  
  _orig_num_cells = mesh.numCells();
  _orig_num_vertices = mesh.numVertices();
#else
  // extract boundary information
  findInteriorBoundaries(mesh);
#endif 

  // build list of cells to coarsen
  findCellsToCoarsen(cell_marker, attempt_count);

  // initialize mapping
#ifndef ____USE_D_MESH____
  _vertex_map.init(mesh.numVertices());
  _cell_map.init(cell_marker.size());
#endif

  // reset request list
  _cells_to_request.clear();
  _vertices_to_request.clear();
}
//-----------------------------------------------------------------------------
#ifndef ____USE_D_MESH____
void CoarseningManager::findInteriorBoundaries(Mesh& mesh)
{
  // set to false on whole domain
  _int_bnd_vertices.resize(mesh.numVertices());
  _int_bnd_vertices = false;
  //_int_bnd_cells.resize(mesh.numCells());
  //_int_bnd_cells = false;

  // no process boundaries if only one process
  if ( MPI::numProcesses() == 1 )
    return;

  // get boundary mesh
  BoundaryMesh boundary;
  boundary.init_interior(mesh);

  // no boundary vertices - nothing to do here
  if ( boundary.numVertices() == 0 )
    return;

  MeshFunction<uint> *bnd_vertex_map = boundary.data().meshFunction("vertex map");
  dolfin_assert(bnd_vertex_map);

  // fill Mesh functions
  for ( VertexIterator v_it(boundary) ; !v_it.end() ; ++v_it )
  {
    Vertex v(mesh, bnd_vertex_map->get(v_it->index()));
    _int_bnd_vertices.at(v.index()) = true;

    // all connected cells are also on the boundary
    //for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
    //  _int_bnd_cells.at(c_it->index()) = true;
  }
}
#endif
//-----------------------------------------------------------------------------
/*void CoarseningManager::findDomainBoundaries(Mesh& mesh)
{
  // set to false on whole domain
  _bnd_vertices.resize(mesh.numVertices());
  _bnd_vertices = false;
  _bnd_cells.resize(mesh.numCells());
  _bnd_cells = false;

  // get boundary mesh
  BoundaryMesh boundary;
  boundary.init(mesh);

  // no boundary vertices - nothing to do here
  if ( boundary.numVertices() == 0 )
    return;

  MeshFunction<uint> *bnd_vertex_map = boundary.data().meshFunction("vertex map");
  dolfin_assert(bnd_vertex_map);

  // fill Mesh functions
  for ( VertexIterator v_it(boundary) ; !v_it.end() ; ++v_it )
  {
    Vertex v(mesh, bnd_vertex_map->get(v_it->index()));
    _bnd_vertices.at(v.index()) = true;

    // all connected cells are also on the boundary
    for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
      _bnd_cells.at(c_it->index()) = true;
  }
}*/
//-----------------------------------------------------------------------------
void CoarseningManager::findIndependentSet(Mesh& mesh, bool coarsen_boundary)
{
  // set to false on whole domain
  _forbidden_vertices.resize(mesh.numVertices());
  _forbidden_vertices = false;

  // if boundary coarsening is forbidden: put boundary vertices into set first
  if ( !coarsen_boundary )
  {
    BoundaryMesh boundary;
    boundary.init(mesh);

    if ( boundary.numVertices() > 0 )
    {
      MeshFunction<uint> *bnd_vertex_map = boundary.data().meshFunction("vertex map");

      // add boundary vertices to set
      for ( VertexIterator v_it(boundary) ; !v_it.end() ; ++v_it )
        _forbidden_vertices.at(bnd_vertex_map->get(v_it->index())) = true;
    }
  }

  // iterate over remaining vertices
  for ( VertexIterator v_it(mesh) ; !v_it.end() ; ++v_it )
  {
    if ( !_forbidden_vertices[v_it->index()] )
      _forbidden_vertices.at(v_it->index()) = isIndependentVertex(*v_it);
  }
}
//-----------------------------------------------------------------------------
bool CoarseningManager::isIndependentVertex(Vertex& v)
{
  // Check all vertices connected to that vertex
  for ( VertexIterator v_it(v) ; !v_it.end() ; ++v_it )
  {
    if ( v_it->index() == v.index() )
      continue;

    if ( _forbidden_vertices.at(v_it->index()) )
      return false;
  }

  // no neighboring point in set: independent vertex
  return true;
}
//-----------------------------------------------------------------------------
#ifdef ____USE_D_MESH____
//-----------------------------------------------------------------------------
template<typename T>
void CoarseningManager::findCellsToCoarsen(MeshFunction<T>& cell_marker,
                                           MeshFunction<uint> * attempt_count)
{
  _cells_to_coarsen.clear();

  for ( std::list<DCell *>::iterator c_it(_dmesh->cells.begin()) ;
        c_it != _dmesh->cells.end() ; ++c_it )
  {
    DCell * dc = *c_it;
    if ( cell_marker.get(dc->id) > 0 )
      _cells_to_coarsen.push_back(std::make_pair(dc, attempt_count->get(dc->id)));
  }
}
//-----------------------------------------------------------------------------
void CoarseningManager::removeErasedCellsFromCoarseningList()
{
  // remove deleted entities from coarsening list
  for ( List< std::pair<DCell *,uint> >::iterator it(_cells_to_coarsen.begin()) ; 
        it != _cells_to_coarsen.end() ; )
  {
    DCell * dc = it->first;
    if ( dc->deleted ) 
      it = _cells_to_coarsen.erase(it);
    else
      ++it;
  }
}
//-----------------------------------------------------------------------------
void CoarseningManager::buildMFArrays(Mesh& mesh, Array<int>& old2new_cells,
                                      Array<int>& old2new_vertices,
                                      Array< std::pair< MeshFunction<uint>*, 
                                      MeshFunction<uint>* > >& cell_functions,
                                      Array< std::pair< MeshFunction<double>*,
                                      MeshFunction<double>*> >& vertex_functions)
{
  // Prepare forbidden vertices for exchange
  MeshFunction<double> * forbidden_vertices = new MeshFunction<double>(mesh, 0);
  MeshFunction<double> * forbidden_vertices_new = new MeshFunction<double>();
  for ( uint old_idx(0) ; old_idx < _forbidden_vertices.size();
        ++old_idx )
  {
    int new_idx = old2new_vertices[old_idx];
    if ( new_idx < 0 )
      continue;
    forbidden_vertices->set(new_idx, double(isForbiddenVertex(old_idx)) );
  }
  vertex_functions.push_back( 
    std::make_pair(forbidden_vertices, forbidden_vertices_new) );

  // Prepare cell_marker and attempt count for exchange
  MeshFunction<uint> * cell_marker = new MeshFunction<uint>(mesh, mesh.topology().dim());
  MeshFunction<uint> * cell_marker_new = new MeshFunction<uint>();
  MeshFunction<uint> * attempts = new MeshFunction<uint>(mesh, mesh.topology().dim());
  MeshFunction<uint> * attempts_new = new MeshFunction<uint>();
  *cell_marker = 0u;
  *attempts = 0u;
  for ( List< std::pair<DCell *, uint> >::iterator it(_cells_to_coarsen.begin()) ; 
        it != _cells_to_coarsen.end() ; ++it )
  {
    DCell * dc = it->first;
    int new_idx = old2new_cells[dc->id];
    if ( new_idx < 0 ) 
      continue;
    cell_marker->set(new_idx, 1u);
    attempts->set(new_idx,it->second);
  }
  cell_functions.push_back( std::make_pair(cell_marker, cell_marker_new) );
  cell_functions.push_back( std::make_pair(attempts, attempts_new) );
}
//-----------------------------------------------------------------------------
void CoarseningManager::cleanupMFArrays(Array< std::pair< MeshFunction<uint>*, 
                                        MeshFunction<uint>* > >& cell_functions,
                                        Array< std::pair< MeshFunction<double>*,
                                        MeshFunction<double>*> >& vertex_functions)
{
  for ( Array< std::pair<MeshFunction<uint>*,MeshFunction<uint>*> >::iterator
        it(cell_functions.begin()) ; it != cell_functions.end() ; ++it )
  {
    delete it->first;
    delete it->second;
  }

  for ( Array< std::pair<MeshFunction<double>*,MeshFunction<double>*> >::iterator
        it(vertex_functions.begin()) ; it != vertex_functions.end() ; ++it )
  {
    delete it->first;
    delete it->second;
  }
}
//-----------------------------------------------------------------------------
void CoarseningManager::updateIndependentSet(Mesh& mesh, MeshFunction<double>& 
                                             forbidden_vertices_new)
{
  _forbidden_vertices.resize(mesh.numVertices());
  _forbidden_vertices = false;
  for ( VertexIterator v_it(mesh) ; !v_it.end() ; ++v_it )
    if ( forbidden_vertices_new.get(v_it->index()) > 0.5 )
      _forbidden_vertices[v_it->index()] = true;
}
//-----------------------------------------------------------------------------
void CoarseningManager::exchangeRequests(Mesh& mesh, Array<int>& old2new_cells,
                                         Array<int>& old2new_vertices,
                                         uint max_num_requested_vertices,
                                         MeshFunction<uint> *& partitions)
{
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();

  // List of vertices to request from owner
  Array<uint> *send_list_requests = new Array<uint>[pe_size];

  // Received requests
  uint * recv_buff_requests = new uint[2 * max_num_requested_vertices];

  // Map of requested vertices (that this process owns) and requesting processes
  std::map<uint,uint> requested_vertices;

  // Build send lists of requests
  for ( List<uint>::iterator it(_vertices_to_request.begin()) ; 
        it != _vertices_to_request.end() ; ++it )
  {
    dolfin_assert( isInteriorBoundaryVertex(*it) );
    dolfin_assert( old2new_vertices[*it] >= 0 );
    uint global_index = mesh.distdata().get_global(old2new_vertices[*it], 0);

    // vertex belongs to other process: request has to be distributed by owner
    if ( mesh.distdata().is_ghost(old2new_vertices[*it], 0) )
    {
      uint owner = mesh.distdata().get_owner(old2new_vertices[*it], 0);
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
  for ( uint i(1) ; i < pe_size ; ++i )
  {
    uint src = (rank - i + pe_size) % pe_size;
    uint dest = (rank + i) % pe_size;

    MPI_Sendrecv( &send_list_requests[dest][0], send_list_requests[dest].size(),
      MPI_UNSIGNED, dest, 0, 
      recv_buff_requests, max_num_requested_vertices, MPI_UNSIGNED, src, 0,
      MPI::DOLFIN_COMM, &status );
    MPI_Get_count( &status, MPI_UNSIGNED, &recv_size );

    // process received requests and puts them into the map
    std::map<uint,uint>::iterator m_it;
    for ( uint i(0) ; i < recv_size ; ++i )
    {
      /*
      uint requested_vertex = recv_buff_requests[i];
      uint local_index = mesh.distdata().get_local(requested_vertex, 0);
      uint lowest_sharing_rank = 
            *( mesh.distdata().get_shared_adj(local_index, 0).begin() );
      requested_vertices[requested_vertex] = lowest_sharing_rank;
      /*/
      // search for this index in the map
      uint requested_vertex = recv_buff_requests[i];
      m_it = requested_vertices.find(requested_vertex);

      // Another process also requested this vertex: lower rank wins
      if ( m_it != requested_vertices.end() )
        m_it->second = std::min(src, m_it->second);
      // no conflict: simply insert into map
      else
        requested_vertices[requested_vertex] = src;
    }
  }

  // clear buffers
  for ( uint i(0) ; i < pe_size ; ++i )
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
  for ( std::map<uint,uint>::iterator m_it(requested_vertices.begin()) ; 
        m_it != requested_vertices.end() ; ++m_it )
  {
    uint local_index = mesh.distdata().get_local(m_it->first, 0);
    uint pe = m_it->second;

    // set of processes that share this vertex
    _set<uint> shared_adj = mesh.distdata().get_shared_adj(local_index, 0);
    for ( _set<uint>::iterator s_it(shared_adj.begin()) ;
          s_it != shared_adj.end() ; ++s_it )
    {
      if ( *s_it == rank ) continue;

      send_list_requests[*s_it].push_back(m_it->first);
      send_list_requests[*s_it].push_back(m_it->second);
    }

    // select lowest process index for all cells around requested vertex
    Vertex v(mesh, local_index);
    uint target_proc = m_it->second;

    /*for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
      target_proc = std::min(target_proc, partitions->get(*c_it));

    // set partitions
    for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
    {
      if ( partitions->get(*c_it) == rank && target_proc != rank ) 
        ++num_send_cells;
      partitions->set(*c_it, target_proc);
    }/*/
    for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
    {
      if ( partitions->get(*c_it) == rank && target_proc != rank )
        ++num_send_cells;
      if ( requested_cells(*c_it) )
        partitions->set(*c_it, std::min(target_proc, partitions->get(*c_it)));
      else
        partitions->set(*c_it, target_proc);
      requested_cells.set(*c_it, true);
    }
  }

  // pairwise communication to exchange requests
  for ( uint i(1) ; i < pe_size ; ++i )
  {
    uint src = (rank - i + pe_size) % pe_size;
    uint dest = (rank + i) % pe_size;

    MPI_Sendrecv( &send_list_requests[dest][0], send_list_requests[dest].size(),
      MPI_UNSIGNED, dest, 0, 
      recv_buff_requests, 2 * max_num_requested_vertices, MPI_UNSIGNED, src, 0,
      MPI::DOLFIN_COMM, &status );
    MPI_Get_count( &status, MPI_UNSIGNED, &recv_size );

    // process received requests and marks cells accordingly
    for ( uint i(0) ; i < recv_size ; i += 2 )
    {
      uint local_index = mesh.distdata().get_local(recv_buff_requests[i], 0);
      Vertex v(mesh, local_index);

      // select lowest process index for all cells around requested vertex
      uint target_proc = recv_buff_requests[i+1];
      for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
        target_proc = std::min(target_proc, partitions->get(*c_it));

      // set partitions
      for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
      {
        if ( partitions->get(*c_it) == rank && target_proc != rank ) 
          ++num_send_cells;
        partitions->set(*c_it, target_proc);
      }
    }
  }

  // clear buffers
  for ( uint i(0) ; i < pe_size ; ++i )
  {
    send_list_requests[i].clear();
  }
  delete[] send_list_requests;
  delete[] recv_buff_requests;
}
//-----------------------------------------------------------------------------
bool CoarseningManager::migrate(uint num_cells_coarsened)
{
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint num_migrations_before_loadbalancing = 15;

  if (rank == 0)
    cout << "Starting migration." << endl;

  // Cleanup Coarsening List
  removeErasedCellsFromCoarseningList();

  if ( pe_size == 1 )
  {
    dolfin_assert(_cells_to_request.size() == 0);
    dolfin_assert(_vertices_to_request.size() == 0);
    return (num_cells_coarsened > 0);
  }

  // determine global numbers of coarsened cells and remaining cells
  uint global_nums[3];
  uint local_nums[3];
  local_nums[0] = num_cells_coarsened;
  local_nums[1] = _cells_to_coarsen.size();
  local_nums[2] = _vertices_to_request.size();
  MPI_Allreduce(local_nums, global_nums, 3, 
                MPI_UNSIGNED, MPI_MAX, MPI::DOLFIN_COMM);

  if ( rank == 0 )
    cout << "Max " << global_nums[0] << " coarsened, " << 
            global_nums[1] << " remaining, " <<
            global_nums[2] << " requested." << endl;

  // No cells left for coarsening: we're done!
  if ( global_nums[1] == 0 )
    return false;
  // Nothing happened since last time: we consider ourselves done!
  else if ( global_nums[0] == _global_max_coarsened_cells &&
            global_nums[1] == _global_max_remaining_cells )
    return false;

  _global_max_coarsened_cells = global_nums[0];
  _global_max_remaining_cells = global_nums[1];
  uint max_num_requested_vertices = global_nums[2];

  // Export DMesh to Mesh
  Array<int> old2new_cells(_orig_num_cells);
  Array<int> old2new_vertices(_orig_num_vertices);
  Mesh omesh;
  _dmesh->expKeepNumbering(omesh, &old2new_cells, &old2new_vertices);

  // Compute cell-vertex connectivity
  omesh.init(0,omesh.topology().dim());

  // Pointer to MeshFunction with new partitions
  MeshFunction<uint> *partitions = 0;

  // Lists of MeshFunctions for exchange
  Array< std::pair< MeshFunction<uint> * , MeshFunction<uint> * > > 
                                                              cell_functions;
  Array< std::pair< MeshFunction<double> * , MeshFunction<double> * > > 
                                                              vertex_functions;
  
  buildMFArrays(omesh, old2new_cells, old2new_vertices, 
                cell_functions, vertex_functions);

  // Exchange requests
  if ( _migrations < num_migrations_before_loadbalancing )
  {
    exchangeRequests(omesh, old2new_cells, old2new_vertices, 
                     max_num_requested_vertices, partitions);
  }
  else
  {
    // obtain new partitions from loadbalancer
    //partitions = new MeshFunction<uint>();
    //omesh.partition(*partitions);

    // obtain new partitions from loadbalancer
    dolfin_set("Load balancer redistribute", false);
    MeshFunction<bool> cell_marker_b;
    MeshFunctionConverter::cast(*(cell_functions.front().first), cell_marker_b);
    //LoadBalancer::balance(omesh, cell_marker_b);
    LoadBalancer::balance(omesh, cell_marker_b, LoadBalancer::EdgeCollapse);
    partitions = omesh.data().meshFunction("partitions");
    dolfin_assert(partitions);
  }

  // distribute partitioning and MeshFunctions
  omesh.distribute( *partitions, cell_functions, vertex_functions );
  omesh.renumber();

  dolfin_assert(omesh.numCells() > 0);

  // Re-initialize
  initCommon(omesh, *(cell_functions[0].second), cell_functions[1].second);

  // Update independent set
  updateIndependentSet(omesh, *(vertex_functions.front().second));
  
  // Cleanup MeshFunction Arrays
  cleanupMFArrays(cell_functions, vertex_functions);

  // increase counter
  if ( _migrations < num_migrations_before_loadbalancing )
  {
    ++_migrations;
    delete partitions;
  }
  else
  {
    _migrations = 0;
    ++_load_balances;
  }

  return true;
}
//-----------------------------------------------------------------------------
#else // ____USE_D_MESH____
//-----------------------------------------------------------------------------
template<typename T>
void CoarseningManager::findCellsToCoarsen(MeshFunction<T>& cell_marker)
{
  _cells_to_coarsen.clear();

  for ( CellIterator c_it(cell_marker.mesh()) ; !c_it.end() ; ++c_it )
  {
    if ( cell_marker.get(c_it->index()) )
      _cells_to_coarsen.push_back(c_it->index());
  }
}
//-----------------------------------------------------------------------------
bool CoarseningManager::migrate(Mesh& mesh, bool repeat)
{
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();

  if ( pe_size == 1 )
  {
    dolfin_assert(_cells_to_request.size() == 0);
    dolfin_assert(_vertices_to_request.size() == 0);
    return repeat;
  }

  // exchange maximum number of requested vertices
  int local_status[3], global_status[3];
  local_status[0] = _vertices_to_request.size(); // local number of vertices
  local_status[1] = repeat;                      // local termination?
  local_status[2] = _migrated_cells;             // number of cells migrated
                                                 // in last migration
  MPI_Allreduce(&local_status, &global_status, 3, MPI_INT, 
                MPI_MAX, MPI::DOLFIN_COMM);

  uint max_num_requested_vertices = global_status[0]; // max number of vertices
  repeat = global_status[1];                          // global termination?
  uint max_num_migrated_cells = global_status[2];     // max num migrated cells

  // migration nowhere necessary
  if ( max_num_requested_vertices == 0 )
  {
    _migrated_cells = 0;
    return repeat;
  }

  // no cells migrated before and nothing happened in the coarsening: we're done
  if ( max_num_migrated_cells == 0 && !repeat )
  {
    _migrated_cells = 0;
    return false;
  }

  // List of vertices to request from owner
  Array<uint> *send_list_requests = new Array<uint>[pe_size];

  // Received requests
  uint * recv_buff_requests = new uint[2 * max_num_requested_vertices];

  // Map of requested vertices (that this process owns) and requesting processes
  std::map<uint,uint> requested_vertices;

  // Build send lists of requests
  for ( List<uint>::iterator it(_vertices_to_request.begin()) ; 
        it != _vertices_to_request.end() ; ++it )
  {
    dolfin_assert( isInteriorBoundaryVertex(*it) );
    int coarse_index = _vertex_map.getCoarseFromFine(*it);
    dolfin_assert( coarse_index >= 0 );
    uint global_index = mesh.distdata().get_global(coarse_index, 0);

    // vertex belongs to other process: request has to be distributed by owner
    if ( mesh.distdata().is_ghost(coarse_index, 0) )
    {
      uint owner = mesh.distdata().get_owner(coarse_index, 0);
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
  for ( uint i(1) ; i < pe_size ; ++i )
  {
    uint src = (rank - i + pe_size) % pe_size;
    uint dest = (rank + i) % pe_size;

    MPI_Sendrecv( &send_list_requests[dest][0], send_list_requests[dest].size(),
      MPI_UNSIGNED, dest, 0, 
      recv_buff_requests, max_num_requested_vertices, MPI_UNSIGNED, src, 0,
      MPI::DOLFIN_COMM, &status );
    MPI_Get_count( &status, MPI_UNSIGNED, &recv_size );

    // process received requests and puts them into the map
    std::map<uint,uint>::iterator m_it;
    for ( uint i(0) ; i < recv_size ; ++i )
    {
      // search for this index in the map
      uint requested_vertex = recv_buff_requests[i];
      m_it = requested_vertices.find(requested_vertex);

      // Another process also requested this vertex: lower rank wins
      if ( m_it != requested_vertices.end() )
        m_it->second = std::min(src, m_it->second);
      // no conflict: simply insert into map
      else
        requested_vertices[requested_vertex] = src;
    }
  }

  // clear buffers
  for ( uint i(0) ; i < pe_size ; ++i )
  {
    send_list_requests[i].clear();
  }

  // New partitions according to requests. Initialized with current partitions
  MeshFunction<uint> partitions(mesh, mesh.topology().dim());
  partitions = rank;
  uint num_send_cells(0);

  // Build send list of vertices with requesting processes
  for ( std::map<uint,uint>::iterator m_it(requested_vertices.begin()) ; 
        m_it != requested_vertices.end() ; ++m_it )
  {
    uint local_index = mesh.distdata().get_local(m_it->first, 0);
    uint pe = m_it->second;

    // set of processes that share this vertex
    _set<uint> shared_adj = mesh.distdata().get_shared_adj(local_index, 0);
    for ( _set<uint>::iterator s_it(shared_adj.begin()) ;
          s_it != shared_adj.end() ; ++s_it )
    {
      if ( *s_it == rank ) continue;

      send_list_requests[*s_it].push_back(m_it->first);
      send_list_requests[*s_it].push_back(m_it->second);
    }

    // select lowest process index for all cells around requested vertex
    Vertex v(mesh, local_index);
    uint target_proc = m_it->second;
    for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
      target_proc = std::min(target_proc, partitions(*c_it));

    // set partitions
    for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
    {
      if ( partitions(*c_it) == rank && target_proc != rank ) 
        ++num_send_cells;
      partitions(*c_it) = target_proc;
    }
  }

  // pairwise communication to exchange requests
  for ( uint i(1) ; i < pe_size ; ++i )
  {
    uint src = (rank - i + pe_size) % pe_size;
    uint dest = (rank + i) % pe_size;

    MPI_Sendrecv( &send_list_requests[dest][0], send_list_requests[dest].size(),
      MPI_UNSIGNED, dest, 0, 
      recv_buff_requests, 2 * max_num_requested_vertices, MPI_UNSIGNED, src, 0,
      MPI::DOLFIN_COMM, &status );
    MPI_Get_count( &status, MPI_UNSIGNED, &recv_size );

    // process received requests and marks cells accordingly
    for ( uint i(0) ; i < recv_size ; i += 2 )
    {
      uint local_index = mesh.distdata().get_local(recv_buff_requests[i], 0);
      Vertex v(mesh, local_index);

      // select lowest process index for all cells around requested vertex
      uint target_proc = recv_buff_requests[i+1];
      for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
        target_proc = std::min(target_proc, partitions(*c_it));

      // set partitions
      for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
      {
        if ( partitions(*c_it) == rank && target_proc != rank ) 
          ++num_send_cells;
        partitions(*c_it) = target_proc;
      }
    }
  }

  // clear buffers
  for ( uint i(0) ; i < pe_size ; ++i )
  {
    send_list_requests[i].clear();
  }
  delete[] send_list_requests;
  delete[] recv_buff_requests;

  // Lists of MeshFunctions for exchange
  Array< std::pair< MeshFunction<uint> * , MeshFunction<uint> * > > 
                                                              cell_functions;
  Array< std::pair< MeshFunction<double> * , MeshFunction<double> * > > 
                                                              vertex_functions;
  
  // Prepare forbidden vertices for exchange
  MeshFunction<double> forbidden_vertices(mesh, 0);
  MeshFunction<double> forbidden_vertices_new;
  for ( VertexIterator v_it(mesh) ; !v_it.end() ; ++v_it )
    forbidden_vertices.set(*v_it, double( 
      isForbiddenVertex( _vertex_map.getFineFromCoarse(v_it->index()) ) ) );
  vertex_functions.push_back( 
    std::make_pair(&forbidden_vertices, &forbidden_vertices_new) );

  // Prepare cell_marker for coarsening for exchange
  MeshFunction<uint> cell_marker(mesh, mesh.topology().dim());
  MeshFunction<uint> cell_marker_new;
  cell_marker = false;
  for ( List<uint>::iterator it(_cells_to_request.begin()) ; 
        it != _cells_to_request.end() ; ++it )
  {
    int coarse_index = _cell_map.getCoarseFromFine(*it);
    if ( coarse_index < 0 ) 
      continue;
    cell_marker.set(coarse_index, 1);
  }
  for ( List<uint>::iterator it(_cells_to_coarsen.begin()) ; 
        it != _cells_to_coarsen.end() ; ++it )
  {
    int coarse_index = _cell_map.getCoarseFromFine(*it);
    if ( coarse_index < 0 ) 
      continue;
    cell_marker.set(coarse_index, 1);
  }
  cell_functions.push_back( std::make_pair(&cell_marker, &cell_marker_new) );

  _migrated_cells = num_send_cells;

  // distribute partitioning and MeshFunctions
  mesh.distribute( partitions, cell_functions, vertex_functions );

  // Re-initialize
  initCommon(mesh, cell_marker_new);

  // Update independent set
  //MeshFunctionConverter::cast(forbidden_vertices_new, _forbidden_vertices);
  _forbidden_vertices.clear();
  _forbidden_vertices.resize(mesh.numVertices(), false);
  for ( VertexIterator v_it(mesh) ; !v_it.end() ; ++v_it )
    if ( forbidden_vertices_new.get(v_it->index()) > 0.5 )
      _forbidden_vertices[v_it->index()] = true;

  return true;
}
//-----------------------------------------------------------------------------
void CoarseningManager::updateDistdata(Mesh& old_mesh, Mesh& new_mesh)
{
  new_mesh.distdata().clear();

  // vertices
  for ( VertexIterator v_it(new_mesh) ; !v_it.end() ; ++v_it )
  {
    // set global index
    uint old_index = _vertex_map.getCoarseFromNewCoarse(v_it->index());
    uint glb_index = old_mesh.distdata().get_global(old_index, 0);
    new_mesh.distdata().set_map(v_it->index(), glb_index, 0);

    // set ghost/shared state
    if (old_mesh.distdata().is_ghost(old_index, 0))
    {
      new_mesh.distdata().set_ghost(v_it->index(), 0);
      new_mesh.distdata().set_ghost_owner(v_it->index(), 
                  old_mesh.distdata().get_owner(old_index, 0), 0);
    }
    else if (old_mesh.distdata().is_shared(old_index, 0))
    {
      new_mesh.distdata().set_shared(v_it->index(), 0);
      new_mesh.distdata().get_shared_adj(v_it->index(), 0) =
                  old_mesh.distdata().get_shared_adj(old_index, 0);
    }
    else
      continue;
  }
}
//-----------------------------------------------------------------------------
#endif // ____USE_D_MESH____
//-----------------------------------------------------------------------------
