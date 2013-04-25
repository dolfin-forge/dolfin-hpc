// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-03-25
// Last changed: 2013-04-03

#include <dolfin/mesh/CoarseningManager.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
//#include <dolfin/mesh/MeshFunctionConverter.h>
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
CoarseningManager::CoarseningManager(MeshFunction<bool>& cell_marker, 
                                     bool coarsen_boundary)
#ifdef ____USE_D_MESH____
: _dmesh(0)
#endif
{
  init(cell_marker, coarsen_boundary);
}
//-----------------------------------------------------------------------------
CoarseningManager::~CoarseningManager()
{
  if ( _dmesh )
    delete _dmesh;
}
//-----------------------------------------------------------------------------
void CoarseningManager::init(MeshFunction<bool>& cell_marker, 
                             bool coarsen_boundary)
{
  // migrated cells are initially set to a value larger zero to match the
  // break conditions
  _migrated_cells = 1;

  initCommon(cell_marker);

  dolfin_assert( _int_bnd_vertices.size() == cell_marker.mesh().numVertices() );
  dolfin_assert( _int_bnd_cells.size() == cell_marker.mesh().numCells() );
  dolfin_assert( _bnd_vertices.size() == cell_marker.mesh().numVertices() );
  dolfin_assert( _bnd_cells.size() == cell_marker.mesh().numCells() );
  dolfin_assert( _cells_to_coarsen.size() >= 0 );

  // find independent set
  findIndependentSet(cell_marker.mesh(), coarsen_boundary);

  dolfin_assert( _forbidden_vertices.size() == cell_marker.mesh().numVertices() );
}
//-----------------------------------------------------------------------------
template<typename T>
void CoarseningManager::initCommon(MeshFunction<T>& cell_marker)
{
#ifdef ____USE_D_MESH____
    /// Delete old dmesh and import new mesh into dmesh
    if ( _dmesh )
      delete _dmesh;
    _dmesh = new DMesh;
    _dmesh->imp(cell_marker.mesh());
#endif 

  // extract boundary information
  findInteriorBoundaries(cell_marker.mesh());
  findDomainBoundaries(cell_marker.mesh());

  // build list of cells to coarsen
  findCellsToCoarsen(cell_marker);

  // initialize mapping
#ifndef ____USE_D_MESH____
  _vertex_map.init(cell_marker.mesh().numVertices());
  _cell_map.init(cell_marker.size());
#endif

  // reset request list
  _cells_to_request.clear();
  _vertices_to_request.clear();
}
//-----------------------------------------------------------------------------
void CoarseningManager::findInteriorBoundaries(Mesh& mesh)
{
  // set to false on whole domain
  _int_bnd_vertices.resize(mesh.numVertices());
  _int_bnd_vertices = false;
  _int_bnd_cells.resize(mesh.numCells());
  _int_bnd_cells = false;

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
    for ( CellIterator c_it(v) ; !c_it.end() ; ++c_it )
      _int_bnd_cells.at(c_it->index()) = true;
  }
}
//-----------------------------------------------------------------------------
void CoarseningManager::findDomainBoundaries(Mesh& mesh)
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
}
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
template<typename T>
void CoarseningManager::findCellsToCoarsen(MeshFunction<T>& cell_marker)
{
  _cells_to_coarsen.clear();

  for ( CellIterator c_it(cell_marker.mesh()) ; !c_it.end() ; ++c_it )
  {
    if ( cell_marker.get(c_it->index()) )
#ifdef ____USE_D_MESH____
      _cells_to_coarsen.push_back( _dmesh->getCell(c_it->index()) );
#else
      _cells_to_coarsen.push_back(c_it->index());
#endif
  }
}
//-----------------------------------------------------------------------------
#ifdef ____USE_D_MESH____
//-----------------------------------------------------------------------------
bool CoarseningManager::checkDCellNumbering(uint max_index)
{
  bool ret = true;
  for ( std::list<DCell *>::iterator c_it(_dmesh->cells.begin()) ; 
        c_it != _dmesh->cells.end() ; ++c_it )
  {
    if ( (*c_it)->id > max_index )
      ret = false;
  }
  return ret;
}
//-----------------------------------------------------------------------------
bool CoarseningManager::migrate(Mesh& /*omesh*/, bool repeat)
{
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();

  dolfin_assert( checkDCellNumbering(_bnd_cells.size() - 1) );

  if ( pe_size == 1 )
  {
    dolfin_assert(_cells_to_request.size() == 0);
    dolfin_assert(_vertices_to_request.size() == 0);
    return repeat;
  }

  message("[%d] Initializing migration", rank);

  // exchange maximum number of requested vertices
  int local_status[3], remote_status[3];
  local_status[0] = _vertices_to_request.size(); // local number of vertices
  local_status[1] = repeat;                      // local termination?
  local_status[2] = _migrated_cells;             // number of cells migrated
                                                 // in last migration
  MPI_Allreduce(&local_status, &remote_status, 3, MPI_INT, 
                MPI_MAX, MPI::DOLFIN_COMM);

  uint max_num_requested_vertices = remote_status[0]; // max number of vertices
  repeat = remote_status[1];                          // global termination?
  uint max_num_migrated_cells = remote_status[2];     // max num migrated cells

  message("[%d] max_num_requested_vertices = %d", rank, max_num_requested_vertices);

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

  dolfin_assert( checkDCellNumbering(_bnd_cells.size() - 1) );

  // Export DMesh to Mesh
  Array<int> old2new_cells(_bnd_cells.size());
  Array<int> old2new_vertices(_bnd_vertices.size());
  Mesh omesh;
  _dmesh->expKeepNumbering(omesh, &old2new_cells, &old2new_vertices);

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
    uint global_index = omesh.distdata().get_global(old2new_vertices[*it], 0);

    // vertex belongs to other process: request has to be distributed by owner
    if ( omesh.distdata().is_ghost(old2new_vertices[*it], 0) )
    {
      uint owner = omesh.distdata().get_owner(old2new_vertices[*it], 0);
      send_list_requests[owner].push_back(global_index);
    }
    // vertex belongs to this process: put request into map
    else
    {
      requested_vertices[global_index] = rank;
    }
  }

  message("[%d] Sending requests to owners", MPI::processNumber());

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
  MeshFunction<uint> partitions(omesh, omesh.topology().dim());
  partitions = rank;
  uint num_send_cells(0);

  // Build send list of vertices with requesting processes
  for ( std::map<uint,uint>::iterator m_it(requested_vertices.begin()) ; 
        m_it != requested_vertices.end() ; ++m_it )
  {
    uint local_index = omesh.distdata().get_local(m_it->first, 0);
    uint pe = m_it->second;

    // set of processes that share this vertex
    _set<uint> shared_adj = omesh.distdata().get_shared_adj(local_index, 0);
    for ( _set<uint>::iterator s_it(shared_adj.begin()) ;
          s_it != shared_adj.end() ; ++s_it )
    {
      if ( *s_it == rank ) continue;

      send_list_requests[*s_it].push_back(m_it->first);
      send_list_requests[*s_it].push_back(m_it->second);
    }

    // select lowest process index for all cells around requested vertex
    Vertex v(omesh, local_index);
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

  message("[%d] Receiving requests from owners", MPI::processNumber());

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
      uint local_index = omesh.distdata().get_local(recv_buff_requests[i], 0);
      Vertex v(omesh, local_index);

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
  MeshFunction<double> forbidden_vertices(omesh, 0);
  MeshFunction<double> forbidden_vertices_new;
  for ( uint old_idx(0) ; old_idx < _forbidden_vertices.size();
        ++old_idx )
  {
    int new_idx = old2new_vertices[old_idx];
    if ( new_idx < 0 )
      continue;
    forbidden_vertices.set(new_idx, double(isForbiddenVertex(old_idx)) );
  }
  vertex_functions.push_back( 
    std::make_pair(&forbidden_vertices, &forbidden_vertices_new) );

  // Prepare cell_marker for coarsening for exchange
  MeshFunction<uint> cell_marker(omesh, omesh.topology().dim());
  MeshFunction<uint> cell_marker_new;
  cell_marker = false;
  for ( List<uint>::iterator it(_cells_to_request.begin()) ; 
        it != _cells_to_request.end() ; ++it )
  {
    int new_idx = old2new_cells[*it];
    if ( new_idx < 0 ) 
      continue;
    cell_marker.set(new_idx, 1);
  }
  for ( List<DCell *>::iterator it(_cells_to_coarsen.begin()) ; 
        it != _cells_to_coarsen.end() ; ++it )
  {
    DCell * dc = *it;
    int new_idx = old2new_cells[dc->id];
    if ( new_idx < 0 ) 
      continue;
    cell_marker.set(new_idx, 1);
  }
  cell_functions.push_back( std::make_pair(&cell_marker, &cell_marker_new) );

  _migrated_cells = num_send_cells;
  message("[%d] Distribution. Sending %d cells", rank, num_send_cells);

  // distribute partitioning and MeshFunctions
  omesh.distribute( partitions, cell_functions, vertex_functions );

  // Re-initialize
  initCommon(cell_marker_new);

  // Update independent set
  //MeshFunctionConverter::cast(forbidden_vertices_new, _forbidden_vertices);
  _forbidden_vertices.resize(omesh.numVertices());
  _forbidden_vertices = false;
  for ( VertexIterator v_it(omesh) ; !v_it.end() ; ++v_it )
    if ( forbidden_vertices_new.get(v_it->index()) > 0.5 )
      _forbidden_vertices[v_it->index()] = true;

  message("[%d] Distribution done! Now %d cells are marked for coarsening!", 
            rank, _cells_to_coarsen.size());

  return true;
}
//-----------------------------------------------------------------------------
#else // ____USE_D_MESH____
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

  message("[%d] Initializing migration", rank);

  // exchange maximum number of requested vertices
  int local_status[3], remote_status[3];
  local_status[0] = _vertices_to_request.size(); // local number of vertices
  local_status[1] = repeat;                      // local termination?
  local_status[2] = _migrated_cells;             // number of cells migrated
                                                 // in last migration
  MPI_Allreduce(&local_status, &remote_status, 3, MPI_INT, 
                MPI_MAX, MPI::DOLFIN_COMM);

  uint max_num_requested_vertices = remote_status[0]; // max number of vertices
  repeat = remote_status[1];                          // global termination?
  uint max_num_migrated_cells = remote_status[2];     // max num migrated cells

  message("[%d] max_num_requested_vertices = %d", rank, max_num_requested_vertices);

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
  message("[%d] Distribution. Sending %d cells", rank, num_send_cells);

  // distribute partitioning and MeshFunctions
  mesh.distribute( partitions, cell_functions, vertex_functions );

  // Re-initialize
  initCommon(cell_marker_new);

  // Update independent set
  //MeshFunctionConverter::cast(forbidden_vertices_new, _forbidden_vertices);
  _forbidden_vertices.clear();
  _forbidden_vertices.resize(mesh.numVertices(), false);
  for ( VertexIterator v_it(mesh) ; !v_it.end() ; ++v_it )
    if ( forbidden_vertices_new.get(v_it->index()) > 0.5 )
      _forbidden_vertices[v_it->index()] = true;

  message("[%d] Distribution done! Now %d cells are marked for coarsening!", 
            rank, _cells_to_coarsen.size());

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