// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrøm, 2007.
// Modified by Anders Logg, 2007.
// Modified by Niclas Jansson, 2008.
// Modified by Balthasar Reuter, 2013.
//
// First added:  2007-05-30
// Last changed: 2008-01-18

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MPIMeshCommunicator.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <map>


#ifdef HAVE_MPI
   #include <mpi.h>
#endif

using namespace dolfin;

#ifdef HAVE_MPI

//-----------------------------------------------------------------------------
MPIMeshCommunicator::MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MPIMeshCommunicator::~MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh,
				     MeshFunction<uint>& distribution)
{
  distributeCommon(mesh, distribution, 0, 0);
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh, 
				       MeshFunction<uint>& distribution, 
				       MeshFunction<bool>& old_cell_marker,
				       MeshFunction<bool>& cell_marker) 
{
  distributeCommon(mesh, distribution, &old_cell_marker, &cell_marker);
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh, 
               MeshFunction<uint>& distribution,
               Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > 
                >& cell_functions)
{
  distributeCommon(mesh, distribution, cell_functions);
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distributeCommon(Mesh& mesh, 
               MeshFunction<uint>& distribution, 
               MeshFunction<bool>* old_cell_marker,
               MeshFunction<bool>* cell_marker) 
{

  MeshDistributedData distdata;  
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint gdim = mesh.geometry().dim();
  uint ndims = mesh.type().numVertices(mesh.topology().dim());

  Array<real> *send_list_vertices = new Array<real>[pe_size];
  Array<uint> *send_list_mappings = new Array<uint>[pe_size];
  Array<uint> *send_list_cells = new Array<uint>[pe_size];
  Array<real> coords;
  Array<uint> cl, shared_buffer;
  Array<bool> cm; 
  uint num_cells, num_vertices, target_proc, glb_index, offset;

  int recv_size,recv_size_cell,send_size;
  int recv_count,recv_count_vertices,recv_count_cells;

  // Process mesh entities according to distribution
  uint vi =0;
  // Distribution defined per vertex
  if(distribution.dim() == 0) {
    for (VertexIterator v(mesh); !v.end(); ++v){
      glb_index = mesh.distdata().get_global(*v); 
      if(!mesh.distdata().is_ghost(v->index(), 0)){
  if(distribution.get(*v) != rank) {
    target_proc = distribution.get(*v);
    send_list_mappings[target_proc].push_back(glb_index);   
    send_list_vertices[target_proc].push_back(v->point().x());
    send_list_vertices[target_proc].push_back(v->point().y());
    if(gdim > 2)
      send_list_vertices[target_proc].push_back(v->point().z());
  } 
  else if(distribution.get(*v) == rank) {
    coords.push_back(v->point().x());
    coords.push_back(v->point().y());
    if(gdim > 2)
      coords.push_back(v->point().z());
    distdata.set_map(vi++, glb_index, 0);
  }
      }
    }
    recv_count_cells = 0;
  }
  // Distribution defined per cell
  else if(distribution.dim() == mesh.topology().dim()) {

    MeshFunction<bool> vertex_used(mesh, 0);
    vertex_used = false;

    for(CellIterator c(mesh); !c.end(); ++c){      
      if(distribution.get(*c) != rank) {
  target_proc = distribution.get(*c);
  for(VertexIterator v(*c); !v.end(); ++v) {
    // Buffer cell global vertex indices 
    glb_index = mesh.distdata().get_global(*v);     
    send_list_cells[target_proc].push_back(glb_index);

    // Buffer all cell vertices that belong to another processor
      if(!mesh.distdata().is_ghost(v->index(), 0) && !vertex_used.get(*v)){
      send_list_mappings[target_proc].push_back(glb_index);    
      send_list_vertices[target_proc].push_back(v->point().x());
      send_list_vertices[target_proc].push_back(v->point().y());
      if(gdim > 2)
        send_list_vertices[target_proc].push_back(v->point().z());
      vertex_used.set(*v, true);
    }
  }

  // Transfer Cell marker  (mesh refinement)
  if( cell_marker ) {
    if( old_cell_marker->get(*c) ) 
      send_list_cells[target_proc].push_back(1);
    else 
      send_list_cells[target_proc].push_back(0);
  }
      }
      else {
  for(VertexIterator v(*c); !v.end(); ++v) {
    if(!vertex_used.get(*v)) {
      glb_index = mesh.distdata().get_global(*v); 
      if(!mesh.distdata().is_ghost(v->index(), 0)){   
        coords.push_back(v->point().x());
        coords.push_back(v->point().y());
        if(gdim > 2)
    coords.push_back(v->point().z());
        distdata.set_map(vi++, glb_index, 0);    
        vertex_used.set(*v, true);        
      }
    }
  }
      }
    }
    
    recv_count_cells = 0;
    for(uint i=0; i<pe_size; i++){
      send_size = send_list_cells[i].size();
      MPI_Reduce(&send_size, &recv_count_cells, 1, MPI_INT, 
     MPI_SUM, i, MPI::DOLFIN_COMM);
    }
  }
  else
    error("Distribution defined on unkown mesh entity");

  // Exchange the processed entities
  recv_count = 0;
  for(uint i=0; i<pe_size; i++){
    send_size = send_list_vertices[i].size();
    MPI_Reduce(&send_size, &recv_count, 1, MPI_INT,MPI_MAX, i, MPI::DOLFIN_COMM);
  }
  recv_count_vertices = recv_count / gdim;
  num_vertices = recv_count;

  double *recv_buff = new double[recv_count];
  uint *recv_buff_map = new uint[recv_count_vertices];

  num_cells = recv_count_cells;
  uint *recv_buff_cell = new uint[recv_count_cells];
  uint *rcp = &recv_buff_cell[0];
  
  MPI_Status status;
  uint src, dest, buff_map;
  for(uint i=1; i<pe_size; i++){
    
    src = (rank - i + pe_size) % pe_size;
    dest = (rank + i) % pe_size;

    MPI_Sendrecv(&send_list_vertices[dest][0], send_list_vertices[dest].size(),
     MPI_DOUBLE, dest, 0, recv_buff, recv_count, MPI_DOUBLE, src, 
     0,MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_DOUBLE,&recv_size);


    MPI_Sendrecv(&send_list_mappings[dest][0], send_list_mappings[dest].size(),
     MPI_UNSIGNED, dest, 1, recv_buff_map, recv_count_vertices, 
     MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_list_cells[dest][0], send_list_cells[dest].size(), 
     MPI_UNSIGNED, dest, 2, rcp, recv_count_cells, MPI_UNSIGNED,
     src, 2, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_size_cell);
    rcp += recv_size_cell;
    recv_count_cells -= recv_size_cell;

    buff_map = 0;
    for(int i=0; i < recv_size; i +=gdim){
      if(!distdata.have_global( recv_buff_map[buff_map], 0)) {
  distdata.set_map(vi++, recv_buff_map[buff_map], 0);
  coords.push_back(recv_buff[i]);
  coords.push_back(recv_buff[i+1]);
  if( gdim > 2)
    coords.push_back(recv_buff[i+2]);
      }
      buff_map++;
    }
  }

  //Clear send buffers
  for(uint i=0; i<pe_size; i++){
    send_list_cells[i].clear();
    send_list_vertices[i].clear();
    send_list_mappings[i].clear();
  }
  delete[] send_list_vertices;
  delete[] send_list_mappings;
  delete[] send_list_cells;
  delete[] recv_buff_map;
  delete[] recv_buff;

  // Process new and old cells if distribution is defined on cells
  if(distribution.dim() == mesh.topology().dim()) {

    // Add old cells
    for(CellIterator c(mesh); !c.end(); ++c){
      if(distribution.get(*c) == rank) {
  for(VertexIterator v(*c); !v.end(); ++v){
    glb_index = mesh.distdata().get_global(*v);   
    if(!distdata.have_global(glb_index, 0)){
      cl.push_back(vi);
      for(uint j=0;j<gdim;j++)
        coords.push_back(0.0);
      distdata.set_map(vi, glb_index, 0);
      distdata.set_ghost(vi++, 0);
      shared_buffer.push_back(glb_index);
    }
    else 
      cl.push_back(distdata.get_local(glb_index, 0));
  }

  // Mark cell for refinement
  if( cell_marker )
    cm.push_back(old_cell_marker->get(*c));      
      }
    }

    // Add new cells
    uint cell_n = 0;
    for(uint i=0; i <num_cells; i++){
      if( cell_marker && cell_n == ndims){
  cm.push_back( recv_buff_cell[i] == 1 );
  cell_n = 0;
      }
      else {
  if(distdata.have_global(recv_buff_cell[i], 0))
    cl.push_back(distdata.get_local(recv_buff_cell[i], 0));
  else{
    cl.push_back(vi);
    for(uint j=0;j<gdim;j++)
      coords.push_back(0.0);
    distdata.set_map(vi, recv_buff_cell[i], 0);
    distdata.set_ghost(vi++, 0);
    shared_buffer.push_back(recv_buff_cell[i]);
  }
  cell_n++;
      }      
    }

    // Exchange ghosted entities
    Array<real> send_buff;
    Array<uint> send_buff_indices, recv_source;
    send_size = shared_buffer.size();
    recv_count_vertices = static_cast<uint>(gdim) * send_size;
    recv_buff = new double[recv_count_vertices];
    double *rp = &recv_buff[0];
    recv_buff_map = new uint[send_size];
    uint *rmp = &recv_buff_map[0];        
    MPI_Allreduce(&send_size, &recv_count, 1, MPI_INT,MPI_MAX, MPI::DOLFIN_COMM);
    uint *shared = new uint[recv_count];
    for(uint i=1; i<pe_size; i++){

      src = (rank - i + pe_size) % pe_size;
      dest = (rank + i) % pe_size;

      MPI_Sendrecv(&shared_buffer[0], shared_buffer.size(), MPI_UNSIGNED, dest,
       1, shared, recv_count, MPI_UNSIGNED, src ,1, MPI::DOLFIN_COMM,
       &status);
      MPI_Get_count(&status,MPI_UNSIGNED,&recv_size);
      
      for(int j=0; j<recv_size; j++)
  if(distdata.have_global(shared[j], 0) &&
     !distdata.is_ghost(distdata.get_local(shared[j], 0), 0)){
    offset = distdata.get_local(shared[j], 0) * gdim;
    send_buff.push_back(coords[offset]);
    send_buff.push_back(coords[offset + 1]);
    if(gdim >2)
      send_buff.push_back(coords[offset + 2]);
    send_buff_indices.push_back(shared[j]);
    distdata.set_shared(distdata.get_local(shared[j], 0), 0);
    distdata.set_shared_adj(distdata.get_local(shared[j], 0), src, 0);
  }

      MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_DOUBLE, src, 2,
       rp, recv_count_vertices, MPI_DOUBLE, dest, 2, 
       MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status,MPI_DOUBLE,&recv_size);

      rp += recv_size;
      recv_count_vertices -= recv_size;

      MPI_Sendrecv(&send_buff_indices[0], send_buff_indices.size(), 
       MPI_UNSIGNED, src, 3, rmp, send_size, MPI_UNSIGNED, 
       dest, 3, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status,MPI_UNSIGNED,&recv_size);  
      
      rmp += recv_size;
      send_size -=recv_size;
      
      for(int k=0; k < recv_size; k++)
  recv_source.push_back(status.MPI_SOURCE);

      send_buff.clear();
      send_buff_indices.clear();
    }        

    uint j=0;
    for(uint i=0; i < shared_buffer.size(); i++){
      offset = distdata.get_local(recv_buff_map[i], 0) * gdim;
      coords[offset] = recv_buff[j];
      coords[offset+1] = recv_buff[j+1];
      if(gdim > 2)
  coords[offset+2] = recv_buff[j+2];
      j += gdim;
      distdata.set_ghost_owner( distdata.get_local(recv_buff_map[i], 0), 
        recv_source[i], 0);
    }
    shared_buffer.clear();
    recv_source.clear();
    delete[] recv_buff_map;
    delete[] recv_buff;
    delete[] shared;
      
  }
  delete[] recv_buff_cell;        


  num_vertices = coords.size() / gdim ;
  num_cells = cl.size() / ndims;

  // Construct new mesh and add all buffered entities
  Mesh new_mesh;
  MeshEditor editor;
  editor.open(new_mesh, mesh.type().cellType(),
        mesh.topology().dim(), mesh.geometry().dim());

  distdata.set_global_numVertices(mesh.distdata().global_numVertices());
  distdata.set_global_numEdges(mesh.distdata().global_numEdges());
  distdata.set_global_numFaces(mesh.distdata().global_numFaces());
  distdata.set_global_numCells(mesh.distdata().global_numCells());

  editor.initVertices(num_vertices);
  editor.initCells(num_cells);

  vi=0;
  for(uint i=0; i<coords.size(); i +=gdim)
    switch(gdim){
    case 2:
      editor.addVertex(vi++,coords[i],coords[i+1]); break;
    case 3:
      editor.addVertex(vi++,coords[i],coords[i+1], coords[i+2]); break;
    }
  coords.clear();
  
  dolfin_assert(cl.size()%ndims == 0);
  uint ci=0;
  for(uint i=0; i <cl.size(); i+=ndims) {
    switch(ndims){
    case 3:
      editor.addCell(ci++,cl[i],cl[i+1],cl[i+2]); break;
    case 4:
      editor.addCell(ci++,cl[i],cl[i+1],cl[i+2],cl[i+3]); break;
    }
  }
  cl.clear();
  editor.close();

  // Overwrite old mesh with new, and invalidate numbering
  new_mesh._distdata = distdata;
  mesh = new_mesh;
  mesh.distdata().invalid_numbering();
  mesh.distdata().invalid_ownership();

  // Mark new cells for refinement
  if( cell_marker ) {
    dolfin_assert( cm.size() == mesh.numCells());
    cell_marker->init(mesh, mesh.topology().dim());
    for(uint i =0; i < cm.size(); i++)
      cell_marker->set( i , cm[i] );
    cm.clear();
  }
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distributeCommon(Mesh& mesh, 
               MeshFunction<uint>& distribution, 
               Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > 
                >& cell_functions) 
{
  typedef Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > > 
    MFArrayType;
  typedef 
    Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > >::iterator 
    MFArrayIteratorType;

  MeshDistributedData distdata;  
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint gdim = mesh.geometry().dim();
  uint ndims = mesh.type().numVertices(mesh.topology().dim());
  uint nfunctions = cell_functions.size();

  Array<real> *send_list_vertices = new Array<real>[pe_size];
  Array<uint> *send_list_mappings = new Array<uint>[pe_size];
  Array<uint> *send_list_cells = new Array<uint>[pe_size];
  Array<real> coords;
  Array<uint> cl, shared_buffer;
  Array<uint> functions[nfunctions];
  uint num_cells, num_vertices, target_proc, glb_index, offset;

  int recv_size,recv_size_cell,send_size;
  int recv_count,recv_count_vertices,recv_count_cells;

  // Process mesh entities according to distribution
  uint vi =0;
  // Distribution defined per vertex
  if(distribution.dim() == 0) {
    for (VertexIterator v(mesh); !v.end(); ++v){
      glb_index = mesh.distdata().get_global(*v); 
      if(!mesh.distdata().is_ghost(v->index(), 0)){
  if(distribution.get(*v) != rank) {
    target_proc = distribution.get(*v);
    send_list_mappings[target_proc].push_back(glb_index);   
    send_list_vertices[target_proc].push_back(v->point().x());
    send_list_vertices[target_proc].push_back(v->point().y());
    if(gdim > 2)
      send_list_vertices[target_proc].push_back(v->point().z());
  } 
  else if(distribution.get(*v) == rank) {
    coords.push_back(v->point().x());
    coords.push_back(v->point().y());
    if(gdim > 2)
      coords.push_back(v->point().z());
    distdata.set_map(vi++, glb_index, 0);
  }
      }
    }
    recv_count_cells = 0;
  }
  // Distribution defined per cell
  else if(distribution.dim() == mesh.topology().dim()) {

    MeshFunction<bool> vertex_used(mesh, 0);
    vertex_used = false;

    for(CellIterator c(mesh); !c.end(); ++c)
    {      
      if(distribution.get(*c) != rank) 
      {
        target_proc = distribution.get(*c);
        for(VertexIterator v(*c); !v.end(); ++v) 
        {
          // Buffer cell global vertex indices 
          glb_index = mesh.distdata().get_global(*v);     
          send_list_cells[target_proc].push_back(glb_index);

          // Buffer all cell vertices that belong to another processor
          if(!mesh.distdata().is_ghost(v->index(), 0) && !vertex_used.get(*v))
          {
            send_list_mappings[target_proc].push_back(glb_index);    
            send_list_vertices[target_proc].push_back(v->point().x());
            send_list_vertices[target_proc].push_back(v->point().y());
            if(gdim > 2)
              send_list_vertices[target_proc].push_back(v->point().z());
            vertex_used.set(*v, true);
          }
        }

        // Transfer functions
        if( nfunctions > 0 ) 
        {
          for ( MFArrayIteratorType f_it(cell_functions.begin()) ; 
                f_it != cell_functions.end() ; ++f_it )
          {
            send_list_cells[target_proc].push_back(f_it->first->get(*c));
          }
        }
      }
      else 
      {
        for(VertexIterator v(*c); !v.end(); ++v) 
        {
          if(!vertex_used.get(*v)) 
          {
            glb_index = mesh.distdata().get_global(*v); 
            if(!mesh.distdata().is_ghost(v->index(), 0))
            {   
              coords.push_back(v->point().x());
              coords.push_back(v->point().y());
              if(gdim > 2)
                coords.push_back(v->point().z());
              distdata.set_map(vi++, glb_index, 0);    
              vertex_used.set(*v, true);        
            }
          }
        }
      }
    }
    
    recv_count_cells = 0;
    for(uint i=0; i<pe_size; i++)
    {
      send_size = send_list_cells[i].size();
      MPI_Reduce(&send_size, &recv_count_cells, 1, MPI_INT, 
      MPI_SUM, i, MPI::DOLFIN_COMM);
    }
  }
  else
    error("Distribution defined on unkown mesh entity");

  // Exchange the processed entities
  recv_count = 0;
  for(uint i=0; i<pe_size; i++)
  {
    send_size = send_list_vertices[i].size();
    MPI_Reduce(&send_size, &recv_count, 1, MPI_INT,MPI_MAX, i, MPI::DOLFIN_COMM);
  }
  recv_count_vertices = recv_count / gdim;
  num_vertices = recv_count;

  double *recv_buff = new double[recv_count];
  uint *recv_buff_map = new uint[recv_count_vertices];

  num_cells = recv_count_cells;
  uint *recv_buff_cell = new uint[recv_count_cells];
  uint *rcp = &recv_buff_cell[0];
  
  MPI_Status status;
  uint src, dest, buff_map;
  for(uint i=1; i<pe_size; i++)
  {
    src = (rank - i + pe_size) % pe_size;
    dest = (rank + i) % pe_size;

    MPI_Sendrecv(&send_list_vertices[dest][0], send_list_vertices[dest].size(),
      MPI_DOUBLE, dest, 0, recv_buff, recv_count, MPI_DOUBLE, src, 
      0,MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_DOUBLE,&recv_size);


    MPI_Sendrecv(&send_list_mappings[dest][0], send_list_mappings[dest].size(),
      MPI_UNSIGNED, dest, 1, recv_buff_map, recv_count_vertices, 
      MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_list_cells[dest][0], send_list_cells[dest].size(), 
      MPI_UNSIGNED, dest, 2, rcp, recv_count_cells, MPI_UNSIGNED,
      src, 2, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_size_cell);
    rcp += recv_size_cell;
    recv_count_cells -= recv_size_cell;

    buff_map = 0;
    for(int i=0; i < recv_size; i +=gdim)
    {
      if(!distdata.have_global( recv_buff_map[buff_map], 0)) 
      {
        distdata.set_map(vi++, recv_buff_map[buff_map], 0);
        coords.push_back(recv_buff[i]);
        coords.push_back(recv_buff[i+1]);
      if( gdim > 2)
        coords.push_back(recv_buff[i+2]);
      }
      buff_map++;
    }
  }

  //Clear send buffers
  for(uint i=0; i<pe_size; i++){
    send_list_cells[i].clear();
    send_list_vertices[i].clear();
    send_list_mappings[i].clear();
  }
  delete[] send_list_vertices;
  delete[] send_list_mappings;
  delete[] send_list_cells;
  delete[] recv_buff_map;
  delete[] recv_buff;

  // Process new and old cells if distribution is defined on cells
  if(distribution.dim() == mesh.topology().dim()) 
  {
    // Add old cells
    for(CellIterator c(mesh); !c.end(); ++c)
    {
      if(distribution.get(*c) == rank) 
      {
        for(VertexIterator v(*c); !v.end(); ++v)
        {
          glb_index = mesh.distdata().get_global(*v);   
          if(!distdata.have_global(glb_index, 0))
          {
            cl.push_back(vi);
            for(uint j=0;j<gdim;j++)
              coords.push_back(0.0);
            distdata.set_map(vi, glb_index, 0);
            distdata.set_ghost(vi++, 0);
            shared_buffer.push_back(glb_index);
          }
          else 
            cl.push_back(distdata.get_local(glb_index, 0));
        }

        // Conserve function values
        if( nfunctions > 0 )
        {
          uint f_id(0);
          for ( MFArrayIteratorType f_it(cell_functions.begin()) ; 
                f_it != cell_functions.end() ; ++f_it, ++f_id )
          {
            functions[f_id].push_back(f_it->first->get(*c));
          }
        }
      }
    }

    // Add new cells
    uint cell_n = 0;
    for(uint i=0; i < num_cells; i++)
    {
      if( nfunctions > 0 && cell_n >= ndims)
      {
        functions[cell_n - ndims].push_back(recv_buff_cell[i]);
        ++cell_n;
        if ( (cell_n - ndims) == nfunctions)
          cell_n = 0;
      }
      else 
      {
        if(distdata.have_global(recv_buff_cell[i], 0))
          cl.push_back(distdata.get_local(recv_buff_cell[i], 0));
        else
        {
          cl.push_back(vi);
          for(uint j=0;j<gdim;j++)
            coords.push_back(0.0);
          distdata.set_map(vi, recv_buff_cell[i], 0);
          distdata.set_ghost(vi++, 0);
          shared_buffer.push_back(recv_buff_cell[i]);
        }
        cell_n++;
      }
    }

    // Exchange ghosted entities
    Array<real> send_buff;
    Array<uint> send_buff_indices, recv_source;
    send_size = shared_buffer.size();
    recv_count_vertices = static_cast<uint>(gdim) * send_size;
    recv_buff = new double[recv_count_vertices];
    double *rp = &recv_buff[0];
    recv_buff_map = new uint[send_size];
    uint *rmp = &recv_buff_map[0];        
    MPI_Allreduce(&send_size, &recv_count, 1, MPI_INT,MPI_MAX, MPI::DOLFIN_COMM);
    uint *shared = new uint[recv_count];
    for(uint i=1; i<pe_size; i++)
    {

      src = (rank - i + pe_size) % pe_size;
      dest = (rank + i) % pe_size;

      MPI_Sendrecv(&shared_buffer[0], shared_buffer.size(), MPI_UNSIGNED, dest,
       1, shared, recv_count, MPI_UNSIGNED, src ,1, MPI::DOLFIN_COMM,
       &status);
      MPI_Get_count(&status,MPI_UNSIGNED,&recv_size);
      
      for(int j=0; j<recv_size; j++)
        if( distdata.have_global(shared[j], 0) &&
           !distdata.is_ghost(distdata.get_local(shared[j], 0), 0) )
        {
          offset = distdata.get_local(shared[j], 0) * gdim;
          send_buff.push_back(coords[offset]);
          send_buff.push_back(coords[offset + 1]);
          if(gdim >2)
            send_buff.push_back(coords[offset + 2]);
          send_buff_indices.push_back(shared[j]);
          distdata.set_shared(distdata.get_local(shared[j], 0), 0);
          distdata.set_shared_adj(distdata.get_local(shared[j], 0), src, 0);
        }

      MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_DOUBLE, src, 2,
       rp, recv_count_vertices, MPI_DOUBLE, dest, 2, 
       MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status,MPI_DOUBLE,&recv_size);

      rp += recv_size;
      recv_count_vertices -= recv_size;

      MPI_Sendrecv(&send_buff_indices[0], send_buff_indices.size(), 
       MPI_UNSIGNED, src, 3, rmp, send_size, MPI_UNSIGNED, 
       dest, 3, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status,MPI_UNSIGNED,&recv_size);  
      
      rmp += recv_size;
      send_size -=recv_size;
      
      for(int k=0; k < recv_size; k++)
        recv_source.push_back(status.MPI_SOURCE);

      send_buff.clear();
      send_buff_indices.clear();
    }

    uint j=0;
    for(uint i=0; i < shared_buffer.size(); i++)
    {
      offset = distdata.get_local(recv_buff_map[i], 0) * gdim;
      coords[offset] = recv_buff[j];
      coords[offset+1] = recv_buff[j+1];
      if(gdim > 2)
        coords[offset+2] = recv_buff[j+2];
      j += gdim;
      distdata.set_ghost_owner( distdata.get_local(recv_buff_map[i], 0), 
                                recv_source[i], 0);
    }
    shared_buffer.clear();
    recv_source.clear();
    delete[] recv_buff_map;
    delete[] recv_buff;
    delete[] shared;
      
  }
  delete[] recv_buff_cell;        


  num_vertices = coords.size() / gdim ;
  num_cells = cl.size() / ndims;

  // Construct new mesh and add all buffered entities
  Mesh new_mesh;
  MeshEditor editor;
  editor.open(new_mesh, mesh.type().cellType(),
              mesh.topology().dim(), mesh.geometry().dim());

  distdata.set_global_numVertices(mesh.distdata().global_numVertices());
  distdata.set_global_numEdges(mesh.distdata().global_numEdges());
  distdata.set_global_numFaces(mesh.distdata().global_numFaces());
  distdata.set_global_numCells(mesh.distdata().global_numCells());

  editor.initVertices(num_vertices);
  editor.initCells(num_cells);

  vi=0;
  for(uint i=0; i<coords.size(); i +=gdim)
    switch(gdim)
    {
    case 2:
      editor.addVertex(vi++,coords[i],coords[i+1]); break;
    case 3:
      editor.addVertex(vi++,coords[i],coords[i+1], coords[i+2]); break;
    }
  coords.clear();
  
  dolfin_assert(cl.size()%ndims == 0);
  uint ci=0;
  for(uint i=0; i <cl.size(); i+=ndims) 
  {
    switch(ndims)
    {
    case 3:
      editor.addCell(ci++,cl[i],cl[i+1],cl[i+2]); break;
    case 4:
      editor.addCell(ci++,cl[i],cl[i+1],cl[i+2],cl[i+3]); break;
    }
  }
  cl.clear();
  editor.close();

  // Overwrite old mesh with new, and invalidate numbering
  new_mesh._distdata = distdata;
  mesh = new_mesh;
  mesh.distdata().invalid_numbering();
  mesh.distdata().invalid_ownership();

  // Fill new cell functions
  if( nfunctions > 0 ) 
  {
    uint f_id(0);
    for ( MFArrayIteratorType f_it(cell_functions.begin()) ; 
          f_it != cell_functions.end() ; ++f_it, ++f_id )
    {
      dolfin_assert( functions[f_id].size() == mesh.numCells() );
      f_it->second->init(mesh, mesh.topology().dim());
      for(uint i =0; i < functions[f_id].size(); ++i)
        f_it->second->set( i , functions[f_id][i] );
      functions[f_id].clear();
    }
  }
}
//-----------------------------------------------------------------------------

#else

//-----------------------------------------------------------------------------
MPIMeshCommunicator::MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MPIMeshCommunicator::~MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh,
				     MeshFunction<uint>& distribution)
{

  error("Cannot distribute mesh without MPI.");  
}
//-----------------------------------------------------------------------------
  void MPIMeshCommunicator::distribute(Mesh& mesh, 
				       MeshFunction<uint>& distribution, 
				       MeshFunction<bool>& old_cell_marker,
				       MeshFunction<bool>& cell_marker) 
{
  error("Cannot distribute mesh without MPI.");  
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh, 
               MeshFunction<uint>& distribution,
               Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > 
                >& cell_functions)
{
  error("Cannot distribute mesh without MPI.");  
}
//-----------------------------------------------------------------------------

#endif
