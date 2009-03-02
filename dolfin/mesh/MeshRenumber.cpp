// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//



#include "MeshRenumber.h"
#include "MeshDistributedData.h"
#include "MeshFunction.h"
#include "BoundaryMesh.h"
#include "MeshData.h"
#include "Edge.h"
#include "Cell.h"
#include "Face.h"
#include "Vertex.h"

#include <dolfin/common/Array.h>
#include <dolfin/main/MPI.h>

#include <mpi.h>
#include <cstdlib> 
#include <ctime> 

#include <map>



using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshRenumber::renumber(Mesh& mesh)
{
  renumber_vertices(mesh);
  renumber_edges(mesh);
  renumber_faces(mesh);
  renumber_cells(mesh);
}       		    
//-----------------------------------------------------------------------------
void MeshRenumber::renumber_vertices(Mesh& mesh)
{
  if(mesh.distdata()._valid_vertex_numbering || MPI::numProcesses() == 1)
    return;

  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();

  // Number of own vertices
  uint offset = 0;
  uint num_vert = mesh.numVertices() - mesh.distdata().num_ghost(0);

#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_vert, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Scan(&num_vert, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  offset -= num_vert;
#endif

  uint num_glb;  
  MPI_Allreduce(&num_vert, &num_glb, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  mesh.distdata().set_global_numVertices(num_glb);
  
  _map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numVertices(); i++){
    if(!mesh.distdata().is_ghost(i, 0)){
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;
    }
  }

  Array<uint> *ghost_buff = new Array<uint>[pe_size];
  for(MeshGhostIterator iter(mesh.distdata(), 0); !iter.end(); ++iter)
    ghost_buff[iter.owner()].push_back(mesh.distdata().get_global(iter.index(), 0)); 
 
  MPI_Status status;
  Array<uint> send_buff;
  uint src,dest;
  uint recv_size = mesh.distdata().num_ghost(0); 
  int recv_count, recv_size_gh, send_size;  
  
  for(uint i = 0; i < pe_size; i++) {
    send_size = ghost_buff[i].size();
    MPI_Reduce(&send_size, &recv_size_gh, 1, 
	       MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
  }

  uint *recv_ghost = new uint[ recv_size_gh];
  uint *recv_buff = new uint[ recv_size ];
  
  for(uint j=1; j < pe_size; j++){
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&ghost_buff[dest][0], ghost_buff[dest].size(),
		 MPI_UNSIGNED, dest, 1, recv_ghost, recv_size_gh, 
		 MPI_UNSIGNED, src, 1, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);

    for(int k=0; k < recv_count; k++)
      send_buff.push_back(new_global[mesh.distdata().get_local(recv_ghost[k], 0)]);

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, src, 2,
		 recv_buff, recv_size , MPI_UNSIGNED, dest, 2, 
		 MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);
    
    for(int j=0; j < recv_count; j++){
      new_global[mesh.distdata().get_local(ghost_buff[dest][j] , 0)] = recv_buff[j];
      new_local[recv_buff[j] ] = mesh.distdata().get_local(ghost_buff[dest][j], 0);
    }    
    send_buff.clear();
  }
  
  // Use new numbering
  mesh.distdata().local_indices[0] = new_local;
  mesh.distdata().global_indices[0] = new_global;  
  mesh.distdata()._valid_vertex_numbering = true;
  mesh.distdata().finialize(0);

  delete[] recv_buff;
  delete[] recv_ghost;
  for(uint i = 0; i < pe_size; i++)
    ghost_buff[i].clear();
  delete[] ghost_buff;
}
//-----------------------------------------------------------------------------
void MeshRenumber::renumber_edges(Mesh& mesh)
{

  if( mesh.distdata()._valid_edge_numbering || MPI::numProcesses() == 1)
    return;  

  // Flush shared/ghosted edges
  mesh.distdata().flush_edges();

  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  std::map<EdgeKey, uint> edge_map, edge_id;  
  EdgeKey key;
  Array<uint> send_buff, send_buff_id;
  _map<uint,uint> send_mapping;
  _set<uint> used_edge;

  srand((uint)time(0) + rank);
  for(MeshSharedIterator sv(mesh.distdata(), 0); !sv.end(); ++sv){
    Vertex v(mesh, sv.index());
    for(EdgeIterator e(v); !e.end(); ++e) {
      if( used_edge.count(e->index()) == 0) {
	const uint *edge_v = e->entities(0);
	key = edge_key(edge_v[0], edge_v[1]);
	edge_map[key] = e->index();
	edge_id[key] = (uint) rand();
	send_buff.push_back( mesh.distdata().get_global(edge_v[0], 0) );
	send_buff.push_back( mesh.distdata().get_global(edge_v[1], 0) );
	send_buff_id.push_back(edge_id[key]);
	used_edge.insert(e->index());
	mesh.distdata().set_shared(*e);
      }
    }    
  }


  uint num_ghost = 0;
  // Assign ownership of shared edges
  MPI_Status status;
  uint src,dest;
  int max_un, num_un, max_id, num_id, recv_count;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);
  num_id = send_buff_id.size();
  MPI_Allreduce(&num_id, &max_id, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);
  uint *recv_buff = new uint[max_un];
  uint *recv_buff_id = new uint[max_id];

  for(int j = 1 ; j < (int) pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    
    
    MPI_Sendrecv(&send_buff_id[0], num_id, MPI_UNSIGNED, dest, 1, 
		 recv_buff_id, max_id, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
    
    for(uint i =0; i < (uint) recv_count ; i += 2){
      // Check if I have the vertices
      if(mesh.distdata().have_global(recv_buff[i], 0) &&
	 mesh.distdata().have_global(recv_buff[i+1], 0)) {

	// Generate edge key
	key = edge_key(mesh.distdata().get_local(recv_buff[i], 0),
		       mesh.distdata().get_local(recv_buff[i+1], 0));

	// Check if I have the corresponding edge
	if(edge_id.count(key)) {
	  if( recv_buff_id[i>>1] < edge_id[key] ||
	      (recv_buff_id[i>>1] == edge_id[key] && status.MPI_SOURCE < rank)){
	    edge_id.erase(key);
	    mesh.distdata().set_ghost( edge_map[key], 1);
	    num_ghost++;
	  }
	}
      }
    }
  }
  
  dolfin_assert(num_ghost == mesh.distdata().num_ghost(1));

  // Number of own edges
  uint offset = 0;
  uint num_edges = mesh.numEdges() - mesh.distdata().num_ghost(1);

#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_edges, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Scan(&num_edges, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  offset -= num_edges;
#endif

  uint num_glb;  
  MPI_Allreduce(&num_edges, &num_glb, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numEdges(num_glb);

  send_buff.clear();
  //  send_buff.reserve(mesh.distdata().num_ghost(1));

  uint num = 0;
  _map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numEdges(); i++){
    if( !mesh.distdata().is_ghost(i, 1) ) { 
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;
    }
    else {
      Edge e(mesh, i);
      const uint *edge_v = e.entities(0);
      send_buff.push_back( mesh.distdata().get_global(edge_v[0], 0) );
      send_buff.push_back( mesh.distdata().get_global(edge_v[1], 0) );
      send_mapping[num++] = e.index();	
    }    

  }
  num_un = send_buff.size();

  //Exchange assigned global numbers
  Array<uint> global_buff;
  for(int j = 1 ; j < pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    
    
    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
    
    for(int i =0; i < recv_count ; i += 2){
      // Check if I have the vertices
      if(mesh.distdata().have_global(recv_buff[i], 0) &&
	 mesh.distdata().have_global(recv_buff[i+1], 0)) {
	
	// Generate edge key
	key = edge_key(mesh.distdata().get_local(recv_buff[i], 0),
		       mesh.distdata().get_local(recv_buff[i+1], 0));

	if(edge_id.count(key)){
	  global_buff.push_back(i>>1);
	  global_buff.push_back( new_global[ edge_map[key] ] );	  	  
	  mesh.distdata().set_shared( edge_map[key], 1);
	}
      }
    }
    
    MPI_Sendrecv(&global_buff[0], global_buff.size(), MPI_UNSIGNED, src, 2,
		 recv_buff, max_un, MPI_UNSIGNED, dest, 2,
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  

    for(int i = 0 ; i < recv_count; i += 2){
      new_global[ send_mapping[ recv_buff[i]] ] = recv_buff[i+1];
      new_local[ recv_buff[i+1] ] =  send_mapping[ recv_buff[i] ];
      mesh.distdata().set_ghost_owner(send_mapping[recv_buff[i]], 
				      status.MPI_SOURCE, 1);
    }
    global_buff.clear();
  }


  // Use new numbering
  mesh.distdata().local_indices[1] = new_local;
  mesh.distdata().global_indices[1] = new_global;
  mesh.distdata()._valid_edge_numbering =  true;
 
  delete[] recv_buff;
  delete[] recv_buff_id;

}
//-----------------------------------------------------------------------------
void MeshRenumber::renumber_faces(Mesh& mesh)
{

  if( mesh.distdata()._valid_face_numbering || 
      MPI::numProcesses() == 1 || mesh.topology().dim() == 2)
    return;  

  mesh.distdata().flush_faces();

  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  BoundaryMesh local_boundary;
  local_boundary.init_interior(mesh);
  //local_boundary.init_local(mesh);
  MeshFunction<uint>* cell_map = local_boundary.data().meshFunction("cell map");
 
  Array<uint> send_buff, send_buff_id;
  std::map<FaceKey, uint> face_map, face_id;
  std::map<uint,uint> send_mapping;
  FaceKey facekey;
  _set<uint> used_face;

  srand((uint)time(0) +  rank);

  for(CellIterator bf(local_boundary); !bf.end(); ++bf){
    Face f(mesh, cell_map->get(*bf));
    send_buffer_face(send_buff, mesh, f);
    facekey = face_key(f);
    face_map[facekey] = f.index();
    face_id[facekey] = (uint) rand();
    send_buff_id.push_back(face_id[facekey]);
    mesh.distdata().set_shared(f);
  }

  uint num_ghost = 0 ;
  Face f(mesh, 0);
  uint inc = 2 * f.numEntities(0);

  // Assign ownership of shared faces
  MPI_Status status;
  uint src,dest;
  int max_un, num_un, max_id, num_id, recv_count, recv_count_id;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);
  num_id = send_buff_id.size();
  MPI_Allreduce(&num_id, &max_id, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);
  uint *recv_buff = new uint[max_un];
  uint *recv_buff_id = new uint[max_id];
  EdgeKey key;
  for(int j = 1 ; j < (int) pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    
    
    MPI_Sendrecv(&send_buff_id[0], num_id, MPI_UNSIGNED, dest, 1, 
		 recv_buff_id, max_id, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count_id);  

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  

    uint ii = 0;
    for(uint i = 0; i < (uint) recv_count ; ii++, i += inc){    
      facekey.clear();
      uint num_ok = 0;
      for(uint k = 0; k < inc; k += 2){
	// Check if I have the vertices     
	if(mesh.distdata().have_global(recv_buff[i+k], 0) &&
	   mesh.distdata().have_global(recv_buff[i+k+1], 0)) {
	  // Generate edge key
	  key = edge_key(mesh.distdata().get_local(recv_buff[i+k], 0),
			 mesh.distdata().get_local(recv_buff[i+k+1], 0));
	  facekey.insert(key);
	  num_ok++;
	}
      }

      if(num_ok < 3)
	continue;
	
      // Check if I have the corresponding edge
      if(face_id.count(facekey)) {
	dolfin_assert(face_id.count(facekey));
	if( recv_buff_id[ii] < face_id[facekey] ||
	    (recv_buff_id[ii] == face_id[facekey] && status.MPI_SOURCE < rank)){
	  face_id.erase(facekey);
	  mesh.distdata().set_ghost( face_map[facekey], 2);
	  num_ghost++;
	}
      }
    }
  }


  dolfin_assert( num_ghost == mesh.distdata().num_ghost(2));


  // Number of own faces
  uint offset = 0;
  uint num_faces = mesh.numFaces() - mesh.distdata().num_ghost(2);

#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_faces, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Scan(&num_faces, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  offset -= num_faces;
#endif

  uint num_glb;  
  MPI_Allreduce(&num_faces, &num_glb, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numFaces(num_glb);  
 
  send_buff.clear();
  //  send_buff.reserve(mesh.distdata().num_ghost(2));

  uint num = 0;
  _map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numFaces(); i++){
    if( !mesh.distdata().is_ghost(i, 2) ){
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;
    }
    else {
      Face f(mesh, i);
      send_mapping[num++] = f.index();        
      send_buffer_face(send_buff, mesh, f);
    }
  }

  num_un = send_buff.size();
  //Exchange assigned global numbers
  Array<uint> global_buff;
  for(int j = 1 ; j < pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    
    
    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
    uint ii = 0;
    for(uint i = 0; i < (uint) recv_count ; ii++, i += inc){    
      // Check if I have the vertices     
      facekey.clear();
      uint num_ok = 0;
      for(uint k = 0; k < inc ; k += 2){
	if(mesh.distdata().have_global(recv_buff[i+k], 0) &&
	   mesh.distdata().have_global(recv_buff[i+k+1], 0)) {
	  // Generate edge key
	  key = edge_key(mesh.distdata().get_local(recv_buff[i+k], 0),
			 mesh.distdata().get_local(recv_buff[i+k+1], 0));

	  facekey.insert(key);
	  num_ok++;

	}
      }
      if(num_ok < 3)
	continue;
      
      if(face_id.count(facekey)){
	global_buff.push_back(ii);
	global_buff.push_back( new_global[ face_map[facekey] ] );
	mesh.distdata().set_shared(face_map[facekey], 2);
      }
    }
    
    MPI_Sendrecv(&global_buff[0], global_buff.size(), MPI_UNSIGNED, src, 2,
		 recv_buff, max_un, MPI_UNSIGNED, dest, 2,
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
    
    for(int i = 0 ; i < recv_count; i += 2){
      new_global[ send_mapping[ recv_buff[i]] ] = recv_buff[i+1];
      new_local[ recv_buff[i+1] ] =  send_mapping[ recv_buff[i] ];
      mesh.distdata().set_ghost_owner(send_mapping[ recv_buff[i] ], 
				      status.MPI_SOURCE, 2 );
    }
    global_buff.clear();
  }
  
  mesh.distdata().local_indices[2] = new_local;
  mesh.distdata().global_indices[2] = new_global;
  mesh.distdata()._valid_face_numbering =  true;

  delete[] recv_buff;
  delete[] recv_buff_id;

}
//-----------------------------------------------------------------------------
void MeshRenumber::renumber_cells(Mesh& mesh)
{
  if( mesh.distdata()._valid_cell_numbering || MPI::numProcesses() == 1)
    return;

  uint offset = 0;
  uint num_cells = mesh.numCells();

#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_cells, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Scan(&num_cells, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  offset -= num_cells;
#endif

  _map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numCells(); i++){
    new_global[i] = offset++;
    new_local [ new_global[i] ] = i;
  }

  uint num_glb;  
  MPI_Allreduce(&num_cells, &num_glb, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numCells(num_glb);
  
  // Use new numbering
  mesh.distdata().local_indices[3] = new_local;
  mesh.distdata().global_indices[3] = new_global;  
  mesh.distdata()._valid_cell_numbering = true;
  mesh.distdata().finialize(3);
 
}
//-----------------------------------------------------------------------------
std::pair<dolfin::uint, dolfin::uint> MeshRenumber::edge_key(uint id1,uint id2)
{
  if(id2 < id1){
    EdgeKey key(id2,id1);    
    return key;
  }
  else {
    EdgeKey key(id1,id2);    
    return key;
  }
  
}
//-----------------------------------------------------------------------------
std::set<std::pair<dolfin::uint,dolfin::uint> > MeshRenumber::face_key(Face& f)
{
  const uint *face_v = f.entities(0);
  FaceKey fk; 
  fk.insert( edge_key(face_v[0], face_v[1]) );
  fk.insert( edge_key(face_v[1], face_v[2]) );
  
  switch(f.numEntities(0)) 
    {
    case 3:
      fk.insert( edge_key(face_v[2], face_v[0]) );
      break;
    case 4:
      fk.insert( edge_key(face_v[2], face_v[3]) );
      fk.insert( edge_key(face_v[3], face_v[0]) );
      break;
    default:
      error("Unkown entity");
    }    
 
  return fk;
}
//-----------------------------------------------------------------------------
void MeshRenumber::send_buffer_face(Array<uint>& send_buff, Mesh& mesh, Face& f)
{
  const uint *face_v = f.entities(0);
  send_buff.push_back( mesh.distdata().get_global(face_v[0], 0) );
  send_buff.push_back( mesh.distdata().get_global(face_v[1], 0) );
    
  send_buff.push_back( mesh.distdata().get_global(face_v[1], 0) );
  send_buff.push_back( mesh.distdata().get_global(face_v[2], 0) );

  switch(f.numEntities(0)) 
    {
    case 3:
      send_buff.push_back( mesh.distdata().get_global(face_v[2], 0) );
      send_buff.push_back( mesh.distdata().get_global(face_v[0], 0) );	
      break;
    case 4:
      send_buff.push_back( mesh.distdata().get_global(face_v[2], 0) );
      send_buff.push_back( mesh.distdata().get_global(face_v[3], 0) );

      send_buff.push_back( mesh.distdata().get_global(face_v[3], 0) );
      send_buff.push_back( mesh.distdata().get_global(face_v[0], 0) );	
      break;
    default:
      error("Unkown entity");
    }          
}
//-----------------------------------------------------------------------------


