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
#include <map>
#include <cstdlib> 
#include <ctime> 

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
  uint* num_vert = new uint[ pe_size ];

  if(mesh.distdata()._size != mesh.distdata().global_vertex_indices.size())
    error("invalid size");
  // Number of own vertices
  num_vert[ rank ] = mesh.distdata()._size - mesh.distdata().num_ghost();
  MPI_Allgather(&num_vert[ rank ], 1, MPI_UNSIGNED,
		num_vert, 1, MPI_UNSIGNED, MPI_COMM_WORLD);  

  uint offset = 0;
  for(uint i = 1; i < rank+1; i++)
    offset += num_vert[i-1];

  std::map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.distdata()._size; i++){
    if(!mesh.distdata().is_ghost(i)){
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;
    }
  }

  Array<uint> *ghost_buff = new Array<uint>[pe_size];
  for(MeshGhostIterator iter(mesh.distdata()); !iter.end(); ++iter)
    ghost_buff[iter.owner()].push_back(mesh.distdata().get_global(iter.index(), 0)); 
 
  MPI_Status status;
  Array<uint> send_buff;
  int gh_count = static_cast<int>(mesh.distdata().num_ghost());
  uint src,dest;
  uint recv_size = gh_count;
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
  mesh.distdata().local_vertex_indices = new_local;
  mesh.distdata().global_vertex_indices = new_global;  
  mesh.distdata()._valid_vertex_numbering = true;
  
  delete[] recv_buff;
  delete[] recv_ghost;
  for(uint i = 0; i < pe_size; i++)
    ghost_buff[i].clear();
  delete[] ghost_buff;
  delete[] num_vert;    
}
//-----------------------------------------------------------------------------
void MeshRenumber::renumber_edges(Mesh& mesh)
{
  if( mesh.distdata()._valid_edge_numbering || MPI::numProcesses() == 1)
    return;
  
  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  MeshFunction<bool> ghosted_edge(mesh, 1);
  //  ghosted_edge.init(mesh, 1);
  ghosted_edge = false;
  uint num_ghosts = 0;
  
  std::map<EdgeKey, uint> edge_map, edge_id, owns_edge;  

  Array<uint> send_buff, send_buff_id;
  std::map<uint,uint> send_mapping;
  uint num = 0;

  srand((uint)time(0));
  for(MeshSharedIterator sv(mesh.distdata()); !sv.end(); ++sv){
    Vertex v(mesh, sv.index());
    for(EdgeIterator e(v); !e.end(); ++e) {
        const uint *edge_v = e->entities(0);
	uint edge_vert[2];
	edge_vert[0] = mesh.distdata().get_global(edge_v[0], 0);
	edge_vert[1] = mesh.distdata().get_global(edge_v[1], 0);
	EdgeKey key = edge_key(edge_v[0], edge_v[1]);
	edge_map[key] = e->index();
	edge_id[key] = (uint) rand + rank;
	send_buff.push_back( edge_vert[0]);
	send_buff.push_back( edge_vert[1]);
	send_buff_id.push_back(edge_id[key]);
	owns_edge[key] = true;
	send_mapping[num++] = e->index();
    }    
  }

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
  EdgeKey key;

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
      if(mesh.distdata().have_global(recv_buff[i]) &&
	 mesh.distdata().have_global(recv_buff[i+1])) {

	// Generate edge key
	key = edge_key(mesh.distdata().get_local(recv_buff[i], 0),
		       mesh.distdata().get_local(recv_buff[i+1], 0));

	// Check if I have the corresponding edge
	if(edge_id.count(key)) {
	  if( recv_buff_id[i>>1] < edge_id[key] ||
	      recv_buff_id[i>>1] == edge_id[key] && status.MPI_SOURCE < rank){
	    edge_id.erase(key);
	    owns_edge[key] = false;
	    ghosted_edge.set( edge_map[key], true);
	    num_ghosts++;
	  }
	}
      }
    }
  }
  
  // Number of own edges
  uint* num_edges= new uint[ pe_size ];
  num_edges[ rank ] = mesh.numEdges() - num_ghosts;
  MPI_Allgather(&num_edges[ rank ], 1, MPI_UNSIGNED, num_edges, 1,
		MPI_UNSIGNED, MPI_COMM_WORLD);
  
  uint num_glb;  
  MPI_Allreduce(&num_edges[ rank ], &num_glb, 1,
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numEdges(num_glb);

  uint offset = 0;
  for(int i = 1; i < rank+1; i++)
      offset += num_edges[i-1];

  std::map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numEdges(); i++){
    if( !ghosted_edge.get(i) ) {
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;
    }
  }

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
      if(mesh.distdata().have_global(recv_buff[i]) &&
	 mesh.distdata().have_global(recv_buff[i+1])) {
	
	// Generate edge key
	key = edge_key(mesh.distdata().get_local(recv_buff[i], 0),
		       mesh.distdata().get_local(recv_buff[i+1], 0));

	if(owns_edge[key] && owns_edge.count(key) == 1) {
	  global_buff.push_back(i>>1);
	  global_buff.push_back( new_global[ edge_map[key] ] );	  	  
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
     }
     global_buff.clear();
  }


  // Use new numbering
 mesh.distdata().local_edge_indices = new_local;
 mesh.distdata().global_edge_indices = new_global;
 mesh.distdata()._valid_edge_numbering =  true;
 
 delete[] num_edges;
 delete[] recv_buff;
 delete[] recv_buff_id;

}
//-----------------------------------------------------------------------------
void MeshRenumber::renumber_faces(Mesh& mesh)
{
  if( mesh.distdata()._valid_face_numbering || 
      MPI::numProcesses() == 1 || mesh.topology().dim() == 2)
    return;  

  std::map<FaceKey, uint> face_map, face_id, owns_face;
  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  MeshFunction<bool> ghosted_face;
  ghosted_face.init(mesh, 2);
  ghosted_face =  false;

  BoundaryMesh local_boundary;
  local_boundary.init_local(mesh);
  MeshFunction<uint>* cell_map = local_boundary.data().meshFunction("cell map");

  uint num_ghosts = 0 ;
  Array<uint> send_buff, send_buff_id;
  std::map<uint,uint> send_mapping;
  uint num = 0;
  FaceKey face_key;
  srand((uint)time(0));
  for(CellIterator bf(local_boundary); !bf.end(); ++bf){
    Face f(mesh, cell_map->get(*bf));
    face_key.clear();
    for(EdgeIterator e(f); !e.end(); ++e) {
      const uint *edge_v = e->entities(0);
      uint edge_vert[2];
      edge_vert[0] = mesh.distdata().get_global(edge_v[0], 0);
      edge_vert[1] = mesh.distdata().get_global(edge_v[1], 0);
      EdgeKey key = edge_key(edge_v[0], edge_v[1]);
      send_buff.push_back( edge_vert[0]);
      send_buff.push_back( edge_vert[1]);
      face_key.insert(key);
    }
    if(face_map.count(face_key) == 0) {
      face_map[face_key] = f.index();
      face_id[face_key] = (uint) rand() +  rank;
      owns_face[face_key] = true;
      send_mapping[num++] = f.index();
      send_buff_id.push_back(face_id[face_key]);
    }
  }

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
      // Check if I have the vertices     
      face_key.clear();

      uint num_ok = 0;
      for(uint k = 0; k < inc; k += 2){
	if(mesh.distdata().have_global(recv_buff[i+k]) &&
	   mesh.distdata().have_global(recv_buff[i+k+1])) {
	  // Generate edge key
	  key = edge_key(mesh.distdata().get_local(recv_buff[i+k], 0),
			 mesh.distdata().get_local(recv_buff[i+k+1], 0));
	  face_key.insert(key);
	  num_ok++;
	}
      }

      if(num_ok < 3)
	continue;
	
      // Check if I have the corresponding edge
      if(face_id.count(face_key)) {
	if( recv_buff_id[ii] < face_id[face_key] ||
	    recv_buff_id[ii] == face_id[face_key] && status.MPI_SOURCE < rank){
	  face_id.erase(face_key);
	  owns_face[face_key] = false;
	  ghosted_face.set( face_map[face_key], true);
	  num_ghosts++;
	}
      }
    }
  }

  // Number of own faces
  uint* num_faces= new uint[ pe_size ];
  num_faces[ rank ] = mesh.numFaces() - num_ghosts;
  MPI_Allgather(&num_faces[ rank ], 1, MPI_UNSIGNED, num_faces, 1,
		MPI_UNSIGNED, MPI_COMM_WORLD);
  
  uint num_glb;  
  MPI_Allreduce(&num_faces[ rank ], &num_glb, 1,
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numFaces(num_glb);  
  
  uint offset = 0;
  for(int i = 1; i < rank+1; i++)
    offset += num_faces[i-1];
  
  std::map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numFaces(); i++){
    if( !ghosted_face.get(i) ) {
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;
    }
  }
  
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
      face_key.clear();
      uint num_ok = 0;
      for(uint k = 0; k < inc ; k += 2){
	if(mesh.distdata().have_global(recv_buff[i+k]) &&
	   mesh.distdata().have_global(recv_buff[i+k+1])) {
	  // Generate edge key
	  key = edge_key(mesh.distdata().get_local(recv_buff[i+k], 0),
			 mesh.distdata().get_local(recv_buff[i+k+1], 0));
	  face_key.insert(key);
	  num_ok++;

	}
      }
      if(num_ok < 3)
	continue;
      
      if(owns_face[face_key] && owns_face.count(face_key) == 1) {
	global_buff.push_back(ii);
	global_buff.push_back( new_global[ face_map[face_key] ] );	  	  
      }
    }
    
    MPI_Sendrecv(&global_buff[0], global_buff.size(), MPI_UNSIGNED, src, 2,
		 recv_buff, max_un, MPI_UNSIGNED, dest, 2,
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
    
    for(int i = 0 ; i < recv_count; i += 2){
      new_global[ send_mapping[ recv_buff[i]] ] = recv_buff[i+1];
      new_local[ recv_buff[i+1] ] =  send_mapping[ recv_buff[i] ];
    }
    global_buff.clear();
  }
  
  mesh.distdata().local_face_indices = new_local;
  mesh.distdata().global_face_indices = new_global;
  mesh.distdata()._valid_face_numbering =  true;

  delete[]  num_faces;

}
//-----------------------------------------------------------------------------
void MeshRenumber::renumber_cells(Mesh& mesh)
{
  if( mesh.distdata()._valid_cell_numbering || MPI::numProcesses() == 1)
    return;

  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint* num_cells = new uint[ pe_size ];
  
  // Number of own cells
  num_cells[ rank ] = mesh.numCells();
  MPI_Allgather(&num_cells[ rank ], 1, MPI_UNSIGNED, num_cells, 1,
		MPI_UNSIGNED, MPI_COMM_WORLD);
  
  uint offset = 0;
  for(uint i = 1; i < rank+1; i++)
    offset += num_cells[i-1];

  std::map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numCells(); i++){
    new_global[i] = offset++;
    new_local [ new_global[i] ] = i;
  }


  uint num_glb;  
  MPI_Allreduce(&num_cells[ rank ], &num_glb, 1,
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numCells(num_glb);
  

  delete[]  num_cells;
  


  // Use new numbering
  mesh.distdata().local_cell_indices = new_local;
  mesh.distdata().global_cell_indices = new_global;
  
  mesh.distdata()._valid_cell_numbering = true;
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

