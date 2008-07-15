// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include "OwnershipComputation.h"
#include "Mesh.h"
#include "MeshData.h"
#include "MeshRenumber.h"
#include "BoundaryMesh.h"
#include "Edge.h"
#include "Cell.h"
#include "Face.h"
#include "Vertex.h"

#include "GlobalFacetMap.h"

#include <dolfin/common/Array.h>
#include <dolfin/main/MPI.h>

#include <mpi.h>
#include <map>
#include <cstdlib>
#include <ctime>

using namespace dolfin;

//-----------------------------------------------------------------------------
void OwnershipComputation::generate_ownership(Mesh& mesh)
{

  BoundaryMesh local_boundary;
  local_boundary.init_local(mesh);

  MeshRenumber::renumber_vertices(mesh);
  init_edge_ownership(mesh, local_boundary);
  init_face_ownership(mesh, local_boundary);
}
//-----------------------------------------------------------------------------
void OwnershipComputation::init_edge_ownership(Mesh& mesh,
					       BoundaryMesh& local_boundary)
{
  if( mesh.distdata()._valid_edge_ownership || MPI::numProcesses() == 1)
    return;  

  mesh.distdata().ghost[1].clear();
  mesh.distdata().shared[1].clear();

  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  std::map<EdgeKey, uint> edge_map, edge_id;  

  Array<uint> send_buff, send_buff_id;
  std::map<uint,uint> send_mapping;
  std::set<uint> used_edge;
  MeshFunction<uint>* cell_map = local_boundary.data().meshFunction("cell map");

  srand((uint)time(0));
  for(CellIterator bf(local_boundary); !bf.end(); ++bf){
    Facet f(mesh, cell_map->get(*bf));
    for(EdgeIterator e(f); !e.end(); ++e) {
      if(used_edge.count(e->index()) == 0) {
	const uint *edge_v = e->entities(0);
	if(mesh.distdata().is_shared( edge_v[0], 0) &&
	   mesh.distdata().is_shared( edge_v[1], 0)) {
	  EdgeKey key = edge_key(edge_v[0], edge_v[1]);
	  edge_map[key] = e->index();
	  edge_id[key] = (uint) rand() + rank;
	  send_buff.push_back(mesh.distdata().get_global(edge_v[0], 0));
	  send_buff.push_back(mesh.distdata().get_global(edge_v[1], 0));
	  send_buff_id.push_back(edge_id[key]);
	  used_edge.insert(e->index());     
	  mesh.distdata().set_shared(*e);
	}       	
      }
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
  std::map<uint, uint> edge_votes;

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
	      recv_buff_id[i>>1] == edge_id[key] && status.MPI_SOURCE < rank ){
	    edge_id.erase(key);
	    edge_votes[ edge_map[key] ] = recv_buff_id[i>>1];
	    mesh.distdata().set_ghost(edge_map[key], 1);
	    mesh.distdata().set_ghost_owner(edge_map[key], status.MPI_SOURCE, 1);
	  }
	}
	else {
	  if(edge_map.count(key)){ 
	    if( recv_buff_id[i>>1] < edge_votes[ edge_map[key] ]) {		
	      edge_votes[ edge_map[key] ]  = recv_buff_id[i>>1];
	      mesh.distdata().set_ghost_owner(edge_map[key], status.MPI_SOURCE, 1);
	    }
	  }
	}
      }
    }
  }

  mesh.distdata()._valid_edge_ownership = true;

}
//-----------------------------------------------------------------------------
void OwnershipComputation::init_face_ownership(Mesh& mesh, 
					       BoundaryMesh& local_boundary)
{
  if( mesh.distdata()._valid_face_ownership ||
      MPI::numProcesses() == 1 || mesh.topology().dim() == 2)
    return;  

  mesh.distdata().ghost[2].clear();
  mesh.distdata().shared[2].clear();

  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  std::map<FaceKey, uint> face_map, face_id;
  std::map<uint, uint> face_votes;
  std::set<uint> used_face;

  Array<uint> send_buff, send_buff_id;
  FaceKey face_key;

  GlobalFacetMap facetmap(mesh);
  facetmap.init();

  MeshFunction<uint>* cell_map = local_boundary.data().meshFunction("cell map");

  srand((uint)time(0));
  for(CellIterator bf(local_boundary); !bf.end(); ++bf){
    Face f(mesh, cell_map->get(*bf));
    Facet ff(mesh, cell_map->get(*bf));    
    if(!facetmap.globalFacet(ff) ) {
      face_key.clear();
      for(EdgeIterator e(f); !e.end(); ++e) {
	const uint *edge_v = e->entities(0);
	face_key.insert(edge_key(edge_v[0], edge_v[1]));
	send_buff.push_back(mesh.distdata().get_global(edge_v[0], 0));
	send_buff.push_back(mesh.distdata().get_global(edge_v[1], 0));	
      }
      if(face_map.count(face_key) == 0) {
	face_map[face_key] = f.index();
	face_id[face_key] = (uint) rand() +  rank;
	send_buff_id.push_back(face_id[face_key]);
	mesh.distdata().set_shared(f.index(), 2);
      }
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
      face_key.clear();
      uint num_ok = 0;
      for(uint k = 0; k < inc; k += 2){
	// Check if I have the vertices     
	if(mesh.distdata().have_global(recv_buff[i+k]) &&
	   mesh.distdata().have_global(recv_buff[i+k+1])) {
	  // Generate edge key
	  face_key.insert(edge_key(mesh.distdata().get_local(recv_buff[i+k], 0),
				   mesh.distdata().get_local(recv_buff[i+k+1], 0)));
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
	  mesh.distdata().set_ghost(face_map[face_key], 2);
	  mesh.distdata().set_ghost_owner(face_map[face_key], 
					  status.MPI_SOURCE, 2);
	  face_votes[ face_map[face_key] ]  = recv_buff_id[ii];
	}
      }
      else {
	if( face_map.count(face_key) )
	  if( recv_buff_id[ii] < face_votes[ face_map[face_key] ]) {		
	    face_votes[ face_map[face_key] ]  = recv_buff_id[ii];
	    mesh.distdata().set_ghost_owner(face_map[face_key],
					    status.MPI_SOURCE, 2);
	    
	  }
      }
    }
  }

  mesh.distdata()._valid_face_ownership = true;
}
//-----------------------------------------------------------------------------
std::pair<dolfin::uint, dolfin::uint> OwnershipComputation::edge_key(uint id1,
								     uint id2)
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
