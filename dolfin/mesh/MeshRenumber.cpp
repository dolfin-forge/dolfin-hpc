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
#include "OwnershipComputation.h"

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
  mesh.init(0, 0);
  
  renumber_vertices(mesh);

  OwnershipComputation::generate_ownership(mesh);

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

  // Number of own vertices
  num_vert[ rank ] = mesh.numVertices() - mesh.distdata().num_ghost();
  MPI_Allgather(&num_vert[ rank ], 1, MPI_UNSIGNED,
		num_vert, 1, MPI_UNSIGNED, MPI_COMM_WORLD);  

  uint offset = 0;
  for(uint i = 1; i < rank+1; i++)
    offset += num_vert[i-1];

  std::map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numVertices(); i++){
    if(!mesh.distdata().is_ghost(i)){
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;
    }
  }

  uint num_glb;  
  MPI_Allreduce(&num_vert[ rank ], &num_glb, 1,
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numVertices(num_glb);

  Array<uint> *ghost_buff = new Array<uint>[pe_size];
  for(MeshGhostIterator iter(mesh.distdata(), 0); !iter.end(); ++iter)
    ghost_buff[iter.owner()].push_back(mesh.distdata().get_global(iter.index(), 0)); 
 
  MPI_Status status;
  Array<uint> send_buff;
  uint src,dest;
  uint recv_size = mesh.distdata().num_ghost(); 
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

  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint* num_edge = new uint[ pe_size ];

  // Number of own vertices
  num_edge[ rank ] = mesh.numEdges() - mesh.distdata().num_ghost(1);
  MPI_Allgather(&num_edge[ rank ], 1, MPI_UNSIGNED,
		num_edge, 1, MPI_UNSIGNED, MPI_COMM_WORLD);  

  uint offset = 0;
  for(uint i = 1; i < rank+1; i++)
    offset += num_edge[i-1];

  std::map<EdgeKey, uint> edge_map;
  std::map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numEdges(); i++){
    if(!mesh.distdata().is_ghost(i, 1)){
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;            
    }
   
    Edge e(mesh, i);
    const uint *edge_v = e.entities(0);      
    EdgeKey key = edge_key(edge_v[0], edge_v[1]);
    edge_map[key] = i;

  }

  uint num_glb;  
  MPI_Allreduce(&num_edge[ rank ], &num_glb, 1,
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numEdges(num_glb);

  Array<uint> *ghost_buff = new Array<uint>[pe_size];
  for(MeshGhostIterator iter(mesh.distdata(), 1); !iter.end(); ++iter) {
    Edge eg(mesh, iter.index());
    const uint *edge_v = eg.entities(0);
    ghost_buff[iter.owner()].push_back(mesh.distdata().get_global(edge_v[0], 0)); 
    ghost_buff[iter.owner()].push_back(mesh.distdata().get_global(edge_v[1], 0)); 
  }
 
  MPI_Status status;
  Array<uint> send_buff;
  uint src,dest;
  uint recv_size = mesh.distdata().num_ghost(1); 
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

    for(int k=0; k < recv_count; k += 2) {
      EdgeKey key = edge_key(mesh.distdata().get_local(recv_ghost[k], 0),
			     mesh.distdata().get_local(recv_ghost[k+1], 0));
      send_buff.push_back(new_global[edge_map[key]]);
    }

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, src, 2,
		 recv_buff, recv_size , MPI_UNSIGNED, dest, 2, 
		 MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);
    
    uint k = 0;
    for(int j=0; j < recv_count; j++){
      EdgeKey key = edge_key(mesh.distdata().get_local(ghost_buff[dest][k] , 0),
			     mesh.distdata().get_local(ghost_buff[dest][k+1] , 0));
      
      new_global[ edge_map[key] ] = recv_buff[j];
      new_local[recv_buff[j] ] = edge_map[key];
      k += 2;
    }    
    send_buff.clear();
  }
  
  // Use new numbering
  mesh.distdata().local_indices[1] = new_local;
  mesh.distdata().global_indices[1] = new_global;  
  mesh.distdata()._valid_edge_numbering = true;
  
  delete[] recv_buff;
  delete[] recv_ghost;
  for(uint i = 0; i < pe_size; i++)
    ghost_buff[i].clear();
  delete[] ghost_buff;
  delete[] num_edge;    
}
//-----------------------------------------------------------------------------
void MeshRenumber::renumber_faces(Mesh& mesh)
{
  if( mesh.distdata()._valid_face_numbering || 
      MPI::numProcesses() == 1 || mesh.topology().dim() == 2)
    return;  

  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint* num_faces = new uint[ pe_size ];
  
  // Number of own vertices
  num_faces[ rank ] = mesh.numFaces() - mesh.distdata().num_ghost(2);
  MPI_Allgather(&num_faces[ rank ], 1, MPI_UNSIGNED,
		num_faces, 1, MPI_UNSIGNED, MPI_COMM_WORLD);  

  uint offset = 0;
  for(uint i = 1; i < rank+1; i++)
    offset += num_faces[i-1];

  std::map<FaceKey, uint> face_map;
  FaceKey face_key;
  std::map<uint,uint> new_local,new_global;  
  for(uint i = 0; i< mesh.numFaces(); i++){
    if(!mesh.distdata().is_ghost(i, 2)){
      new_global[i] = offset++;
      new_local [ new_global[i] ] = i;
    }
    Face f(mesh, i);
    face_key.clear();
    for(EdgeIterator e(f); !e.end(); ++e){
      const uint *edge_v = e->entities(0);      
      face_key.insert(edge_key(edge_v[0], edge_v[1]));
    }
    face_map[face_key] = f.index();
  }

  uint num_glb;  
  MPI_Allreduce(&num_faces[ rank ], &num_glb, 1,
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  mesh.distdata().set_global_numFaces(num_glb);

  Array<uint> *ghost_buff = new Array<uint>[pe_size];
  for(MeshGhostIterator iter(mesh.distdata(), 2); !iter.end(); ++iter) {
    Face f(mesh, iter.index());
    for(EdgeIterator e(f); !e.end(); ++e)
      ghost_buff[iter.owner()].push_back(mesh.distdata().get_global(e->index(), 1)); 
  }
 
  MPI_Status status;
  Array<uint> send_buff;
  uint src,dest;
  uint recv_size = mesh.distdata().num_ghost(2); 
  int recv_count, recv_size_gh, send_size;  
  Face f(mesh, 0);  
  uint num_edges =  f.numEntities(1);
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
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for(int k=0; k < recv_count; k +=num_edges) {
      face_key.clear();
      for(uint i = 0; i < num_edges; i++) {
	Edge e(mesh, mesh.distdata().get_local(recv_ghost[k+i], 1));
	const uint *edge_v = e.entities(0);      
	face_key.insert(edge_key(edge_v[0], edge_v[1]));      
      }
      send_buff.push_back(new_global[ face_map[face_key] ]);
    }

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, src, 2,
		 recv_buff, recv_size , MPI_UNSIGNED, dest, 2, 
		 MPI_COMM_WORLD,&status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
    
    uint k = 0;
    for(int j=0; j < recv_count; j++){
      face_key.clear();
      for(uint i = 0; i < num_edges; i++) {
	Edge e(mesh, mesh.distdata().get_local(ghost_buff[dest][i + k], 1));
	const uint *edge_v = e.entities(0);      
	face_key.insert(edge_key(edge_v[0], edge_v[1]));      
      }
      new_global[ face_map[face_key] ] = recv_buff[j];
      new_local[ recv_buff[j] ] =  face_map[face_key];
      k += num_edges;
    }    
    send_buff.clear();
  }
  
  // Use new numbering
  mesh.distdata().local_indices[2] = new_local;
  mesh.distdata().global_indices[2] = new_global;  
  mesh.distdata()._valid_face_numbering = true;
  
  delete[] recv_buff;
  delete[] recv_ghost;
  for(uint i = 0; i < pe_size; i++)
    ghost_buff[i].clear();
  delete[] ghost_buff;
  delete[] num_faces;    
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
  mesh.distdata().local_indices[3] = new_local;
  mesh.distdata().global_indices[3] = new_global;
  
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

