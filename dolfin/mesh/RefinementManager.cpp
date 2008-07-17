// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-01-21
// Last changed: 2008-03-07

#include "RefinementManager.h"
#include <dolfin/main/MPI.h>
#include <dolfin/common/types.h>
#include "MeshFunction.h"
#include "MeshData.h"
#include "BoundaryMesh.h"

#include "Facet.h"

#include <cstdlib> 
#include <ctime> 
#include <map>
#include <set>
#include <mpi.h>


using namespace dolfin;
//-----------------------------------------------------------------------------
RefinementManager::RefinementManager() :_start_offset(0), _refm_init(false)
{

}
//-----------------------------------------------------------------------------
RefinementManager::RefinementManager(Mesh& mesh) :_start_offset(0), 
						  _refm_init(false)
{
  init(mesh);
}
//-----------------------------------------------------------------------------
RefinementManager::~RefinementManager() 
{

}
//-----------------------------------------------------------------------------
void RefinementManager::init(Mesh& mesh)
{
  uint pe_size = MPI::numProcesses();
  if(pe_size == 1)
    return;
  uint rank = MPI::processNumber();


  uint max_index = std::max(mesh.distdata().global_numVertices(),
			    mesh.distdata().max_index());

  // Assume uniform refinement
  uint num_new_vertices = mesh.size(1);

  // Find out maximum global index assigned
  uint glb_max;
  MPI_Allreduce(&max_index, &glb_max, 1, MPI_UNSIGNED, 
		MPI_MAX, MPI_COMM_WORLD);

  // Gather number of new vertices for each processor
  uint *num_new = new uint[pe_size];
  num_new[rank] = num_new_vertices;
  MPI_Allgather(&num_new[rank], 1, MPI_UNSIGNED,
  		num_new, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
  
  // Assign a safe range for each processor
  _start_offset = glb_max;
  for(uint i=0; i < rank; i++)
    _start_offset += num_new[i];

  // Generate cell - edge connectivity if not generated
  mesh.init(mesh.topology().dim(), 1);
  
  // Generate edge - vertex connectivity if not generated
  mesh.init(1, 0);

  // Generate facet - cell connectivity if not generated
  mesh.init(mesh.topology().dim() - 1, mesh.topology().dim());
  
  // Generate facet - edge connectivity if not generated
  mesh.init(mesh.topology().dim() - 1, 1);

  // Initialize datastructures for processor boundary
  BoundaryMesh local_boundary;
  local_boundary.init_local(mesh);
  MeshFunction<uint>* cell_map = local_boundary.data().meshFunction("cell map");

  cell_forbidden.init(mesh, mesh.topology().dim());
  cell_forbidden = false;

  edge_forbidden.init(mesh, 1);
  edge_forbidden = false;

  boundary_edge.init(mesh, 1);
  boundary_edge = false;

  // Check if facets are faces and not edges
  if(mesh.topology().dim() > 1) {
    for(CellIterator bf(local_boundary); !bf.end(); ++bf){
      Facet f(mesh, cell_map->get(*bf));

      for(CellIterator c(f); !c.end(); ++c) {
      	boundary_set.insert(c->index());
	
	for(EdgeIterator e(*c); !e.end(); ++e) {
	  const uint *edge_v = e->entities(0);
	  if(mesh.distdata().is_shared(edge_v[0]) &&
	     mesh.distdata().is_shared(edge_v[1])) {
	    uint edge_vert[2];
	    edge_vert[0] = mesh.distdata().get_global(edge_v[0], 0);
	    edge_vert[1] = mesh.distdata().get_global(edge_v[1], 0);
	    EdgeKey key = edge_key(edge_vert[0], edge_vert[1]);
	    refined_edge[key] = false;
	    edge_cell_map[key] = e->index();          	
	    boundary_edge.set(*e, true);               
	  }
	}
      }
    }
  }
  else {
    for(EdgeIterator be(local_boundary); !be.end(); ++be) {
      Edge e(mesh, cell_map->get(be->index()));
      const uint *edge_v = e.entities(0);
      if(mesh.distdata().is_shared(edge_v[0]) &&
	 mesh.distdata().is_shared(edge_v[1])) {
	uint edge_vert[2];
	edge_vert[0] = mesh.distdata().get_global(edge_v[0], 0);
	edge_vert[1] = mesh.distdata().get_global(edge_v[1], 0);
	EdgeKey key = edge_key(edge_vert[0], edge_vert[1]);
	refined_edge[key] = false;
	edge_cell_map[key] = e.index();          	
	boundary_edge.set(e, true);               
      }
      
      for(CellIterator bc(e); !bc.end(); ++bc) 
	boundary_set.insert(bc->index());
    }
  }

  _refm_init = true;

  delete[] num_new;

}
//-----------------------------------------------------------------------------
void RefinementManager::add_new_vertex(uint* edge,uint vertex, 
				       Mesh& mesh, bool shared)
{
  // Invalidate global numbering
  mesh.distdata().invalid_numbering();
  
  // Invalidate mesh entity ownership
  mesh.distdata().invalid_ownership();

  // Store edge key in shared list
  if(shared) {
    EdgeKey key = edge_key(edge[0], edge[1]);  
    new_global[key] = _start_offset;
    new_vertex[key] = vertex;
  }

  // Assign a new unique global number
  mesh.distdata().set_map(vertex, _start_offset++, 0);  
}
//-----------------------------------------------------------------------------
void RefinementManager::map_new_vertices(Array<uint> shared_edge,
					 Mesh& oldmesh, Mesh& newmesh)
{

  newmesh.distdata().invalid_numbering();
  if(!_refm_init)
    error("RefinementManager not initialized");
  
  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  srand((uint)time(0)); 
  Array<uint> send_buff, send_buff_id;
  std::map<EdgeKey, uint> edge_id;  
  std::map<EdgeKey, bool> owns_edge;  
  EdgeKey  key;
  for(uint i = 0; i < shared_edge.size(); i +=3) {

    key = edge_key(shared_edge[i], shared_edge[i+1]);
    edge_id[key] = (uint)rand() + rank;
    owns_edge[key] = true;
        
    send_buff.push_back(oldmesh.distdata().get_global(shared_edge[i], 0));
    send_buff.push_back(oldmesh.distdata().get_global(shared_edge[i+1], 0));
    send_buff_id.push_back(edge_id[key]);
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
  for(int j = 1 ; j < pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    
    
    MPI_Sendrecv(&send_buff_id[0], num_id, MPI_UNSIGNED, dest, 1, 
		 recv_buff_id, max_id, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
    
    for(int i =0; i < recv_count ; i += 2){
      // Check if I have the vertices
      if(oldmesh.distdata().have_global(recv_buff[i]) &&
	 oldmesh.distdata().have_global(recv_buff[i+1])) {

	// Generate edge key
	key = edge_key(oldmesh.distdata().get_local(recv_buff[i], 0),
		       oldmesh.distdata().get_local(recv_buff[i+1], 0));

	// Check if I have the corresponding edge
	if(edge_id.count(key)) {
	  if( recv_buff_id[i>>1] < edge_id[key] || 
	      recv_buff_id[i>>1] == edge_id[key] && status.MPI_SOURCE < rank){
	    owns_edge[key] = false;
	    new_global.erase(key);	    
	    edge_id.erase(key);
	  }
	}	
      }
    }
  }

  //Exchange assigned global numbers
  Array<uint> global_buff;
  uint index;  
  for(int j = 1; j < pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    
    
    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
    
    for(int i =0; i < recv_count ; i += 2){      
      if(oldmesh.distdata().have_global(recv_buff[i]) &&
	 oldmesh.distdata().have_global(recv_buff[i+1])) {
	
	key = edge_key(oldmesh.distdata().get_local(recv_buff[i], 0),
		       oldmesh.distdata().get_local(recv_buff[i+1], 0));

	if(owns_edge[key] && owns_edge.count(key) == 1) {
	  global_buff.push_back(i);
	  global_buff.push_back( new_global[key] ); 
	  newmesh.distdata().set_shared(new_vertex[key]);
	}
	else if(owns_edge.count(key) > 1)
	  error("this shouldnt happen");
      }
    }
    
    MPI_Sendrecv(&global_buff[0], global_buff.size(), MPI_UNSIGNED, src, 2,
		 recv_buff, max_un, MPI_UNSIGNED, dest, 2,
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
    
    for(int i = 0 ; i < recv_count; i += 2){
      index = shared_edge[(recv_buff[i]>>1) * 3 + 2];
      dolfin_assert(!newmesh.distdata().have_global(recv_buff[i+1]));
      newmesh.distdata().set_map(index, recv_buff[i+1], 0);
      newmesh.distdata().set_ghost(index);
      newmesh.distdata().set_ghost_owner(index ,status.MPI_SOURCE);
    }
    global_buff.clear();
  }
  
  // MPI aliasing 
  uint tmp =  newmesh.numVertices() - newmesh.distdata().num_ghost();
  uint num_glb;  
  MPI_Allreduce(&tmp, &num_glb, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  newmesh.distdata().set_global_numVertices(num_glb);
  delete[] recv_buff;
  delete[] recv_buff_id;
}
//-----------------------------------------------------------------------------
void RefinementManager::mark_localboundary(Mesh& oldmesh,
					   MeshFunction<bool>& cell_marker,
					   uint& num_new_vertices,
					   uint& num_propagated, 
					   uint& num_new_cells)
{
  
  if(!_refm_init)
    error("RefinementManager not initialized");

  srand((uint)time(0));   
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();

  Array<uint> send_buff;
  uint edge[2];
  EdgeKey key;
  // Set of cell indices with longest edge
  std::map<EdgeKey, uint> num_ref, removed, edge_id;
  std::map<uint, EdgeKey> cell_keymap;
  std::set<uint> cell_forbidden_edges;
  std::set<EdgeKey> forbidden_propagation;
  std::map<uint, uint> edge_vote;

  num_propagated = 0;
  num_new_cells = 0;
  num_new_vertices = 0;
  real max, l;
  uint index = 0;

  // Process cells between processors
  std::set<uint>::iterator bc;

  // Reset forbidden edges and cells
  edge_forbidden = false;
  cell_forbidden = false;

  for(bc = boundary_set.begin(); bc != boundary_set.end(); bc++){ 
    Cell c(oldmesh, *bc);    
    if(cell_marker.get(c) && !cell_forbidden.get(c)) {
      max = 0.0;
      for(EdgeIterator e(c); !e.end(); ++e){
	if(edge_forbidden.get(*e))
	  continue;
	edge_vote[e->index()] = 0;	
	l = e->length();
	if( max < l ) {
	  max = l ;
	  index = e->index();
	}	
      }
      
      if( max > 0.0 ) {
	Edge longest_edge(oldmesh, index);	
	if( on_boundary(longest_edge) ){
	  const uint *edge_v = longest_edge.entities(0);
	  edge[0] = oldmesh.distdata().get_global(edge_v[0], 0);
	  edge[1] = oldmesh.distdata().get_global(edge_v[1], 0);
	  key = edge_key(edge[0], edge[1]);
	  edge_vote[longest_edge.index()] += (uint) rand() + rank; 

	  send_buff.push_back(edge[0]);
	  send_buff.push_back(edge[1]);
	  send_buff.push_back(edge_vote[longest_edge.index()]);


	  for(CellIterator nc(longest_edge); !nc.end(); ++nc){	 
	    dolfin_assert(!cell_forbidden.get(*nc));
	    cell_forbidden.set(*nc, true);
	    for(EdgeIterator e(*nc); !e.end(); ++e){
	      edge_forbidden.set(*e, true);	    
	    }
	  }
	}
      }
    }
  }

  // Decide ownership of refienment
  MPI_Status status;
  uint src,dest;
  int max_un, num_un, recv_count;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);
  uint *recv_buff = new uint[max_un];
  Array<uint> forbidden;
  
  for(uint j=1; j<pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    
    
    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);  
    
    for(int i = 0; i < recv_count; i +=3) {
      key = edge_key(recv_buff[i], recv_buff[i+1]);
      // If rank has edge
      if(edge_cell_map.count(key)) {
	Edge e(oldmesh, edge_cell_map[key]);
	edge_vote[e.index()] +=recv_buff[i+2];
      }
    }
  }

  send_buff.clear();
  cell_forbidden = false;
  edge_forbidden = false;

  for(bc = boundary_set.begin(); bc != boundary_set.end(); bc++){ 
    Cell c(oldmesh, *bc);    
    if(cell_marker.get(c) && !cell_forbidden.get(c)) {
      max = 0.0;
      for(EdgeIterator e(c); !e.end(); ++e){
	if(edge_forbidden.get(*e))
	  continue;
	l = edge_vote[e->index()];
	if( max < l ) {
	  max = l ;
	  index = e->index();
	}
      }
      
      if( max > 0.0 ) {
	Edge longest_edge(oldmesh, index);	
	if( on_boundary(longest_edge) ){
	  const uint *edge_v = longest_edge.entities(0);
	  edge[0] = oldmesh.distdata().get_global(edge_v[0], 0);
	  edge[1] = oldmesh.distdata().get_global(edge_v[1], 0);
	  key = edge_key(edge[0], edge[1]);
	  refined_edge[key] = true;
	  num_new_vertices++;	  	  
	  send_buff.push_back(edge[0]);
	  send_buff.push_back(edge[1]);

	  for(CellIterator nc(longest_edge); !nc.end(); ++nc){	 
	    dolfin_assert(!cell_forbidden.get(*nc));
	    cell_forbidden.set(*nc, true);
	    num_ref[key]++;
            num_new_cells++;
	    cell_refedge[nc->index()] =  longest_edge.index();
            cell_keymap[nc->index()] = key;
	    for(EdgeIterator e(*nc); !e.end(); ++e){
	      edge_forbidden.set(*e, true);	    
	      cell_forbidden_edges.insert(e->index());
	    }
	  }
	}
      }
    }
  }

  Array<uint> terminated;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);
  delete[] recv_buff;
  recv_buff = new uint[max_un];

  for(uint j=1; j<pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);  

    for(int i = 0; i < recv_count; i +=2) {
      key = edge_key(recv_buff[i], recv_buff[i+1]);
      // If rank has edge
      if(edge_cell_map.count(key)) {
	if(!refined_edge[key]) {
	  Edge e(oldmesh, edge_cell_map[key]);
	  if(cell_forbidden_edges.count(e.index()) == 0) {
	    num_new_vertices++;
	    refined_edge[key] = true;
	    for(CellIterator c(e); !c.end(); ++c){
	      dolfin_assert(!cell_forbidden.get(*c));
	      num_new_cells++;
	      cell_forbidden.set(*c, true);
	      cell_refedge[c->index()] = e.index();
	      for(EdgeIterator ce(*c); !ce.end(); ++ce){
		edge_forbidden.set(*ce, true);
		cell_forbidden_edges.insert(ce->index());
	      }
	    }	  
	  }
	  else {
	    forbidden.push_back(recv_buff[i]);
	    forbidden.push_back(recv_buff[i+1]);
	  }	  
	}
      }
    }
    MPI_Sendrecv(&forbidden[0], forbidden.size(), MPI_UNSIGNED, src, 2, 
		 recv_buff, max_un, MPI_UNSIGNED, dest, 2, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);  
    
    for(int i = 0; i < recv_count; i +=2) {
      key = edge_key(recv_buff[i], recv_buff[i+1]);
      if(removed.count(key) == 0) {
	num_ref[key]--;
	if(num_ref[key] == 0)
	  num_new_vertices--;      
	Edge e(oldmesh, edge_cell_map[key]);
	for(CellIterator c(e); !c.end(); ++c){
	  if(cell_forbidden.get(*c))
	    num_new_cells--;
	  cell_forbidden.set(*c, false);
	  cell_refedge.erase(c->index());
	}
	terminated.push_back(recv_buff[i]);
	terminated.push_back(recv_buff[i+1]);
      }
      else
	removed[key] = 1;
    }        
    forbidden.clear();
  }

  num_un = terminated.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);
  delete[] recv_buff;
  recv_buff = new uint[max_un];
  
  for(uint j=1; j<pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;    

    MPI_Sendrecv(&terminated[0], terminated.size(), MPI_UNSIGNED, dest, 1, 
		 recv_buff, max_un, MPI_UNSIGNED, src, 1, 
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);  

    for(int i = 0; i < recv_count; i +=2) {
      key = edge_key(recv_buff[i], recv_buff[i+1]);
      // If rank has edge
      if(edge_cell_map.count(key)) {
	// Remove refinement if edge is refined with terminated refinment
	if(refined_edge[key]) {
	  if(removed.count(key) == 0) {
	    num_ref[key]--;
	    if(num_ref[key] == 0)
	      num_new_vertices--;      
	    Edge e(oldmesh, edge_cell_map[key]);
	    for(CellIterator c(e); !c.end(); ++c){
	      if(cell_forbidden.get(*c))
		num_new_cells--;
	      cell_forbidden.set(*c, false);
	      cell_refedge.erase(c->index());
	    }
	    refined_edge.erase(key);
	  }
	  else
	    removed[key] = 1;
	}
      }
    }
  }


  // Mark unrefined cells shared edges as forbidden
  for(bc = boundary_set.begin(); bc != boundary_set.end(); bc++) {
    Cell c(oldmesh, *bc);
    if(!cell_forbidden.get(c)) {	
      cell_refedge.erase(c.index());
      for(EdgeIterator e(c); !e.end(); ++e) {
       	if(on_boundary(*e)) 
	  edge_forbidden.set(*e, true);
	else {
	  bool ok_to_remove = true;
	  for(CellIterator ec(*e); !ec.end(); ++ec) {
	    if(cell_forbidden.get(*ec)) {
	      ok_to_remove = false;
	      break;
	    }
	  }
	  if(ok_to_remove)
	    edge_forbidden.set(*e, false);
	}
      }
    }	
  }

  delete[] recv_buff;

}
//-----------------------------------------------------------------------------
std::pair<dolfin::uint, dolfin::uint> RefinementManager::edge_key(uint id1,
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
