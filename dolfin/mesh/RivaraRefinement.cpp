// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
//


#include <dolfin/main/MPI.h>
#include <dolfin/common/constants.h>
#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include "Mesh.h"
#include "MeshData.h"
#include "MeshTopology.h"
#include "MeshGeometry.h"
#include "MeshConnectivity.h"
#include "MeshEditor.h"
#include "MeshFunction.h"
#include "MeshRenumber.h"
#include "Vertex.h"
#include "Facet.h"
#include "Edge.h"
#include "Cell.h"
#include "BoundaryMesh.h"
#include "LoadBalancer.h"
#include "RivaraRefinement.h"

#include <algorithm>
#include <cmath> 

#ifdef HAS_MPI
#include <mpi.h>
#endif 

using namespace dolfin;
//-----------------------------------------------------------------------------
void RivaraRefinement::refine(Mesh& mesh, 
			      MeshFunction<bool>& cell_marker,
			      real tf, real tb, real ts)
{
  message("Refining simplicial mesh by recursive Rivara bisection.");

  // Start Loadbalancer
  if(MPI::numProcesses() > 1) {
    begin("Load balancing");
    // Tune loadbalancer using machine specific parameters, if available
    if( tf > 0.0 && tb > 0.0 && ts > 0.0)
      LoadBalancer::balance(mesh, cell_marker, tf, tb, ts, LoadBalancer::LEPP);
    else
      LoadBalancer::balance(mesh, cell_marker, LoadBalancer::LEPP);
    end();
  }

  DMesh dmesh;
  dmesh.imp(mesh);
  
  std::vector<bool> dmarked(mesh.numCells());
  for (CellIterator ci(mesh); !ci.end(); ++ci)
  {
    if(cell_marker.get(*ci) == true)
    {
      dmarked[ci->index()] = true;
    }  
    else
    {
      dmarked[ci->index()] = false;
    }
  }
  
  dmesh.bisectMarked(dmarked);

  // Remove deleted cells from global list
  for(std::list<DCell* >::iterator it = dmesh.cells.begin();
      it != dmesh.cells.end(); )
  {
    
    DCell* dc = *it;
    
    if(dc->deleted)
      it = dmesh.cells.erase(it);
    else
      it++;
  } 
  
  Mesh omesh;    
  dmesh.exp(omesh);
  
  mesh = omesh;
  mesh.distdata().invalid_numbering();
  mesh.distdata().invalid_ownership();
  mesh.renumber();
}
//-----------------------------------------------------------------------------
DVertex::DVertex() : id(0), glb_id(-1), cells(0), p(0.0, 0.0, 0.0), 
		     on_boundary(false), shared(false),
		     ghosted(false), owner(-1)
{
}
//-----------------------------------------------------------------------------
DCell::DCell() : id(0), parent_id(0), vertices(0), deleted(false), nref(0)
{
}
//-----------------------------------------------------------------------------
bool DCell::has_edge(DVertex *v1, DVertex *v2)
{
  uint found = 0;
  for(std::vector<DVertex*>::iterator it = vertices.begin();
      it != vertices.end(); ++it)
    if(*it == v1 || *it == v2)
      found++;
  return (found == 2);
}
//-----------------------------------------------------------------------------
//DMesh::DMesh() : vertices(0), cells(0)
DMesh::DMesh() :  cells(0)
{
}
//-----------------------------------------------------------------------------
DMesh::~DMesh()
{
  // Delete allocated DVertices
  for(std::set<DVertex* >::iterator it = vertices.begin();
      it != vertices.end(); ++it)
    delete *it;
  
  // Delete allocated DCells
  for(std::list<DCell* >::iterator it = cells.begin();
       it != cells.end(); ++it)
     delete *it;
}
//-----------------------------------------------------------------------------
void DMesh::imp(Mesh& mesh)
{
  cell_type = &(mesh.type());
  d = mesh.topology().dim();

  vertices.clear();
  cells.clear();

  MeshRenumber::renumber_vertices(mesh);

  std::vector<DVertex *> vertexvec;

  BoundaryMesh boundary;
  boundary.init_interior(mesh);
  MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");  
  MeshFunction<bool> boundary_cell(mesh, mesh.topology().dim());
  boundary_cell = false;

  // Generate facet - cell connectivity if not generated
  mesh.init(mesh.topology().dim() - 1, mesh.topology().dim());
  for (CellIterator bf(boundary); !bf.end(); ++bf) 
  {
    Facet f(mesh, cell_map->get(*bf));    
    for (CellIterator c(f); !c.end(); ++c) 
      boundary_cell.set(*c, true);	  
  }

  // Assume uniform refinement
  uint num_new = mesh.size(1);

  // Since the mesh is linear numbered, the maximum global index assigned is
  // the number of vertices in the mesh
  uint max_index = mesh.distdata().global_numVertices();  
#ifdef HAS_MPI
  MPI_Allreduce(&max_index, &glb_max, 1, MPI_UNSIGNED, MPI_MAX, MPI::DOLFIN_COMM);
#endif

  Cell c(mesh, 0);
  _salt = c.numEntities(1) * mesh.distdata().global_numCells();

  // Assign a safe range for each processor
  _start_offset = 0;
#ifdef HAS_MPI
#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_new, &_start_offset, 1,
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
  MPI_Scan(&num_new, &_start_offset, 1,
	   MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
  _start_offset -= num_new;
#endif
  _start_offset += glb_max;

  MPI_Allreduce(&_start_offset, &_max, 1, MPI_UNSIGNED, MPI_MAX, MPI::DOLFIN_COMM);
  #endif

  uint counter = 1;
  for (VertexIterator vi(mesh); !vi.end(); ++vi)
  {
    DVertex* dv = new DVertex;    
    dv->p = vi->point();
    dv->glb_id = mesh.distdata().get_global(vi->index(), 0);
    dv->on_boundary = mesh.distdata().is_shared(vi->index(), 0);
    dv->shared = mesh.distdata().is_shared(vi->index(), 0);
    dv->ghosted = mesh.distdata().is_ghost(vi->index(), 0);
    if (dv->ghosted)
      dv->owner = mesh.distdata().get_owner(*vi);    

    if(dv->on_boundary)     
      bc_dvs[dv->glb_id] = dv;
    
    vertices.insert(dv);
    vertexvec.push_back(dv);
    counter++;
  }

  for (CellIterator ci(mesh); !ci.end(); ++ci)
  {
    DCell* dc = new DCell;

    std::vector<DVertex*> vs(ci->numEntities(0));
    uint i = 0;
    for (VertexIterator vi(*ci); !vi.end(); ++vi)
    {
      DVertex* dv = vertexvec[vi->index()];

      vs[i] = dv;
      i++;
    }

    addCell(dc, vs, ci->index());
    // Define the same cell numbering
    dc->id = ci->index();
    
    // Add dynamic cell to list of boundary cells
    if ( boundary_cell.get(*ci) ) 
    {
      for (EdgeIterator e(*ci); !e.end(); ++e) 
      {
	const uint *edge_v = e->entities(0);
	if( mesh.distdata().is_shared(edge_v[0], 0) || mesh.distdata().is_shared(edge_v[1], 0)) 
	{
	  EdgeKey key = edge_key(mesh.distdata().get_global(edge_v[0], 0),
				 mesh.distdata().get_global(edge_v[1], 0));
	}
      }
    }
  }
}
//-----------------------------------------------------------------------------
void DMesh::exp(Mesh& mesh)
{
  number();

  MeshEditor editor;
  editor.open(mesh, cell_type->cellType(), d, d);
  
  editor.initVertices(vertices.size());
  editor.initCells(cells.size());

  // Add old vertices
  uint current_vertex = 0;
  for(std::set<DVertex* >::iterator it = vertices.begin();
      it != vertices.end(); ++it)
  {
    DVertex* dv = *it;
    editor.addVertex(current_vertex, dv->p);
    if(dv->ghosted) {
      mesh.distdata().set_ghost(current_vertex, 0);
      mesh.distdata().set_ghost_owner(current_vertex, dv->owner, 0);
    }
    else if(dv->shared)
      mesh.distdata().set_shared(current_vertex, 0);
    mesh.distdata().set_map(current_vertex++, dv->glb_id, 0);      
  }

  Array<uint> cell_vertices(cell_type->numEntities(0));
  uint current_cell = 0;
  for(std::list<DCell* >::iterator it = cells.begin();
      it != cells.end(); ++it)
  {
    DCell* dc = *it;

    for(uint j = 0; j < dc->vertices.size(); j++)
    {
      DVertex* dv = dc->vertices[j];
      cell_vertices[j] = dv->id;
    }
    editor.addCell(current_cell, cell_vertices);

    current_cell++;
  }
  editor.close();

}
//-----------------------------------------------------------------------------
void DMesh::number()
{
  uint i = 0;
  for(std::set<DVertex* >::iterator it = vertices.begin();
      it != vertices.end(); ++it)
  {
    DVertex* dv = *it;
    dv->id = i;
    i++;
  }

  i = 0;
  for(std::list<DCell* >::iterator it = cells.begin();
      it != cells.end(); ++it)
  {
    DCell* dc = *it;
    dc->id = i;
    i++;
  }  
}
//-----------------------------------------------------------------------------
void DMesh::bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1)
{


  bool closing = false;

  // Find longest edge
  real lmax = 0.0;
  int ptmax = 0;
  uint ii = 0;
  uint jj = 0;
  for(uint i = 0; i < dcell->vertices.size(); i++)
  {
    for(uint j = 0; j < dcell->vertices.size(); j++)
    {
      if(i != j)
      {
	DVertex* v0 = dcell->vertices[i];
	DVertex* v1 = dcell->vertices[j];
	
	real l = 0.0;
	if( v0->glb_id > v1->glb_id)
	  l = v0->p.distance(v1->p); 
	else
	  l = v1->p.distance(v0->p); 
	
	if(fabs(l - lmax) < DOLFIN_EPS)
	{
	  int ptsum = (v0->glb_id) + (v1->glb_id) ;
	  if(ptsum > ptmax)
	  {
	    ii = i;
	    jj = j;
	    lmax = l;
	    ptmax = (v0->glb_id + v1->glb_id);
	  }
	}
	else if(l >= lmax)
 	{
 	  ii = i;
 	  jj = j;
 	  lmax = l;
	  ptmax = (v0->glb_id + v1->glb_id) ;
	}
      }
    }
  }
  
  DVertex* v0 = dcell->vertices[ii];
  DVertex* v1 = dcell->vertices[jj];
  DVertex* mv = 0;

  // Check if no hanging vertices remain, otherwise create hanging
  // vertex and continue refinement
  if((v0 == hv0 || v0 == hv1) && (v1 == hv0 || v1 == hv1))
  {

    mv = hangv;
    closing = true;

    if( v0->on_boundary && v1->on_boundary ) 
    {
      mv->on_boundary = true;
      mv->shared = true;
      bc_dvs[mv->glb_id] = mv;		      
      dolfin_assert(v0->glb_id != v1->glb_id);
      dolfin_assert( ref_edge.find(edge_key(v0->glb_id, v1->glb_id)) != ref_edge.end());
    }
  }
  else
  {
    mv = new DVertex;
    addVertex(mv);
    if( v0->glb_id < v1->glb_id )
      mv->glb_id = (((v0->glb_id * _salt) + (v1->glb_id))) + glb_max;
    else
      mv->glb_id = (((v1->glb_id * _salt) + (v0->glb_id))) + glb_max;
    mv->p = (dcell->vertices[ii]->p + dcell->vertices[jj]->p) / 2.0; 

    // Add hanging node on shared edges to propagation buffer
    if( v0->on_boundary && v1->on_boundary) 
    {
      prop_edge node;
      node.mv = mv->glb_id;
      node.v1 = v0->glb_id;
      node.v2 = v1->glb_id;
      node.owner = MPI::processNumber();
      std::pair<uint, prop_edge> _prop_( dcell->nref, node);
      propagate.push_back( _prop_ );
      dcell->nref++;
      bc_dvs[mv->glb_id] = mv;		
      mv->on_boundary = true;
      mv->shared = true;
      mv->ghosted = false;
      mv->owner = MPI::processNumber();
      ref_edge[edge_key(v0->glb_id, v1->glb_id)] = mv;
    }
    
    closing = false;
  }

  // Create new cells
  DCell* c0 = new DCell;
  DCell* c1 = new DCell;
  c0->nref = dcell->nref;
  c1->nref = dcell->nref;
  std::vector<DVertex*> vs0(0);
  std::vector<DVertex*> vs1(0);
  for(uint i = 0; i < dcell->vertices.size(); i++)
  {
    if(i != ii)
      vs0.push_back(dcell->vertices[i]);

    if(i != jj)
      vs1.push_back(dcell->vertices[i]);
  } 
  vs0.push_back(mv);
  vs1.push_back(mv);

  addCell(c0, vs0, dcell->parent_id);
  addCell(c1, vs1, dcell->parent_id);    

  removeCell(dcell);

  // Continue refinement
  if(!closing)
  {
    // Bisect opposite cell of edge with hanging node
    for(;;)
    {
      DCell* copp = opposite(dcell, v0, v1);
      if(copp != 0)
      {
	bisect(copp, mv, v0, v1);
      }
      else
      {
	break;
      }
    }
  }
}
//-----------------------------------------------------------------------------
DCell* DMesh::opposite(DCell* dcell, DVertex* v1, DVertex* v2)
{
  for(std::list<DCell* >::iterator it = v1->cells.begin();
      it != v1->cells.end(); ++it)
  {
    DCell* c = *it;

    if(c != dcell && !c->deleted)
    {
      int matches = 0;
      for(uint i = 0; i < c->vertices.size(); i++)
      {
	if(c->vertices[i] == v1 || c->vertices[i] == v2)
	{
	  matches++;
	}
      }

      if(matches == 2)
      {
	// Found opposite cell
	return c;
      }
    }
  }  
  return 0;
}
//-----------------------------------------------------------------------------
void DMesh::addVertex(DVertex* v)
{
  //  vertices.push_back(v);
  vertices.insert(v);
  //  if(v->glb_id < 0)
  //    v->glb_id = _start_offset++;
}
//-----------------------------------------------------------------------------
void DMesh::addCell(DCell* c, std::vector<DVertex*> vs, int parent_id)
{
  for(uint i = 0; i < vs.size(); i++)
  {
    DVertex* v = vs[i];
    c->vertices.push_back(v);
    v->cells.push_back(c);
  }

  cells.push_back(c);
  c->parent_id = parent_id;
}
//-----------------------------------------------------------------------------
void DMesh::removeCell(DCell* c)
{
  /*
  for(uint i = 0; i < c->vertices.size(); ++i)
  {
    DVertex* v = c->vertices[i];
    v->cells.remove(c);
  } 
*/ 
  //  cells.remove(c);
  c->deleted = true;
  //  delete c;
}
//-----------------------------------------------------------------------------
void DMesh::bisectMarked(std::vector<bool> marked_ids)
{
  std::list<DCell*> marked_cells;
  for(std::list<DCell* >::iterator it = cells.begin();
      it != cells.end(); ++it)
  {
    DCell* c = *it;

    if(marked_ids[c->id])
    {
      marked_cells.push_back(c);
    }
  }

  for(std::list<DCell* >::iterator it = marked_cells.begin();
      it != marked_cells.end(); ++it)
  {
    DCell* c = *it;

    if(!c->deleted)
    {
      bisect(c, 0, 0, 0);
    }
  }

  std::vector<Propagation> propagated;
  std::list<Propagation> leftovers;

  bool empty = false;

  while(!empty) { 
    
    if(MPI::processNumber() == 0 && propagate.size() > 0)
      begin("Propagate refinement...");

    propagate_refinement(propagated, empty);

    if( empty && propagated.size() == 0) break;         
    propagate.clear();

    for(std::vector<Propagation>::iterator it = propagated.begin(); 
	it != propagated.end(); ++it) {

      DVertex* mv = 0;
      dolfin_assert(it->second.v1 != it->second.v2);
      if(ref_edge.find(edge_key(it->second.v1, it->second.v2)) != ref_edge.end()) {
        mv = ref_edge[edge_key(it->second.v1, it->second.v2)];

	if( mv->owner > (int) it->second.owner) {
	  mv->ghosted = true;
	  mv->shared = true;
          mv->owner = it->second.owner;
	}

	continue;
      }
      
      DVertex* v1 = 0;
      DVertex* v2 = 0;
      
      if(!v1 && bc_dvs.find(it->second.v1) != bc_dvs.end()) 
      {
	dolfin_assert(bc_dvs.find(it->second.v1) != bc_dvs.end());
	v1 = bc_dvs[it->second.v1];
      }
      if(!v2 && bc_dvs.find(it->second.v2) != bc_dvs.end())
      { 
	dolfin_assert(bc_dvs.find(it->second.v2) != bc_dvs.end());
	v2 = bc_dvs[it->second.v2];
      }

      if(!v1 || !v2)
      {
	leftovers.push_back(*it);    
	continue;
      }

      for(std::list<DCell* >::iterator ic = v1->cells.begin();
	  ic != v1->cells.end(); ++ic)       
      {
	if(!(*ic)->deleted) {
	  if((*ic)->has_edge(v1, v2))  {
	    dolfin_assert((*ic)->vertices.size() > 0);
	    if(mv == 0) 
	    {
	      mv = new DVertex;
	      mv->shared = true;
	      mv->glb_id = it->second.mv;				
	      vertices.insert(mv);	      

	      if( MPI::processNumber() < it->second.owner)
	      {
		mv->ghosted = false;
		mv->owner = MPI::processNumber();		  
		prop_edge node;
		node.mv = mv->glb_id;
		node.v1 = it->second.v1;
		node.v2 = it->second.v2;
		node.owner = mv->owner;
		std::pair<uint, prop_edge> prop(0, node);
		propagate.push_back(prop);
	      }
	      else 
	      {
		mv->ghosted = true;
		mv->owner = it->second.owner;
	      }
	      
	      mv->p = (v1->p + v2->p) / 2.0;
	      mv->on_boundary = true;
	      bc_dvs[mv->glb_id] = mv;
	      ref_edge[edge_key(v1->glb_id, v2->glb_id)] = mv;
	    }
	    dolfin_assert((*ic) > 0);
	    bisect((*ic), mv, v1, v2);
	  }
	}
      }
    }
    
    propagated.clear();
    for(std::list<Propagation>::iterator it = leftovers.begin(); 
	it != leftovers.end(); ++it)
      propagated.push_back(*it);
    leftovers.clear();

    if(MPI::processNumber() == 0)
      end();    
    
  }  
}
//-----------------------------------------------------------------------------
#ifdef HAS_MPI
//-----------------------------------------------------------------------------
void DMesh::propagate_naive(std::vector<Propagation>& propagated, bool& empty)
{
  // Allocate receive buffer
  int num_prop = propagate.size() * 5;
  int max_prop, recv_count;
  MPI_Allreduce(&num_prop, &max_prop, 1, MPI_INTEGER, MPI_MAX, MPI::DOLFIN_COMM);

  int *recv_buff = new int[max_prop];
  int *send_buff = new int[num_prop];
  int *sp = &send_buff[0];
  
  for(std::vector<Propagation>::iterator it = propagate.begin();
      it != propagate.end(); ++it) 
  {
    *(sp++) = it->first;
    *(sp++) = it->second.mv;
    *(sp++) = it->second.v1;
    *(sp++) = it->second.v2;
    *(sp++) = it->second.owner;
  }

  MPI_Status status;
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint dest, src;

  empty = true;
  for (uint j = 1; j < pe_size; j++) 
  {
    src = (rank -j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], num_prop, MPI_INTEGER, dest, 1,
		 recv_buff, max_prop, MPI_INTEGER, src, 1,
		 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_INTEGER, &recv_count);    

    if (recv_count > 0)
      empty = false;

    dolfin_assert(recv_count%5 == 0);
    for(int k = 0; k < recv_count; k += 5) 
    {   

      prop_edge node;
      node.mv = recv_buff[k+1];
      node.v1 = recv_buff[k+2];
      node.v2 = recv_buff[k+3];
      node.owner = recv_buff[k+4];
      
      Propagation prop(recv_buff[k], node);
      propagated.push_back(prop);
    }

  }

  less_pair comp;
  std::sort(propagated.begin(), propagated.end(), comp);
  
  short prop, gprop;
  prop = (empty == false);
  MPI_Allreduce(&prop, &gprop, 1, MPI_SHORT, MPI_SUM, MPI::DOLFIN_COMM);
  empty = (gprop == 0);
  
  delete[] send_buff;
  delete[] recv_buff;
}
//-----------------------------------------------------------------------------
void DMesh::propagate_hypercube(std::vector<Propagation>& propagated, 
				bool& empty)
{
 
  // Allocate receive buffer
  int num_prop = propagate.size() * 5;
  int total_prop, recv_count;
  MPI_Allreduce(&num_prop, &total_prop, 1,
		MPI_INTEGER, MPI_SUM, MPI::DOLFIN_COMM);

  int *recv_buff = new int[total_prop];
  int *state = new int[total_prop];
  int *sp = &state[0];
  uint state_size = 0;
  
  for(std::vector<Propagation>::iterator it = propagate.begin();
      it != propagate.end(); ++it) 
  {
    *(sp++) = it->first;
    *(sp++) = it->second.mv;
    *(sp++) = it->second.v1;
    *(sp++) = it->second.v2;
    *(sp++) = it->second.owner;
    state_size += 5;
  }

  MPI_Status status;
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint dest;
  uint D = 1;
#ifdef __sgi
  uint _log2, x;
  x = pe_size;
  _log2 = 0;
  while(x > 1){ _log2++; x>>=1;}
  for(uint j = 0; j < _log2; j++)
#else
  for(uint j = 0; j < log2(pe_size) ; j++)
#endif
  {
    dest = rank^(D<<j);

    MPI_Sendrecv(state, state_size, MPI_INTEGER, dest, 1, 
		 recv_buff, total_prop, MPI_INTEGER, dest, 1,
		 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_INTEGER, &recv_count);
    
    dolfin_assert(recv_count%5 == 0);
    for(int k = 0; k < recv_count; k += 5) 
    {   

      prop_edge node;
      node.mv = recv_buff[k+1];
      node.v1 = recv_buff[k+2];
      node.v2 = recv_buff[k+3];
      node.owner = recv_buff[k+4];
      
      Propagation prop(recv_buff[k], node);
      propagated.push_back(prop);
    }
    memcpy(sp, recv_buff, recv_count*sizeof(int));
    sp += recv_count;
    state_size += recv_count;

  }

  less_pair comp;
  std::sort(propagated.begin(), propagated.end(), comp);
  empty = (state_size == 0);

  delete[] recv_buff;
  delete[] state;
}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
void DMesh::propagate_naive(std::vector<Propagation>& propagated, bool& empty) 
{
  error("Rivara needs MPI");
}
//-----------------------------------------------------------------------------
void DMesh::propagate_hypercube(std::vector<Propagation>& propagated, bool& empty)
{
  error("Rivara needs MPI");
}
//-----------------------------------------------------------------------------
#endif
