// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009.
//
#include "RivaraRefinement.h"
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
#include "Vertex.h"
#include "Facet.h"
#include "Edge.h"
#include "Cell.h"
#include "BoundaryMesh.h"
#include "RivaraRefinement.h"
#include "LoadBalancer.h"



#ifdef HAS_MPI
#include <mpi.h>
#endif 

using namespace dolfin;
//-----------------------------------------------------------------------------
void RivaraRefinement::refine(Mesh& mesh, 
			      MeshFunction<bool>& cell_marker,
			      MeshFunction<uint>& cell_map,
			      real tf, real tb, real ts)
{
  message("Refining simplicial mesh by recursive Rivara bisection.");

  // Start Loadbalancer
  if(MPI::numProcesses() > 1) {
    begin("Load balancing");
    //    Tune loadbalancer using machine specific parameters, if available
    if( tf > 0.0 && tb > 0.0 && ts > 0.0)
      LoadBalancer::balance(mesh, cell_marker, tf, tb, ts);
    else
      LoadBalancer::balance(mesh, cell_marker,LoadBalancer::Rivara);
    end();
  }
  mesh.renumber();

  int d = mesh.topology().dim();

  // Dynamic mesh test
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
  
  std::vector<int> new2old_cell_arr;
  
  Mesh omesh;
  
  dmesh.exp(omesh, new2old_cell_arr);
  
  mesh = omesh;
  
  // Generate mesh function map
  cell_map.init(mesh, d);
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    cell_map.set(c->index(), new2old_cell_arr[c->index()]);
  }
}
//-----------------------------------------------------------------------------
DVertex::DVertex() : id(0), glb_id(-1), cells(0), p(0.0, 0.0, 0.0), 
		     on_boundary(false), shared(false),
		     ghosted(false), owner(0)
{
}
//-----------------------------------------------------------------------------
DCell::DCell() : id(0), parent_id(0), vertices(0), deleted(false)
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

  std::vector<DVertex *> vertexvec;

  BoundaryMesh boundary;
  boundary.init_interior(mesh);
  File bc_m("boundary.pvd");
  bc_m << boundary;
  MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");
  
  MeshFunction<bool> on_boundary(mesh, 0);
  on_boundary = false;
  MeshFunction<bool> boundary_cell(mesh, mesh.topology().dim());
  boundary_cell = false;
 // Generate facet - cell connectivity if not generated
  mesh.init(mesh.topology().dim() - 1, mesh.topology().dim());
  for (CellIterator bf(boundary); !bf.end(); ++bf) 
  {
    Facet f(mesh, cell_map->get(*bf));    
    for (CellIterator c(f); !c.end(); ++c) 
    {
      boundary_cell.set(*c, true);	  
      for(EdgeIterator e(*c); !e.end(); ++e) 
      {
	const uint *edge_v = e->entities(0);
	if(mesh.distdata().is_shared(edge_v[0], 0) ||
	   mesh.distdata().is_shared(edge_v[1], 0)) 
	{
	  on_boundary.set(edge_v[0], true);
	  on_boundary.set(edge_v[1], true);
	}
      }
    }     
  }

  // Assume uniform refinement
  uint num_new = mesh.size(1);
  num_new *= 5;
  // Find maximum global index assigned
  uint max_index = std::max(mesh.distdata().global_numVertices(),
			    mesh.distdata().max_index());
  uint glb_max;
  MPI_Allreduce(&max_index, &glb_max, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  
  // Assign a safe range for each processor
  _start_offset = 0;
#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_new, &_start_offset, 1,
	     MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Scan(&num_new, &_start_offset, 1,
	   MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  _start_offset -= num_new;
#endif
  _start_offset += glb_max;
  
  uint counter = 1;
  for (VertexIterator vi(mesh); !vi.end(); ++vi)
  {
    DVertex* dv = new DVertex;    
    dv->p = vi->point();
    dv->glb_id = mesh.distdata().get_global(vi->index(), 0);
    dv->on_boundary = mesh.distdata().is_shared(vi->index(), 0);
    //    dv->on_boundary = on_boundary.get(vi->index());
    dv->shared = mesh.distdata().is_shared(vi->index(), 0);
    dv->ghosted = mesh.distdata().is_ghost(vi->index(), 0);
    if (dv->ghosted)
      dv->owner = mesh.distdata().get_owner(*vi);    

    glb_ids.insert(dv->glb_id);    
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
void DMesh::exp(Mesh& mesh, std::vector<int>& new2old_cell)
{
  number();

  new2old_cell.resize(cells.size());

  mesh.clear();
  mesh.distdata().clear();

  MeshEditor editor;
  Mesh newmesh;
  editor.open(newmesh, cell_type->cellType(), d, d);
  
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
      newmesh.distdata().set_ghost(current_vertex, 0);
      newmesh.distdata().set_ghost_owner(current_vertex, dv->owner, 0);
    }
    if(dv->shared)
      newmesh.distdata().set_shared(current_vertex, 0);
    newmesh.distdata().set_map(current_vertex++, dv->glb_id, 0);

      
  }

  Array<uint> cell_vertices(cell_type->numEntities(0));
  uint current_cell = 0;
  for(std::list<DCell* >::iterator it = cells.begin();
      it != cells.end(); ++it)
  {
    DCell* dc = *it;
    if(dc->deleted)
      error("Deleted");

    for(uint j = 0; j < dc->vertices.size(); j++)
    {
      DVertex* dv = dc->vertices[j];
      cell_vertices[j] = dv->id;
    }
    editor.addCell(current_cell, cell_vertices);
    new2old_cell[current_cell] = dc->parent_id;
    
    current_cell++;
  }
  editor.close();

  mesh = newmesh;
  mesh.distdata().invalid_numbering();
  mesh.distdata().invalid_ownership();
  mesh.renumber();
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

	real l = v0->p.distance(v1->p);

	if(fabs(l - lmax) < DOLFIN_EPS)
	{
	  int ptsum = ((long)v0) + ((long)v1);
	  if(ptsum > ptmax)
	  {
	    ii = i;
	    jj = j;
	    lmax = l;
	    ptmax = ((long)v0 + (long)v1);
	  }
	}
	else if(l >= lmax)
 	{
 	  ii = i;
 	  jj = j;
 	  lmax = l;
	  ptmax = ((long)dcell->vertices[i]) + ((long)dcell->vertices[j]);
 	}
       }
     }
  }

  dolfin_assert(dcell->vertices.size() > 0);
  dolfin_assert(dcell->vertices.size() > 0);
  dolfin_assert(dcell->vertices.size() > ii);
  dolfin_assert(dcell->vertices.size() > jj);
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
      
      dolfin_assert(v0->glb_id != v1->glb_id);
      if( ref_edge.find(edge_key(v0->glb_id, v1->glb_id)) == ref_edge.end())
	ref_edge[edge_key(v0->glb_id, v1->glb_id)] = mv;
    }
  }
  else
  {
    mv = new DVertex;
    mv->p = (dcell->vertices[ii]->p + dcell->vertices[jj]->p) / 2.0;

    addVertex(mv);

    // Add hanging node on shared edges to propagation buffer
    if( v0->on_boundary && v1->on_boundary) 
    {
	propagate.push_back(mv->glb_id);
	propagate.push_back(v0->glb_id);
	propagate.push_back(v1->glb_id);
	propagate.push_back(MPI::processNumber());
	bc_dvs[mv->glb_id] = mv;		
	mv->on_boundary = true;
	mv->shared = true;
	mv->owner = MPI::processNumber();
    }
    
    closing = false;
  }

  // Create new cells
  DCell* c0 = new DCell;
  DCell* c1 = new DCell;
  std::vector<DVertex*> vs0(0);
  std::vector<DVertex*> vs1(0);
  std::vector<uint> sh0;
  std::vector<uint> sh1;
  for(uint i = 0; i < dcell->vertices.size(); i++)
  {
    if(i != ii)
    {
      vs0.push_back(dcell->vertices[i]);
      if( dcell->vertices[i]->on_boundary )
	sh0.push_back(dcell->vertices[i]->glb_id);
    }
    if(i != jj)
    {
      vs1.push_back(dcell->vertices[i]);
      if( dcell->vertices[i]->on_boundary )
	sh1.push_back(dcell->vertices[i]->glb_id);
    }
  } 
  vs0.push_back(mv);
  vs1.push_back(mv);

  addCell(c0, vs0, dcell->parent_id);
  addCell(c1, vs1, dcell->parent_id);
  if( v0->on_boundary && v1->on_boundary) 
  {
    EdgeKey key = edge_key(v0->glb_id, v1->glb_id);
    ref_edge[key] = mv;
  } 

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
  {// dolfin_assert((*it)); dolfin_assert(!(*it)->deleted); 
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
  if(v->glb_id < 0)
    v->glb_id = _start_offset++;
  glb_ids.insert(v->glb_id);
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

  uint pre_num_cells = cells.size();
  std::vector<uint> propagated;
  bool empty = false;
  MPI_Barrier(MPI_COMM_WORLD);
  while(!empty) { 
    
    if(MPI::processNumber() == 0 && propagate.size() > 0)
      begin("Propagate refinement...");
    propagated.clear();
    propagate_naive( propagated, empty);
    
    propagate.clear();
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPI::processNumber() == 0 && propagated.size() > 0) {
      printf("Bisecting...");
      fflush(stdout);
    }
    uint cc = 0;
    dolfin_assert(propagated.size() % 4 == 0);
    for(uint i = 0; i < propagated.size();  i += 4) {

      if(ref_edge.find(edge_key(propagated[i+1], propagated[i+2])) != ref_edge.end()) {
	DVertex *mv = ref_edge[edge_key(propagated[i+1], propagated[i+2])];
	
	//	if ( propagated[i+3] < MPI::processNumber() ||
	//	     (mv->ghosted && (mv->owner > propagated[i+3]))) {
	//	if( mv->shared && mv->owner > propagated[i+3]) {
	if( mv->owner >= propagated[i+3]) {
		  
	  bc_dvs.erase(mv->glb_id);
	  glb_ids.erase(mv->glb_id);	  
	  mv->glb_id = propagated[i];
	  glb_ids.insert(mv->glb_id);
	  bc_dvs[propagated[i]] = mv;
	  mv->ghosted = true;	
	  mv->shared = true;
	  mv->owner = propagated[i+3];

	}
	continue;
      }

      if(MPI::processNumber() == 0) {
	switch(cc)
	{
	case 0:
	  putchar('-'); break;
	case 1:
	  putchar('\\'); break;
	case 2:
	  putchar('|'); break;
	case 3:
	  putchar('/'); break;	  
	}
	fflush(stdout);
	usleep(100000);
	putchar('\b');
      }
      cc = (cc + 1)%4;
      
      DVertex* mv = 0;

      dolfin_assert(bc_dvs.find(propagated[i+1]) != bc_dvs.end());
      dolfin_assert(bc_dvs.find(propagated[i+2]) != bc_dvs.end());
      DVertex* v1 = bc_dvs[propagated[i+1]];
      DVertex* v2 = bc_dvs[propagated[i+2]];

      for(std::list<DCell* >::iterator it = v1->cells.begin();
	  it != v1->cells.end(); ++it)       
      {
	if(!(*it)->deleted) {
	  if((*it)->has_edge(v1, v2))  {
	    dolfin_assert((*it)->vertices.size() > 0);
	    if(mv == 0) 
	    {
	      mv = new DVertex;
	      mv->glb_id = propagated[i];
	      mv->p = (bc_dvs[propagated[i+1]]->p + bc_dvs[propagated[i+2]]->p) / 2.0;
	      mv->on_boundary = true;
	      mv->ghosted = true;
	      mv->shared = true;
	      mv->owner = propagated[i+3];
	      glb_ids.insert(mv->glb_id);
	      bc_dvs[mv->glb_id] = mv;
	      vertices.insert(mv);
	      ref_edge[edge_key(propagated[i+1], propagated[i+2])] = mv;
	    }
	    fflush(stdout);
	    dolfin_assert((*it) > 0);
	    bisect((*it), mv, v1, v2);
	  }
	}
      }
    }
    
    if(MPI::processNumber() == 0 && propagated.size() > 0)
      putchar('\n');

    if(MPI::processNumber() == 0){
      message("Done");
      end();    
    }

    uint num_prop = propagate.size();
    uint num_gprop = 0;
    MPI_Allreduce(&num_prop, &num_gprop, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    empty = (num_gprop == 0);
  }
    
  MPI_Barrier(MPI_COMM_WORLD);
  message("Propagated refinements: %d", (cells.size() - pre_num_cells) / 2);
}
//-----------------------------------------------------------------------------
void DMesh::propagate_naive(std::vector<uint>& propagated, bool& empty)
{

  empty = true;
  // Allocate receive buffer
  int num_prop = propagate.size();
  int max_prop, recv_count;
  MPI_Allreduce(&num_prop, &max_prop, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

  int *recv_buff = new int[max_prop];

  MPI_Status status;
  int src,dest;
  int pe_size = MPI::numProcesses();
  int rank = MPI::processNumber();

  for (int j = 1; j < pe_size; j++) 
  {
    src = (rank -j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;
    
    /*
    MPI_Sendrecv(&propagate[0], propagate.size(), MPI_UNSIGNED, dest, 1,
		 recv_buff, max_prop, MPI_UNSIGNED, src, 1,
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);    
    */
    MPI_Sendrecv(&propagate[0], propagate.size(), MPI_INTEGER, dest, 1,
		 recv_buff, max_prop, MPI_INTEGER, src, 1,
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_INTEGER, &recv_count);    
    if (recv_count > 0)
      empty = false;
    
    for (int k = 0; k < recv_count; k += 4) 
    {
      if( glb_ids.find(recv_buff[k+1]) != glb_ids.end() &&
	  glb_ids.find(recv_buff[k+2]) != glb_ids.end())
      {
	for (int i = 0; i < 4; i++)
	  propagated.push_back(recv_buff[k+i]);
	//	propagated.push_back(status.MPI_SOURCE);
      }
    }
  }

  delete[] recv_buff;
}
//-----------------------------------------------------------------------------
void DMesh::propagate_hypercube(std::vector<uint>& propagated)
{
  // Implement hypercube exchange
  error("Not implemented");
}
//-----------------------------------------------------------------------------
