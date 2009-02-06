// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#include "RivaraRefinement.h"

#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/BoundaryMesh.h>

#include <mpi.h>

using namespace dolfin;


//-----------------------------------------------------------------------------
void RivaraRefinement::refine(Mesh& mesh, 
			      MeshFunction<bool>& cell_marker,
			      MeshFunction<uint>& cell_map)
{
  message("Refining simplicial mesh by recursive Rivara bisection.");

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
DVertex::DVertex() : id(0), glb_id(0), cells(0), p(0.0, 0.0, 0.0), 
		     on_boundary(false)
{
}
//-----------------------------------------------------------------------------
DCell::DCell() : id(0), parent_id(0), vertices(0), deleted(false)
{
}
//-----------------------------------------------------------------------------
DMesh::DMesh() : vertices(0), cells(0)
{
}
//-----------------------------------------------------------------------------
DMesh::~DMesh()
{
  // Delete allocated DVertices
  for(std::list<DVertex* >::iterator it = vertices.begin();
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
  MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");
  
  MeshFunction<bool> on_boundary(mesh, 0);
  on_boundary = false;
  MeshFunction<bool> boundary_cell(mesh, 2);
  boundary_cell = false;

  for (CellIterator bf(boundary); !bf.end(); ++bf) 
  {
    Facet f(mesh, cell_map->get(*bf));    
    for (CellIterator c(f); !c.end(); ++c) 
    {
      boundary_cell.set(*c, true);
      for(EdgeIterator e(*c); !e.end(); ++e) 
      {
	const uint *edge_v = e->entities(0);
	if(mesh.distdata().is_shared(edge_v[0], 0) &&
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

  for (VertexIterator vi(mesh); !vi.end(); ++vi)
  {
    DVertex* dv = new DVertex;    
    dv->p = vi->point();
    dv->glb_id = mesh.distdata().get_global(*vi);
    dv->on_boundary = on_boundary.get(*vi);
    glb_ids.insert(dv->glb_id);    
    if(dv->on_boundary) {

      bc_dvs[dv->glb_id] = dv;
    }

    addVertex(dv);
    vertexvec.push_back(dv);
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
	if( on_boundary.get(edge_v[0]) && on_boundary.get(edge_v[1])) 
	{
	  EdgeKey key = edge_key(mesh.distdata().get_global(edge_v[0], 0),
				 mesh.distdata().get_global(edge_v[1], 0));
	  bc_dcs[key] = dc;
	  break;
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

  MeshEditor editor;
  editor.open(mesh, cell_type->cellType(),
	      d, d);
  
  editor.initVertices(vertices.size());
  editor.initCells(cells.size());

  // Add old vertices
  uint current_vertex = 0;
  for(std::list<DVertex* >::iterator it = vertices.begin();
      it != vertices.end(); ++it)
  {
    DVertex* dv = *it;
    editor.addVertex(current_vertex++, dv->p);
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
    new2old_cell[current_cell] = dc->parent_id;
    
    current_cell++;
  }
  editor.close();
}
//-----------------------------------------------------------------------------
void DMesh::number()
{
  uint i = 0;
  for(std::list<DVertex* >::iterator it = vertices.begin();
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
void DMesh::bisect(DCell* dcell, DVertex* hangv,
		   DVertex* hv0, DVertex* hv1, bool prop)
{
  //cout << "Refining cell: " << dcell->id << endl;

  bool closing = false;

  // Find longest edge
  real lmax = 0.0;
  uint ii = 0;
  uint jj = 0;
  for(uint i = 0; i < dcell->vertices.size(); i++)
  {
    for(uint j = 0; j < dcell->vertices.size(); j++)
    {
      if(i != j)
      {
	real l = dcell->vertices[i]->p.distance(dcell->vertices[j]->p);
	if(l >= lmax)
	{
	  ii = i;
	  jj = j;
	  lmax = l;
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
      mv->on_boundary = true;
      bc_dvs[mv->glb_id] = mv;
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
    bc_dcs.erase(edge_key(v0->glb_id, v1->glb_id));
    bc_dcs[edge_key(v0->glb_id, mv->glb_id)] = c0;
    bc_dcs[edge_key(v1->glb_id, mv->glb_id)] = c1;

    ref_edge.insert(edge_key(v0->glb_id, mv->glb_id));
    ref_edge.insert(edge_key(v1->glb_id, mv->glb_id));
  } 
  
  
  for(uint i = 0; i < vs0.size(); i++)
    for(uint j = 0; j < vs0.size(); j++)
      if ( i != j) {
	if(vs0[i]->on_boundary && vs0[j]->on_boundary) {
	  EdgeKey key = edge_key(vs0[i]->glb_id, vs0[j]->glb_id);
	  bc_dcs.erase(key);
	  bc_dcs[key] = c0;
	}
      }

  
  for(uint i = 0; i < vs1.size(); i++)
    for(uint j = 0; j < vs1.size(); j++)
      if ( i != j) {
	if(vs1[i]->on_boundary && vs1[j]->on_boundary) {
	  EdgeKey key = edge_key(vs1[i]->glb_id, vs1[j]->glb_id);
	  bc_dcs.erase(key);
	  bc_dcs[key] = c1;
	}
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
  {
    DCell* c = *it;

    if(c != dcell)
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
  vertices.push_back(v);
  if(v->glb_id == 0)
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
  for(uint i = 0; i < c->vertices.size(); ++i)
  {
    DVertex* v = c->vertices[i];
    v->cells.remove(c);
  }  
  //cells.remove(c);
  c->deleted = true;
  //delete c;
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
  bool empty;
  do {    
    if(MPI::processNumber() == 0 && propagate.size() > 0)
      begin("Propagate refinement...");
    propagated.clear();
    propagate_naive( propagated, empty);
    propagate.clear();
    if(MPI::processNumber() == 0 && propagated.size() > 0) {
      printf("Bisecting...");
      fflush(stdout);
    }
    uint cc = 0;
    for(uint i = 0; i < propagated.size();  i += 3) {
      if( bc_dcs.find(edge_key(propagated[i+1], propagated[i+2])) == bc_dcs.end())
	continue;      

      DVertex *mv = new DVertex;
      vertices.push_back(mv);
      mv->glb_id = propagated[i];
      mv->p = (bc_dvs[propagated[i+1]]->p + bc_dvs[propagated[i+2]]->p) / 2.0;
      mv->on_boundary = true;
      glb_ids.insert(mv->glb_id);
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
      cc++;
      cc = cc%4;
      //cout<< "."; 
      //      if(ref_edge.find(edge_key(propagated[i+1], propagated[i+2])) == ref_edge.end())
      bisect(bc_dcs[edge_key(propagated[i+1], propagated[i+2])], 
	     mv, bc_dvs[propagated[i+1]], bc_dvs[propagated[i+2]]);
      
      //	bisect(bc_dcs[edge_key(propagated[i+1], propagated[i+2])], 
      //	       mv, 0, 0);
    }
    if(MPI::processNumber() == 0 && propagated.size() > 0) {
      putchar('\n');
    }
    if(MPI::processNumber() == 0)
      end();    
  } while(propagated.size() > 0 && empty);
  MPI_Barrier(MPI_COMM_WORLD);
  message("Propagated refinements: %d", (cells.size() - pre_num_cells) / 2);
}
//-----------------------------------------------------------------------------
void DMesh::propagate_naive(std::vector<uint>& propagated, bool empty)
{

  empty = true;
  // Allocate receive buffer
  int num_prop = propagate.size();
  int max_prop, recv_count;
  MPI_Allreduce(&num_prop, &max_prop, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

  uint *recv_buff = new uint[max_prop];

  MPI_Status status;
  int src,dest;
  int pe_size = MPI::numProcesses();
  int rank = MPI::processNumber();

  for (int j = 1; j < pe_size; j++) 
  {
    src = (rank -j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&propagate[0], propagate.size(), MPI_UNSIGNED, dest, 1,
		 recv_buff, max_prop, MPI_UNSIGNED, src, 1,
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);    
    if (recv_count > 0)
      empty = false;
    
    for (int k = 0; k < recv_count; k += 3) 
    {
      if( glb_ids.find(recv_buff[k+1]) != glb_ids.end() &&
	  glb_ids.find(recv_buff[k+2]) != glb_ids.end()) 
	for (int i = 0; i < 3; i++)
	  propagated.push_back(recv_buff[k+i]);
    }
  }
}
//-----------------------------------------------------------------------------
void DMesh::propagate_hypercube(std::vector<uint>& propagated)
{
  // Implement hypercube exchange
  error("Not implemented");
}
//-----------------------------------------------------------------------------
