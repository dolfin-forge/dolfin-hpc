// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-09
// Last changed: 2009-09-13

#include <fstream>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/function/Function.h>
#include "Checkpoint.h"


using namespace dolfin;
//-----------------------------------------------------------------------------
Checkpoint::Checkpoint() : state(CHECKPOINT), restart_state(OPEN)
{
}
//-----------------------------------------------------------------------------
Checkpoint::~Checkpoint()
{
  in.close();
}
//-----------------------------------------------------------------------------
void Checkpoint::write(real t, Mesh& mesh, std::vector<Function *> func)
{

  message("Writing checkpoint at time %g", t);
  std::ostringstream fname;
  fname << "checkpoint";
  if( MPI::numProcesses() > 1) 
    fname << "_" <<  MPI::processNumber();
  fname << ".chkp";

  std::ofstream out(fname.str().c_str(), std::ofstream::binary);

  out.write((char *) &t, sizeof(real));
  write_mesh(mesh, out);

  write_func(func, out);

  out.close();

  state = RESTART;

}
//-----------------------------------------------------------------------------
void Checkpoint::restart(std::string fname)
{

  if (restart_state != OPEN)
    error("Shut her down, Scotty, she's sucking mud again!");

  std::ostringstream _fname;
  if( MPI::numProcesses() > 1) 
    _fname << fname << "_" << MPI::processNumber() << ".chkp";
  else
    _fname << fname << ".chkp";
  
  in.open(_fname.str().c_str(), std::ifstream::binary);
  in.read((char *) &_t, sizeof(real));
  message("Restarting from time %g", _t);
  state = RESTART;
  restart_state = MESH;
      
}
//-----------------------------------------------------------------------------
void Checkpoint::load(Mesh& mesh)
{ 
  if (restart_state != MESH)
    error("Shut her down, Scotty, she's sucking mud again!");
  
  CellType::Type type;
  uint tdim, gdim, num_vertices, num_cells, num_entities;
  in.read((char *)&type, sizeof(CellType::Type));
  in.read((char *)&tdim, sizeof(uint));
  in.read((char *)&gdim, sizeof(uint));
  in.read((char *)&num_vertices, sizeof(uint));
  in.read((char *)&num_cells, sizeof(uint));
  in.read((char *)&num_entities, sizeof(uint));
  
  real *coords = new real[gdim * num_vertices];
  in.read((char *)coords, (gdim * num_vertices) * sizeof(real));

  MeshEditor editor;
  editor.open(mesh,type, tdim, gdim);
  editor.initVertices(num_vertices);

  uint vi = 0;
  for(uint i = 0 ; i < gdim * num_vertices; i += gdim)
  {
    switch(gdim)
    {      
    case 2:
      editor.addVertex(vi++, coords[i], coords[i+1]); break;
    case 3:
      editor.addVertex(vi++, coords[i], coords[i+1], coords[i+2]); break;
    }
  }
  
  delete[] coords;

  editor.initCells(num_cells);

  uint *cells = new uint[num_entities * num_cells];
  in.read((char *)cells, (num_entities * num_cells) * sizeof(uint));


  Array<uint> v;
  uint ci = 0;
  for (uint i = 0; i < num_entities * num_cells; i += num_entities)
  {
    v.clear();
    for (uint j = 0; j < num_entities; j++)
      v.push_back(cells[i + j]);
    editor.addCell(ci++, v);    
  }  
  editor.close();
  delete[] cells;


  if (MPI::numProcesses() > 1) 
  {
    uint *mapping = new uint[mesh.numVertices()];
    in.read((char *)mapping, mesh.numVertices() * sizeof(uint));
    for (VertexIterator v(mesh); !v.end(); ++v)
      mesh.distdata().set_map(v->index(), mapping[v->index()], 0);
    delete[] mapping;

    uint num_ghost;
    in.read((char *)&num_ghost, sizeof(uint));
    uint *ghosts = new uint[2 * num_ghost];
    in.read((char *)ghosts, 2*num_ghost * sizeof(uint));
    for (uint i = 0; i < 2 * num_ghost; i += 2)
    {
      mesh.distdata().set_ghost(ghosts[i], 0);
      mesh.distdata().set_ghost_owner(ghosts[i], ghosts[i+1], 0);
    }
    delete[] ghosts;


    uint num_shared;
    in.read((char *)&num_shared, sizeof(uint));
    
    uint *shared = new uint[num_shared];
    in.read((char *)shared, num_shared * sizeof(uint));
    for(uint i = 0; i < num_shared; i++)
      mesh.distdata().set_shared(shared[i], 0);
    delete[] shared;

  }
  
  

  restart_state = FUNC;
}
//-----------------------------------------------------------------------------
void Checkpoint::load(std::vector<Function *> func)
{
  if (restart_state != FUNC)
    error("Shut her down, Scotty, she's sucking mud again!");


  std::vector<Function *>::iterator it;
  uint local_size;
  // FIXME store max(local_size) 
  for (it = func.begin(); it != func.end(); ++it)
  {
    in.read((char *)&local_size, sizeof(uint));
    real *values = new real[local_size];
    in.read((char *)values, local_size * sizeof(real));    
    (*it)->vector().set(values);
    (*it)->vector().apply();
    delete[] values;
  }

}
//-----------------------------------------------------------------------------
void Checkpoint::write_mesh(Mesh& mesh, std::ofstream& out)
{

  uint num_coords = mesh.numVertices() * mesh.geometry().dim();
  uint num_entities = mesh.type().numEntities(0);
  uint num_centities = mesh.numCells() * num_entities;
  CellType::Type type = mesh.type().cellType();
  uint tdim = mesh.topology().dim();
  uint gdim = mesh.geometry().dim();
  uint num_vertices = mesh.numVertices();
  uint num_cells = mesh.numCells();
  
  out.write((char *)&type, sizeof(CellType::Type));
  out.write((char *)&tdim, sizeof(uint));
  out.write((char *)&gdim, sizeof(uint));
  out.write((char *)&num_vertices, sizeof(uint));
  out.write((char *)&num_cells, sizeof(uint));  
  out.write((char *)&num_entities, sizeof(uint));
  out.write((char *)mesh.coordinates(), num_coords * sizeof(real));
  out.write((char *)mesh.cells(), num_centities * sizeof(uint));

  if (MPI::numProcesses() > 1) 
  {
    uint *mapping = new uint[mesh.numVertices()];
    for (VertexIterator v(mesh); !v.end(); ++v)
      mapping[v->index()] = mesh.distdata().get_global(*v);
    out.write((char *)mapping, mesh.numVertices() * sizeof(uint));
    delete[] mapping;


    uint num_ghost = mesh.distdata().num_ghost(0);
    out.write((char *)&num_ghost, sizeof(uint));
    uint *ghosts = new uint[2 * num_ghost];
    uint *gp = &ghosts[0];
    for (MeshGhostIterator g(mesh.distdata(), 0); !g.end(); ++g)
    {
      *gp++ = g.index();
      *gp++ = g.owner();
    }
    out.write((char *)ghosts, 2 * num_ghost * sizeof(uint));
    delete[] ghosts;

    uint num_shared = mesh.distdata().num_shared(0);
    out.write((char *)&num_shared, sizeof(uint));

    uint *shared = new uint[num_shared];
    uint *sp = &shared[0];
    for (MeshSharedIterator s(mesh.distdata(), 0); !s.end(); ++s)
      *sp++ = s.index();
    out.write((char *)shared, mesh.distdata().num_shared(0) * sizeof(uint));
    delete[] shared;

  }

}
//-----------------------------------------------------------------------------
void Checkpoint::write_func(std::vector<Function *> func, std::ofstream& out)
{
  std::vector<Function *>::iterator it;

  uint max_size = 0;
  for (it = func.begin(); it != func.end(); ++it)
  {
    if((*it)->type() != Function::discrete)
      error("Checkpoint restart only implemented for discrete functions");    
    max_size = std::max(max_size, (*it)->vector().local_size());
  }
  
  message("Max local_size: %d", max_size);
  real *values = new real[max_size];
  for (it = func.begin(); it != func.end(); ++it)
  {
    uint local_size = (*it)->vector().local_size();
    (*it)->vector().get(values);
    out.write((char *)&local_size, sizeof(uint));
    out.write((char *)values, (*it)->vector().local_size() * sizeof(real));
  }
  delete[] values;
  
}
//-----------------------------------------------------------------------------
