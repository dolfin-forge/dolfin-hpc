// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-09
// Last changed: 2010-06-20

#include <sstream>
#include <fstream>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/Checkpoint.h>

#ifdef ENABLE_MPIIO
#include <mpi.h>
#endif


using namespace dolfin;
//-----------------------------------------------------------------------------
Checkpoint::Checkpoint() : state(CHECKPOINT), restart_state(OPEN), n(0)
{
}
//-----------------------------------------------------------------------------
Checkpoint::~Checkpoint()
{
}
//-----------------------------------------------------------------------------
void Checkpoint::write(std::string fname, uint id, real t, Mesh& mesh, 
		       std::vector<Function *> func,
		       std::vector<Vector *> vec)
{

  message("Writing checkpoint (%s%d) at time %g", fname.c_str(), n%2, t);
  std::ostringstream _fname;
#ifndef ENABLE_MPIIO
  if( MPI::numProcesses() > 1) 
    _fname << fname << (n++)%2 << "_" <<  MPI::processNumber() << ".chkp";
  else
    _fname << fname << (n++)%2 << ".chkp";

  std::ofstream out(_fname.str().c_str(), std::ofstream::binary);

  out.write((char *) &id, sizeof(uint));
  out.write((char *) &t, sizeof(real));

#else
  _fname << fname << (n++)%2 << ".chkp";
    
  MPI_File out;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) _fname.str().c_str(),
		MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out);
  
  byte_offset = 0;
  MPI_File_write_all(out, &id, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  MPI_File_write_all(out, &t, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

  byte_offset += sizeof(uint);
  byte_offset += sizeof(real);
#endif

  write(mesh, out);
  write(func, out);
  write(vec, out);
  
#ifdef ENABLE_MPIIO
  MPI_File_close(&out);
#else
  out.close();
#endif


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
  in.read((char *) &_id, sizeof(uint));
  in.read((char *) &_t, sizeof(real));
  message("Restarting from time %g checkpoint id %d", _t, _id);

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

  Mesh _mesh;
  MeshEditor editor;  
  editor.open(_mesh,type, tdim, gdim);
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
    uint *mapping = new uint[_mesh.numVertices()];
    in.read((char *)mapping, _mesh.numVertices() * sizeof(uint));
    for (VertexIterator v(_mesh); !v.end(); ++v)
      _mesh.distdata().set_map(v->index(), mapping[v->index()], 0);
    delete[] mapping;

    uint num_ghost;
    in.read((char *)&num_ghost, sizeof(uint));
    uint *ghosts = new uint[2 * num_ghost];
    in.read((char *)ghosts, 2*num_ghost * sizeof(uint));
    for (uint i = 0; i < 2 * num_ghost; i += 2)
    {
      _mesh.distdata().set_ghost(ghosts[i], 0);
      _mesh.distdata().set_ghost_owner(ghosts[i], ghosts[i+1], 0);
    }
    delete[] ghosts;


    uint num_shared;
    in.read((char *)&num_shared, sizeof(uint));
    
    uint *shared = new uint[num_shared];
    in.read((char *)shared, num_shared * sizeof(uint));
    for(uint i = 0; i < num_shared; i++)
      _mesh.distdata().set_shared(shared[i], 0);
    delete[] shared;    
  }
  
  mesh = _mesh;
  mesh.distdata().invalid_numbering();
  mesh.renumber();

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

  restart_state = VEC;
}
//-----------------------------------------------------------------------------
void Checkpoint::load(std::vector<Vector *> vec)
{
  if (restart_state != VEC)
    error("Shut her down, Scotty, she's sucking mud again!");
  
  std::vector<Vector *>::iterator it;
  uint local_size;
  for (it = vec.begin(); it != vec.end(); ++it)
  {
    in.read((char *)&local_size, sizeof(uint));
    real *values = new real[local_size];
    in.read((char *)values, local_size * sizeof(real));    
    (*it)->set(values);
    (*it)->apply();
    delete[] values;
  }
  
  in.close();
}
//-----------------------------------------------------------------------------
void Checkpoint::write(Mesh& mesh, chkp_outstream& out)
{

  uint num_coords = mesh.numVertices() * mesh.geometry().dim();
  uint num_entities = mesh.type().numEntities(0);
  uint num_centities = mesh.numCells() * num_entities;
  CellType::Type type = mesh.type().cellType();
  uint tdim = mesh.topology().dim();
  uint gdim = mesh.geometry().dim();
  uint num_vertices = mesh.numVertices();
  uint num_cells = mesh.numCells();
  
#ifdef ENABLE_MPIIO
#else
  out.write((char *)&type, sizeof(CellType::Type));
  out.write((char *)&tdim, sizeof(uint));
  out.write((char *)&gdim, sizeof(uint));
  out.write((char *)&num_vertices, sizeof(uint));
  out.write((char *)&num_cells, sizeof(uint));  
  out.write((char *)&num_entities, sizeof(uint));
  out.write((char *)mesh.coordinates(), num_coords * sizeof(real));
  out.write((char *)mesh.cells(), num_centities * sizeof(uint));
#endif

  if (MPI::numProcesses() > 1) 
  {
    uint *mapping = new uint[mesh.numVertices()];
    for (VertexIterator v(mesh); !v.end(); ++v)
      mapping[v->index()] = mesh.distdata().get_global(*v);
    
#ifdef ENABLE_MPIIO
#else
    out.write((char *)mapping, mesh.numVertices() * sizeof(uint));
#endif
    delete[] mapping;


    uint num_ghost = mesh.distdata().num_ghost(0);
#ifdef ENABLE_MPIIO
#else
    out.write((char *)&num_ghost, sizeof(uint));
#endif
    uint *ghosts = new uint[2 * num_ghost];
    uint *gp = &ghosts[0];
    for (MeshGhostIterator g(mesh.distdata(), 0); !g.end(); ++g)
    {
      *gp++ = g.index();
      *gp++ = g.owner();
    }
#ifdef ENABLE_MPIIO
#else
    out.write((char *)ghosts, 2 * num_ghost * sizeof(uint));
#endif
    delete[] ghosts;

    uint num_shared = mesh.distdata().num_shared(0);
#ifdef ENABLE_MPIIO
#else
    out.write((char *)&num_shared, sizeof(uint));
#endif

    uint *shared = new uint[num_shared];
    uint *sp = &shared[0];
    for (MeshSharedIterator s(mesh.distdata(), 0); !s.end(); ++s)
      *sp++ = s.index();
#ifdef ENABLE_MPIIO
#else
    out.write((char *)shared, mesh.distdata().num_shared(0) * sizeof(uint));
#endif
    delete[] shared;

  }

}
//-----------------------------------------------------------------------------
void Checkpoint::write(std::vector<Function *> func, chkp_outstream& out)
{
  std::vector<Function *>::iterator it;

  uint max_size = 0;
  for (it = func.begin(); it != func.end(); ++it)
  {
    if((*it)->type() != Function::discrete)
      error("Checkpoint restart only implemented for discrete functions");    
    max_size = std::max(max_size, (*it)->vector().local_size());
  }
  
  real *values = new real[max_size];
#ifdef ENABLE_MPIIO
  uint vector_offset[2];
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();
  MPI_Offset tmp_offset;  
#endif

  for (it = func.begin(); it != func.end(); ++it)
  {
    uint local_size = (*it)->vector().local_size();
    (*it)->vector().get(values);

#ifdef ENABLE_MPIIO
    vector_offset[0] = (*it)->vector().offset();
    vector_offset[1] = (*it)->vector().local_size();
    
    tmp_offset = byte_offset + pe_rank * 2 * sizeof(uint);
    MPI_File_write_at_all(out, tmp_offset, &vector_offset[0], 2, 
			  MPI_UNSIGNED, MPI_STATUS_IGNORE);
    tmp_offset = byte_offset + pe_size * 2 * sizeof(uint) + vector_offset[0] * sizeof(real);

    MPI_File_write_at_all(out, tmp_offset, values, vector_offset[1],
			  MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += pe_size * 2 * sizeof(uint) + (*it)->vector().size() * sizeof(real);
#else
    out.write((char *)&local_size, sizeof(uint));
    out.write((char *)values, (*it)->vector().local_size() * sizeof(real));
#endif
  }
  delete[] values;
  
}
//-----------------------------------------------------------------------------
void Checkpoint::write(std::vector<Vector *> vec, chkp_outstream& out)
{
  std::vector<Vector *>::iterator it;

  uint max_size = 0;
  for (it = vec.begin(); it != vec.end(); ++it)
    max_size = std::max(max_size, (*it)->local_size());
  
  real *values = new real[max_size];
#ifdef ENABLE_MPIIO
  uint vector_offset[2];
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();
  MPI_Offset tmp_offset;  
#endif
  for (it = vec.begin(); it != vec.end(); ++it)
  {

    (*it)->get(values);

#ifdef ENABLE_MPIIO
    vector_offset[0] = (*it)->offset();
    vector_offset[1] = (*it)->local_size();
    
    tmp_offset = byte_offset + pe_rank * 2 * sizeof(uint);
    MPI_File_write_at_all(out, tmp_offset, &vector_offset[0], 2, 
			  MPI_UNSIGNED, MPI_STATUS_IGNORE);
    tmp_offset = byte_offset + pe_size * 2 * sizeof(uint) + vector_offset[0] * sizeof(real);

    MPI_File_write_at_all(out, tmp_offset, values, vector_offset[1],
			  MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += pe_size * 2 * sizeof(uint) + (*it)->size() * sizeof(real);
#else
    uint local_size = (*it)->local_size();
    out.write((char *)&local_size, sizeof(uint));
    out.write((char *)values, (*it)->local_size() * sizeof(real));
#endif
  }
  delete[] values;
  
}
//-----------------------------------------------------------------------------

