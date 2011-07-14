// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-09
// Last changed: 2011-06-27

#include <cstring>
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
Checkpoint::Checkpoint() : state(CHECKPOINT), restart_state(OPEN), n(0), 
			   hdr_initialized(false), disp_initialized(false)
{
}
//-----------------------------------------------------------------------------
Checkpoint::~Checkpoint()
{
}
//-----------------------------------------------------------------------------
void Checkpoint::hdr_init(Mesh& mesh, bool static_mesh)
{

  if(!hdr_initialized || !static_mesh)
  {    
    hdr.num_coords = mesh.numVertices() * mesh.geometry().dim();
    hdr.num_entities = mesh.type().numEntities(0);
    hdr.num_centities = mesh.numCells() * hdr.num_entities;
    hdr.type = mesh.type().cellType();
    hdr.tdim = mesh.topology().dim();
    hdr.gdim = mesh.geometry().dim();
    hdr.num_vertices = mesh.numVertices();
    hdr.num_cells = mesh.numCells();  
    hdr.num_ghosts = mesh.distdata().num_ghost(0);
    hdr.num_shared = mesh.distdata().num_shared(0);
  
    
#ifdef ENABLE_MPIIO
    uint local_data[5] = {hdr.num_coords, 
			  hdr.num_centities,
			  hdr.num_vertices,
			  (2*hdr.num_ghosts),
			  hdr.num_shared};

    memset(&hdr.offsets[0], 0, 5 * sizeof(uint));
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&local_data[0], &hdr.offsets[0], 5,
	       MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
    MPI_Scan(&local_data[0], &hdr.offsets[0], 5, 
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
    
    hdr.offsets[0] -= local_data[0];
    hdr.offsets[1] -= local_data[1];
    hdr.offsets[2] -= local_data[2];
    hdr.offsets[3] -= local_data[3];
    hdr.offsets[4] -= local_data[4];
#endif
    
    if (!disp_initialized || !static_mesh) 
    {
      memset(&hdr.disp[0], 0, 5 * sizeof(uint));
      MPI_Allreduce(&local_data[0], &hdr.disp[0], 5, MPI_UNSIGNED,
		    MPI_SUM, MPI::DOLFIN_COMM);
      disp_initialized = true;
    }

#endif
  }
  
  hdr_initialized = true;

}
//-----------------------------------------------------------------------------
void Checkpoint::write(std::string fname, uint id, real t, Mesh& mesh, 
		       std::vector<Function *> func,
		       std::vector<Vector *> vec, bool static_mesh)
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

  hdr_init(mesh, static_mesh);
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
#ifndef ENABLE_MPIIO
  if( MPI::numProcesses() > 1) 
    _fname << fname << "_" << MPI::processNumber() << ".chkp";
  else
    _fname << fname << ".chkp";
#else
  _fname << fname << ".chkp";
#endif

#ifdef ENABLE_MPIIO
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *)  _fname.str().c_str(),
		MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
  MPI_File_read_all(in, &_id, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  MPI_File_read_all(in, &_t, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
  byte_offset = sizeof(uint) + sizeof(real);
#else
  in.open(_fname.str().c_str(), std::ifstream::binary);
  in.read((char *) &_id, sizeof(uint));
  in.read((char *) &_t, sizeof(real));
#endif
  message("Restarting from time %g checkpoint id %d", _t, _id);

  state = RESTART;
  restart_state = MESH;
      
}
//-----------------------------------------------------------------------------
void Checkpoint::load(Mesh& mesh)
{ 
  if (restart_state != MESH)
    error("Shut her down, Scotty, she's sucking mud again!");

#ifdef ENABLE_MPIIO
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();

  MPI_File_read_at_all(in, byte_offset + pe_rank * sizeof(chkp_mesh_hdr),
		       &hdr, sizeof(chkp_mesh_hdr), MPI_BYTE, MPI_STATUS_IGNORE);
  byte_offset += pe_size * sizeof(chkp_mesh_hdr);
#else
  in.read((char *)&hdr, sizeof(chkp_mesh_hdr));
#endif  
  real *coords = new real[hdr.num_coords];

#ifdef ENABLE_MPIIO
  MPI_File_read_at_all(in, byte_offset + hdr.offsets[0] * sizeof(real),
		       coords, hdr.num_coords, MPI_DOUBLE, MPI_STATUS_IGNORE);
  byte_offset += hdr.disp[0] * sizeof(real);
#else
  in.read((char *)coords, (hdr.num_coords) * sizeof(real));
#endif

  Mesh _mesh;
  MeshEditor editor;  
  editor.open(_mesh, hdr.type, hdr.tdim, hdr.gdim);
  editor.initVertices(hdr.num_vertices);

  uint vi = 0;
  for(uint i = 0 ; i < hdr.num_coords; i += hdr.gdim)
  {
    switch(hdr.gdim)
    {      
    case 2:
      editor.addVertex(vi++, coords[i], coords[i+1]); break;
    case 3:
      editor.addVertex(vi++, coords[i], coords[i+1], coords[i+2]); break;
    }
  }
  
  delete[] coords;

  editor.initCells(hdr.num_cells);

  uint *cells = new uint[hdr.num_centities];
#ifdef ENABLE_MPIIO
  MPI_File_read_at_all(in, byte_offset + hdr.offsets[1] * sizeof(uint), 
		       cells, hdr.num_centities, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  byte_offset += hdr.disp[1] * sizeof(uint);
#else
  in.read((char *)cells, (hdr.num_centities) * sizeof(uint));
#endif

  Array<uint> v;
  uint ci = 0;
  for (uint i = 0; i < hdr.num_centities; i += hdr.num_entities)
  {
    v.clear();
    for (uint j = 0; j < hdr.num_entities; j++)
      v.push_back(cells[i + j]);
    editor.addCell(ci++, v);    
  }  
  editor.close();
  delete[] cells;

  if (MPI::numProcesses() > 1) 
  {
    uint *mapping = new uint[_mesh.numVertices()];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in, byte_offset + hdr.offsets[2] * sizeof(uint),
			 mapping, hdr.num_vertices, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += hdr.disp[2] * sizeof(uint);
#else
    in.read((char *)mapping, hdr.num_vertices * sizeof(uint));
#endif
    for (VertexIterator v(_mesh); !v.end(); ++v)
      _mesh.distdata().set_map(v->index(), mapping[v->index()], 0);
    delete[] mapping;

    uint *ghosts = new uint[2 * hdr.num_ghosts];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in, byte_offset + hdr.offsets[3] * sizeof(uint),
			 ghosts, 2*hdr.num_ghosts, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += hdr.disp[3] * sizeof(uint);
#else
    in.read((char *)ghosts, 2*hdr.num_ghosts * sizeof(uint));
#endif
    for (uint i = 0; i < 2 * hdr.num_ghosts; i += 2)
    {
      _mesh.distdata().set_ghost(ghosts[i], 0);
      _mesh.distdata().set_ghost_owner(ghosts[i], ghosts[i+1], 0);
    }
    delete[] ghosts;
    
    uint *shared = new uint[hdr.num_shared];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in, byte_offset + hdr.offsets[4] * sizeof(uint),
			 shared, hdr.num_shared, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += hdr.disp[4] * sizeof(uint);
#else
    in.read((char *)shared, hdr.num_shared * sizeof(uint));
#endif
    for(uint i = 0; i < hdr.num_shared; i++)
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
  uint pe_rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint vector_offset[3];
  // FIXME store max(local_size) 
  for (it = func.begin(); it != func.end(); ++it)
  {
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in, byte_offset + pe_rank * 3 * sizeof(uint),
			 &vector_offset[0], 3, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += pe_size * 3 * sizeof(uint);
    local_size = vector_offset[1];   
#else
    in.read((char *)&local_size, sizeof(uint));
#endif
    real *values = new real[local_size];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in, byte_offset + vector_offset[0] * sizeof(real),
			 values, vector_offset[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += vector_offset[2] * sizeof(real);
#else
    in.read((char *)values, local_size * sizeof(real));    
#endif

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
  uint pe_rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint vector_offset[3];
  for (it = vec.begin(); it != vec.end(); ++it)
  {
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in, byte_offset + pe_rank * 3 * sizeof(uint),
			 &vector_offset[0], 3, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += pe_size * 3 * sizeof(uint);
    local_size = vector_offset[1];
#else
    in.read((char *)&local_size, sizeof(uint));
#endif
    real *values = new real[local_size];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in, byte_offset + vector_offset[0] * sizeof(real),
			 values, vector_offset[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += vector_offset[2] * sizeof(double);
#else
    in.read((char *)values, local_size * sizeof(real));
#endif
    (*it)->set(values);
    (*it)->apply();
    delete[] values;
  }
  
#ifdef ENABLE_MPIIO
  MPI_File_close(&in);
#else
  in.close();
#endif
}
//-----------------------------------------------------------------------------
void Checkpoint::write(Mesh& mesh, chkp_outstream& out)
{

#ifdef ENABLE_MPIIO
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();
  
  MPI_File_write_at_all(out, byte_offset + pe_rank * sizeof(chkp_mesh_hdr),
			&hdr, sizeof(chkp_mesh_hdr), MPI_BYTE, MPI_STATUS_IGNORE);
  byte_offset += pe_size * sizeof(chkp_mesh_hdr);
  
  MPI_File_write_at_all(out, byte_offset + hdr.offsets[0] * sizeof(real),
			mesh.coordinates(), hdr.num_coords, MPI_DOUBLE, MPI_STATUS_IGNORE);
  byte_offset += hdr.disp[0] * sizeof(real);
  
  MPI_File_write_at_all(out, byte_offset + hdr.offsets[1] * sizeof(uint),
			mesh.cells(), hdr.num_centities, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  byte_offset += hdr.disp[1] * sizeof(uint);

#else
  out.write((char *)&hdr, sizeof(chkp_mesh_hdr));
  out.write((char *)mesh.coordinates(), hdr.num_coords * sizeof(real));
  out.write((char *)mesh.cells(), hdr.num_centities * sizeof(uint));
#endif

  if (MPI::numProcesses() > 1) 
  {
    uint *mapping = new uint[mesh.numVertices()];
    for (VertexIterator v(mesh); !v.end(); ++v)
      mapping[v->index()] = mesh.distdata().get_global(*v);    
#ifdef ENABLE_MPIIO
    MPI_File_write_at_all(out, byte_offset + hdr.offsets[2] * sizeof(uint),
			  mapping, hdr.num_vertices, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += hdr.disp[2] * sizeof(uint);   
#else
    out.write((char *)mapping, hdr.num_vertices * sizeof(uint));
#endif
    delete[] mapping;

    uint *ghosts = new uint[2 * hdr.num_ghosts];
    uint *gp = &ghosts[0];
    for (MeshGhostIterator g(mesh.distdata(), 0); !g.end(); ++g)
    {
      *gp++ = g.index();
      *gp++ = g.owner();
    }
#ifdef ENABLE_MPIIO
    MPI_File_write_at_all(out, byte_offset + hdr.offsets[3] * sizeof(uint),
			  ghosts, 2 * hdr.num_ghosts, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += hdr.disp[3] * sizeof(uint);
#else
    out.write((char *)ghosts, 2 * hdr.num_ghosts * sizeof(uint));
#endif
    delete[] ghosts;

    uint *shared = new uint[hdr.num_shared];
    uint *sp = &shared[0];
    for (MeshSharedIterator s(mesh.distdata(), 0); !s.end(); ++s)
      *sp++ = s.index();
#ifdef ENABLE_MPIIO
    MPI_File_write_at_all(out, byte_offset + hdr.offsets[4] * sizeof(uint),
			  shared, hdr.num_shared, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += hdr.disp[4] * sizeof(uint);
#else
    out.write((char *)shared, hdr.num_shared * sizeof(uint));
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
  uint vector_offset[3];
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();
#endif

  for (it = func.begin(); it != func.end(); ++it)
  {
    uint local_size = (*it)->vector().local_size();
    (*it)->vector().get(values);

#ifdef ENABLE_MPIIO
    vector_offset[0] = (*it)->vector().offset();
    vector_offset[1] = (*it)->vector().local_size();
    vector_offset[2] = (*it)->vector().size();
    
    MPI_File_write_at_all(out, byte_offset + pe_rank * 3 * sizeof(uint),
			  &vector_offset[0], 3, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += pe_size * 3 * sizeof(uint);

    MPI_File_write_at_all(out, byte_offset + vector_offset[0] * sizeof(real),
			  values, vector_offset[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += (*it)->vector().size() * sizeof(real);
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
  uint vector_offset[3];
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();
#endif
  for (it = vec.begin(); it != vec.end(); ++it)
  {

    (*it)->get(values);

#ifdef ENABLE_MPIIO
    vector_offset[0] = (*it)->offset();
    vector_offset[1] = (*it)->local_size();
    vector_offset[2] = (*it)->size();
    
    MPI_File_write_at_all(out, byte_offset + pe_rank * 3 * sizeof(uint), 
			  &vector_offset[0], 3, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += pe_size * 3 * sizeof(uint);

    MPI_File_write_at_all(out, byte_offset + vector_offset[0] * sizeof(real),
			  values, vector_offset[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += vector_offset[2] * sizeof(real);
#else
    uint local_size = (*it)->local_size();
    out.write((char *)&local_size, sizeof(uint));
    out.write((char *)values, (*it)->local_size() * sizeof(real));
#endif
  }
  delete[] values;
  
}
//-----------------------------------------------------------------------------
