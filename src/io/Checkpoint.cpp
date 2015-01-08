// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-09
// Last changed: 2013-01-07

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

namespace dolfin
{

//-----------------------------------------------------------------------------
Checkpoint::Checkpoint() :
    state_(CHECKPOINT),
    restart_state_(OPEN),
#ifdef ENABLE_MPIIO
    byte_offset_(0),
#endif
    in_(),
    n_(0),
    id_(0),
    t_(0.0),
    hdr_initialized_(false),
    disp_initialized_(false)
{
}
//-----------------------------------------------------------------------------
Checkpoint::~Checkpoint()
{
}
//-----------------------------------------------------------------------------
void Checkpoint::hdr_init(Mesh& mesh, bool static_mesh)
{

  if (!hdr_initialized_ || !static_mesh)
  {
    hdr_.num_coords = mesh.numVertices() * mesh.geometry().dim();
    hdr_.num_entities = mesh.type().numEntities(0);
    hdr_.num_centities = mesh.numCells() * hdr_.num_entities;
    hdr_.type = mesh.type().cellType();
    hdr_.tdim = mesh.topology().dim();
    hdr_.gdim = mesh.geometry().dim();
    hdr_.num_vertices = mesh.numVertices();
    hdr_.num_cells = mesh.numCells();
    hdr_.num_ghosts = mesh.distdata().num_ghost(0);
    hdr_.num_shared = mesh.distdata().num_shared(0);

#ifdef ENABLE_MPIIO
    uint local_data[5] = { hdr_.num_coords, hdr_.num_centities, hdr_.num_vertices,
                           (2 * hdr_.num_ghosts), hdr_.num_shared };

    memset(&hdr_.offsets[0], 0, 5 * sizeof(uint));
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&local_data[0], &hdr_.offsets[0], 5, MPI_UNSIGNED, MPI_SUM,
               MPI::DOLFIN_COMM);
#else
    MPI_Scan(&local_data[0], &hdr_.offsets[0], 5,
        MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);

    hdr_.offsets[0] -= local_data[0];
    hdr_.offsets[1] -= local_data[1];
    hdr_.offsets[2] -= local_data[2];
    hdr_.offsets[3] -= local_data[3];
    hdr_.offsets[4] -= local_data[4];
#endif

    if (!disp_initialized_ || !static_mesh)
    {
      memset(&hdr_.disp[0], 0, 5 * sizeof(uint));
      MPI_Allreduce(&local_data[0], &hdr_.disp[0], 5, MPI_UNSIGNED, MPI_SUM,
                    MPI::DOLFIN_COMM);
      disp_initialized_ = true;
    }

#endif
  }

  hdr_initialized_ = true;

}
//-----------------------------------------------------------------------------
void Checkpoint::write(std::string fname, uint id, real t, Mesh& mesh,
                       std::vector<Function *> func, std::vector<Vector *> vec,
                       bool static_mesh)
{

  message("Writing checkpoint (%s%d) at time %g", fname.c_str(), n_ % 2, t);
  std::ostringstream _fname;
#ifndef ENABLE_MPIIO
  if( MPI::numProcesses() > 1)
  _fname << fname << (n_++)%2 << "_" << MPI::processNumber() << ".chkp";
  else
  _fname << fname << (n_++)%2 << ".chkp";

  std::ofstream out(_fname.str().c_str(), std::ofstream::binary);

  out.write((char *) &id, sizeof(uint));
  out.write((char *) &t, sizeof(real));

#else
  _fname << fname << (n_++) % 2 << ".chkp";

  MPI_File out;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) _fname.str().c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out);

  byte_offset_ = 0;
  MPI_File_write_all(out, &id, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  MPI_File_write_all(out, &t, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

  byte_offset_ += sizeof(uint);
  byte_offset_ += sizeof(real);
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

}
//-----------------------------------------------------------------------------
void Checkpoint::restart(std::string fname)
{

  if (restart_state_ != OPEN)
  {
    error("Shut her down, Scotty, she's sucking mud again!");
  }

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
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) _fname.str().c_str(),
                MPI_MODE_RDONLY, MPI_INFO_NULL, &in_);
  MPI_File_read_all(in_, &id_, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  MPI_File_read_all(in_, &t_, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
  byte_offset_ = sizeof(uint) + sizeof(real);
#else
  in_.open(_fname.str().c_str(), std::ifstream::binary);
  in_.read((char *) &id_, sizeof(uint));
  in_.read((char *) &t_, sizeof(real));
#endif
  message("Restarting from time %g checkpoint id %d", t_, id_);

  state_ = RESTART;
  restart_state_ = MESH;

}
//-----------------------------------------------------------------------------
void Checkpoint::load(Mesh& mesh)
{
  if (restart_state_ != MESH)
  {
    error("Shut her down, Scotty, she's sucking mud again!");
  }

#ifdef ENABLE_MPIIO
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();

  MPI_File_read_at_all(in_, byte_offset_ + pe_rank * sizeof(chkp_mesh_hdr), &hdr_,
                       sizeof(chkp_mesh_hdr), MPI_BYTE, MPI_STATUS_IGNORE);
  byte_offset_ += pe_size * sizeof(chkp_mesh_hdr);
#else
  in_.read((char *)&hdr_, sizeof(chkp_mesh_hdr));
#endif
  real * coords = new real[hdr_.num_coords];

#ifdef ENABLE_MPIIO
  MPI_File_read_at_all(in_, byte_offset_ + hdr_.offsets[0] * sizeof(real), coords,
                       hdr_.num_coords, MPI_DOUBLE, MPI_STATUS_IGNORE);
  byte_offset_ += hdr_.disp[0] * sizeof(real);
#else
  in_.read((char *)coords, (hdr_.num_coords) * sizeof(real));
#endif

  Mesh _mesh;
  MeshEditor editor(_mesh, hdr_.type, hdr_.tdim, hdr_.gdim);
  editor.initVertices(hdr_.num_vertices);

  uint vi = 0;
  for (uint i = 0; i < hdr_.num_coords; i += hdr_.gdim)
  {
    editor.addVertex(vi++, &coords[i]);
  }

  delete[] coords;

  editor.initCells(hdr_.num_cells);

  uint *cells = new uint[hdr_.num_centities];
#ifdef ENABLE_MPIIO
  MPI_File_read_at_all(in_, byte_offset_ + hdr_.offsets[1] * sizeof(uint), cells,
                       hdr_.num_centities, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  byte_offset_ += hdr_.disp[1] * sizeof(uint);
#else
  in_.read((char *)cells, (hdr_.num_centities) * sizeof(uint));
#endif

  Array<uint> v;
  uint ci = 0;
  for (uint i = 0; i < hdr_.num_centities; i += hdr_.num_entities)
  {
    v.clear();
    for (uint j = 0; j < hdr_.num_entities; ++j)
    {
      v.push_back(cells[i + j]);
    }
    editor.addCell(ci++, v);
  }
  editor.close();
  delete[] cells;

  if (MPI::numProcesses() > 1)
  {
    uint *mapping = new uint[_mesh.numVertices()];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in_, byte_offset_ + hdr_.offsets[2] * sizeof(uint),
                         mapping, hdr_.num_vertices, MPI_UNSIGNED,
                         MPI_STATUS_IGNORE);
    byte_offset_ += hdr_.disp[2] * sizeof(uint);
#else
    in_.read((char *)mapping, hdr_.num_vertices * sizeof(uint));
#endif
    for (VertexIterator v(_mesh); !v.end(); ++v)
      _mesh.distdata().set_map(v->index(), mapping[v->index()], 0);
    delete[] mapping;

    uint *ghosts = new uint[2 * hdr_.num_ghosts];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in_, byte_offset_ + hdr_.offsets[3] * sizeof(uint),
                         ghosts, 2 * hdr_.num_ghosts, MPI_UNSIGNED,
                         MPI_STATUS_IGNORE);
    byte_offset_ += hdr_.disp[3] * sizeof(uint);
#else
    in_.read((char *)ghosts, 2*hdr_.num_ghosts * sizeof(uint));
#endif
    for (uint i = 0; i < 2 * hdr_.num_ghosts; i += 2)
    {
      _mesh.distdata().set_ghost(ghosts[i], 0);
      _mesh.distdata().set_ghost_owner(ghosts[i], ghosts[i + 1], 0);
    }
    delete[] ghosts;

    uint *shared = new uint[hdr_.num_shared];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in_, byte_offset_ + hdr_.offsets[4] * sizeof(uint),
                         shared, hdr_.num_shared, MPI_UNSIGNED,
                         MPI_STATUS_IGNORE);
    byte_offset_ += hdr_.disp[4] * sizeof(uint);
#else
    in_.read((char *)shared, hdr_.num_shared * sizeof(uint));
#endif
    for (uint i = 0; i < hdr_.num_shared; i++)
    {
      _mesh.distdata().set_shared(shared[i], 0);
    }
    delete[] shared;
  }

  mesh = _mesh;
  mesh.distdata().set_invalid_numbering();
  mesh.renumber();

  restart_state_ = FUNC;

}
//-----------------------------------------------------------------------------
void Checkpoint::load(std::vector<Function *> func)
{
  if (restart_state_ != FUNC)
  {
    error("Shut her down, Scotty, she's sucking mud again!");
  }

  std::vector<Function *>::iterator it;
  uint local_size;
  uint pe_rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint vector_offset[3];
  // FIXME store max(local_size)
  for (it = func.begin(); it != func.end(); ++it)
  {
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in_, byte_offset_ + pe_rank * 3 * sizeof(uint),
                         &vector_offset[0], 3, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset_ += pe_size * 3 * sizeof(uint);
    local_size = vector_offset[1];
#else
    in_.read((char *)&local_size, sizeof(uint));
#endif
    real *values = new real[local_size];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in_, byte_offset_ + vector_offset[0] * sizeof(real),
                         values, vector_offset[1], MPI_DOUBLE,
                         MPI_STATUS_IGNORE);
    byte_offset_ += vector_offset[2] * sizeof(real);
#else
    in_.read((char *)values, local_size * sizeof(real));
#endif

    (*it)->vector().set(values);
    (*it)->vector().apply();

    delete[] values;
  }

  restart_state_ = VEC;
}
//-----------------------------------------------------------------------------
void Checkpoint::load(std::vector<Vector *> vec)
{
  if (restart_state_ != VEC)
  {
    error("Shut her down, Scotty, she's sucking mud again!");
  }

  std::vector<Vector *>::iterator it;
  uint local_size;
  uint pe_rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint vector_offset[3];
  for (it = vec.begin(); it != vec.end(); ++it)
  {
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in_, byte_offset_ + pe_rank * 3 * sizeof(uint),
                         &vector_offset[0], 3, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset_ += pe_size * 3 * sizeof(uint);
    local_size = vector_offset[1];
#else
    in_.read((char *)&local_size, sizeof(uint));
#endif
    real *values = new real[local_size];
#ifdef ENABLE_MPIIO
    MPI_File_read_at_all(in_, byte_offset_ + vector_offset[0] * sizeof(real),
                         values, vector_offset[1], MPI_DOUBLE,
                         MPI_STATUS_IGNORE);
    byte_offset_ += vector_offset[2] * sizeof(double);
#else
    in_.read((char *)values, local_size * sizeof(real));
#endif
    (*it)->set(values);
    (*it)->apply();
    delete[] values;
  }

#ifdef ENABLE_MPIIO
  MPI_File_close(&in_);
#else
  in_.close();
#endif
}
//-----------------------------------------------------------------------------
void Checkpoint::write(Mesh& mesh, chkp_outstream& out)
{

#ifdef ENABLE_MPIIO
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();

  MPI_File_write_at_all(out, byte_offset_ + pe_rank * sizeof(chkp_mesh_hdr),
                        &hdr_, sizeof(chkp_mesh_hdr), MPI_BYTE,
                        MPI_STATUS_IGNORE);
  byte_offset_ += pe_size * sizeof(chkp_mesh_hdr);

  MPI_File_write_at_all(out, byte_offset_ + hdr_.offsets[0] * sizeof(real),
                        mesh.geometry().coordinates(), hdr_.num_coords,
                        MPI_DOUBLE, MPI_STATUS_IGNORE);
  byte_offset_ += hdr_.disp[0] * sizeof(real);

  MPI_File_write_at_all(out, byte_offset_ + hdr_.offsets[1] * sizeof(uint),
                        mesh.cells(), hdr_.num_centities, MPI_UNSIGNED,
                        MPI_STATUS_IGNORE);
  byte_offset_ += hdr_.disp[1] * sizeof(uint);

#else
  out.write((char *)&hdr_, sizeof(chkp_mesh_hdr));
  out.write((char *)mesh.geometry().coordinates(), 
	    hdr_.num_coords * sizeof(real));
  out.write((char *)mesh.cells(), hdr_.num_centities * sizeof(uint));
#endif

  if (MPI::numProcesses() > 1)
  {
    uint *mapping = new uint[mesh.numVertices()];
    for (VertexIterator v(mesh); !v.end(); ++v)
      mapping[v->index()] = mesh.distdata().get_global(*v);
#ifdef ENABLE_MPIIO
    MPI_File_write_at_all(out, byte_offset_ + hdr_.offsets[2] * sizeof(uint),
                          mapping, hdr_.num_vertices, MPI_UNSIGNED,
                          MPI_STATUS_IGNORE);
    byte_offset_ += hdr_.disp[2] * sizeof(uint);
#else
    out.write((char *)mapping, hdr_.num_vertices * sizeof(uint));
#endif
    delete[] mapping;

    uint *ghosts = new uint[2 * hdr_.num_ghosts];
    uint *gp = &ghosts[0];
    for (MeshGhostIterator g(mesh.distdata(), 0); !g.end(); ++g)
    {
      *gp++ = g.index();
      *gp++ = g.owner();
    }
#ifdef ENABLE_MPIIO
    MPI_File_write_at_all(out, byte_offset_ + hdr_.offsets[3] * sizeof(uint),
                          ghosts, 2 * hdr_.num_ghosts, MPI_UNSIGNED,
                          MPI_STATUS_IGNORE);
    byte_offset_ += hdr_.disp[3] * sizeof(uint);
#else
    out.write((char *)ghosts, 2 * hdr_.num_ghosts * sizeof(uint));
#endif
    delete[] ghosts;

    uint *shared = new uint[hdr_.num_shared];
    uint *sp = &shared[0];
    for (MeshSharedIterator s(mesh.distdata(), 0); !s.end(); ++s)
      *sp++ = s.index();
#ifdef ENABLE_MPIIO
    MPI_File_write_at_all(out, byte_offset_ + hdr_.offsets[4] * sizeof(uint),
                          shared, hdr_.num_shared, MPI_UNSIGNED,
                          MPI_STATUS_IGNORE);
    byte_offset_ += hdr_.disp[4] * sizeof(uint);
#else
    out.write((char *)shared, hdr_.num_shared * sizeof(uint));
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
    if ((*it)->type() != Function::discrete)
    {
      error("Checkpoint restart only implemented for discrete functions");
    }
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
    (*it)->vector().get(values);

#ifdef ENABLE_MPIIO
    vector_offset[0] = (*it)->vector().offset();
    vector_offset[1] = (*it)->vector().local_size();
    vector_offset[2] = (*it)->vector().size();

    MPI_File_write_at_all(out, byte_offset_ + pe_rank * 3 * sizeof(uint),
                          &vector_offset[0], 3, MPI_UNSIGNED,
                          MPI_STATUS_IGNORE);
    byte_offset_ += pe_size * 3 * sizeof(uint);

    MPI_File_write_at_all(out, byte_offset_ + vector_offset[0] * sizeof(real),
                          values, vector_offset[1], MPI_DOUBLE,
                          MPI_STATUS_IGNORE);
    byte_offset_ += (*it)->vector().size() * sizeof(real);
#else
    uint local_size = (*it)->vector().local_size();
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
  {
    max_size = std::max(max_size, (*it)->local_size());
  }

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

    MPI_File_write_at_all(out, byte_offset_ + pe_rank * 3 * sizeof(uint),
                          &vector_offset[0], 3, MPI_UNSIGNED,
                          MPI_STATUS_IGNORE);
    byte_offset_ += pe_size * 3 * sizeof(uint);

    MPI_File_write_at_all(out, byte_offset_ + vector_offset[0] * sizeof(real),
                          values, vector_offset[1], MPI_DOUBLE,
                          MPI_STATUS_IGNORE);
    byte_offset_ += vector_offset[2] * sizeof(real);
#else
    uint local_size = (*it)->local_size();
    out.write((char *)&local_size, sizeof(uint));
    out.write((char *)values, (*it)->local_size() * sizeof(real));
#endif
  }
  delete[] values;

}
//-----------------------------------------------------------------------------

}
