// Copyright (C) 2009-2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2014.
//
// First  added: 2009
// Last changed: 2014-11-05

#include <dolfin/common/types.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/parameters.h>
#include <typeinfo>

#ifdef ENABLE_MPIIO
#include <mpi.h>
#endif

#include <algorithm>
#include <cstring>
#include <fstream>

namespace dolfin
{

//----------------------------------------------------------------------------
BinaryFile::BinaryFile(const std::string filename) :
    GenericFile(filename),
    t_(0),
    version_(2)
{
  type = "Binary";
}
//----------------------------------------------------------------------------
BinaryFile::BinaryFile(const std::string filename, real const& t) :
    GenericFile(filename),
    t_(&t),
    version_(2)
{
  type = "Binary";
}
//----------------------------------------------------------------------------
BinaryFile::~BinaryFile()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void BinaryFile::read()
{
  opened_read = true;
}
//-----------------------------------------------------------------------------
void BinaryFile::write()
{
  opened_write = true;
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(GenericVector& x)
{

  uint size;

#ifdef ENABLE_MPIIO
  uint offset[2];
  uint pe_rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  
  MPI_File fh;
  MPI_Offset byte_offset;
  BinaryFileHeader hdr;
  bool byteswap;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                    MPI_STATUS_IGNORE);

  byteswap = hdr_check(hdr, BINARY_VECTOR_DATA, pe_size);

  byte_offset = sizeof(BinaryFileHeader);
  MPI_File_read_at_all(fh, byte_offset + pe_rank * 2 * sizeof(uint), &offset[0],
                       2, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  if (byteswap)
  {
    offset[0] = bswap(offset[0]);
    offset[1] = bswap(offset[1]);
  }
  size = offset[1];
  byte_offset += pe_size * 2 * sizeof(uint);

#else
  std::ifstream fp(filename.c_str(), std::ifstream::binary);
  fp.read((char *)&size, sizeof(uint));
#endif

  real *values = new real[size];

#ifdef ENABLE_MPIIO
  MPI_File_read_at_all(fh, byte_offset + offset[0] * sizeof(real), values,
                       offset[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
  if (byteswap)
  {
    for(uint i = 0; i < offset[1]; ++i)
    {
      values[i] = bswap(values[i]);
    }
  }
  MPI_File_close(&fh);
#else
  fp.read((char *)values, size * sizeof(real));
  fp.close();
#endif

  x.init(size);
  x.set(values);
  delete[] values;
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(GenericVector& x)
{
  uint size = x.local_size();
  real *values = new real[size];

#ifdef ENABLE_MPIIO
  uint offset[2] = { 0, 0 };
  offset[0] = x.offset();
  offset[1] = size;
#endif

  x.get(values);

#ifdef ENABLE_MPIIO
  BinaryFileHeader hdr;
  uint pe_rank = MPI::processNumber();
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = MPI::numProcesses();
  hdr.type = BINARY_VECTOR_DATA;
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  MPI_File fh;
  MPI_Offset byte_offset;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  MPI_File_write_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                     MPI_STATUS_IGNORE);
  byte_offset = sizeof(BinaryFileHeader);
  MPI_File_write_at_all(fh, byte_offset + pe_rank * 2 * sizeof(uint),
                        &offset[0], 2, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  byte_offset += hdr.pe_size * 2 * sizeof(uint);
  MPI_File_write_at_all(fh, byte_offset + offset[0] * sizeof(real), values,
                        offset[1], MPI_DOUBLE, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
#else
  std::ofstream fp(filename.c_str(), std::ofstream::binary);
  fp.write((char *)&size, sizeof(uint));
  fp.write((char *)values, x.local_size() * sizeof(real));
  fp.close();
#endif

  delete[] values;

  message(1, "Saved vector to file %s in binary format.", filename.c_str());
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(Function & f)
{

#ifdef ENABLE_MPIIO
  uint pe_size = MPI::numProcesses();

  MPI_File fh;
  MPI_Offset byte_offset;
  BinaryFileHeader hdr;
  bool byteswap;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                    MPI_STATUS_IGNORE);

  byteswap = hdr_check(hdr, BINARY_FUNCTION_DATA, pe_size);

  byte_offset = sizeof(BinaryFileHeader);

  uint nfunc;
  MPI_File_read_at_all(fh, byte_offset, &nfunc, sizeof(uint), MPI_BYTE,
                       MPI_STATUS_IGNORE);
  if (byteswap) nfunc = bswap(nfunc);
  byte_offset += sizeof(uint);
  if (nfunc > 1)
  {
    warning("File contains %d functions, using first with matching dim.",
            nfunc);
  }

  BinaryFunctionHeader f_hdr;
  for (uint i = 0; i < nfunc; ++i)
  {
    MPI_File_read_at_all(fh, byte_offset, &f_hdr, sizeof(BinaryFunctionHeader),
                         MPI_BYTE, MPI_STATUS_IGNORE);
    if (byteswap) bswap_func_hdr(f_hdr);
    byte_offset += sizeof(BinaryFunctionHeader);

    /* Load function if dimension match */
    if (f_hdr.dim == f.value_size())
    {

      uint size = f.value_size() * f.mesh().distdata().num_owned(0);
      real *values = new real[size];
      MPI_File_read_at_all(fh, byte_offset + f.vector().offset() * sizeof(real),
                           values, size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      if (byteswap)
      {
        for (uint j = 0; j < size; ++j)
        {
          values[j] = bswap(values[j]);
        }
      }
      f.vector().set(values);
      delete[] values;

      MPI_File_close(&fh);

      return;
    }

    /* Otherwise continue searching*/
    byte_offset += f_hdr.dim * f.mesh().global_numVertices() * sizeof(real);
  }

  MPI_File_close(&fh);

  error("No matching functions found in binary file");
#else
  error("MPI I/O required for loading functions written in  binary");
#endif
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(std::vector<std::pair<Function*, std::string> >& f)
{
#ifdef ENABLE_MPIIO
  uint pe_size = MPI::numProcesses();

  MPI_File fh;
  MPI_Offset byte_offset;
  BinaryFileHeader hdr;
  bool byteswap;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                    MPI_STATUS_IGNORE);

  byteswap = hdr_check(hdr, BINARY_FUNCTION_DATA, pe_size);

  byte_offset = sizeof(BinaryFileHeader);

  uint nfunc;
  MPI_File_read_at_all(fh, byte_offset, &nfunc, sizeof(uint), MPI_BYTE,
                       MPI_STATUS_IGNORE);
  if (byteswap) nfunc = bswap(nfunc);
  byte_offset += sizeof(uint);

  if (nfunc != f.size())
  {
    error("Number of functions mismatch between set and file");
  }

  BinaryFunctionHeader f_hdr;
  for (std::vector<std::pair<Function*, std::string> >::iterator it = f.begin();
      it != f.end(); it++)
  {

    MPI_File_read_at_all(fh, byte_offset, &f_hdr, sizeof(BinaryFunctionHeader),
                         MPI_BYTE, MPI_STATUS_IGNORE);
    if (byteswap) bswap_func_hdr(f_hdr);
    byte_offset += sizeof(BinaryFunctionHeader);

    Function * u = it->first;

    if (f_hdr.dim != u->value_size())
    {
      error("Dimension of file and function set does not match");
    }

    uint size = u->value_size() * u->mesh().distdata().num_owned(0);
    real *values = new real[size];
    MPI_File_read_at_all(fh, byte_offset + u->vector().offset() * sizeof(real),
                         values, size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    if(byteswap)
    {
      for (uint i = 0; i < size; ++i)
      {
        values[i] = bswap(values[i]);
      }
    }
    u->vector().set(values);
    delete[] values;

    byte_offset += f_hdr.dim * u->mesh().global_numVertices() * sizeof(real);
  }

  MPI_File_close(&fh);

#else
  error("MPI I/O required for loading functions written in  binary");
#endif
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(Function & u)
{
  std::pair<Function*, std::string> f(&u, "U");
  std::vector<std::pair<Function*, std::string> > tmp;
  tmp.push_back(f);
  write_function(tmp);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(std::vector<std::pair<Function*, std::string> >& f)
{
  write_function(f);
}
//----------------------------------------------------------------------------
void BinaryFile::write_function(
    std::vector<std::pair<Function*, std::string> >& f)
{

#ifdef ENABLE_MPIIO

  nameUpdate(counter);

  BinaryFileHeader hdr;
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = MPI::numProcesses();
  hdr.type = BINARY_FUNCTION_DATA;
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  MPI_File fh;
  MPI_Offset byte_offset;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) bin_filename_.c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  MPI_File_write_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                     MPI_STATUS_IGNORE);
  byte_offset = sizeof(BinaryFileHeader);

  uint n_func = f.size();
  MPI_File_write_all(fh, &n_func, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  byte_offset += sizeof(uint);

  // Assume same mesh for all data arrays
  Mesh& mesh = f[0].first->mesh();

  for (std::vector<std::pair<Function*, std::string> >::iterator it = f.begin();
      it != f.end(); it++)
  {
    Function* u = it->first;

    std::string& name = it->second;

    // Get number of components
    uint const value_dim = u->value_size();

    // Allocate memory for function values at vertices
    uint size = mesh.numVertices() * u->value_size();
    real *values = new real[size];
    uint offset = u->vector().offset();

    if ((u->vector().local_size() / value_dim) != mesh.distdata().num_owned(0))
    {
      real *interp_values = new real[size];

      // Get function values at vertices
      u->interpolate_vertex_values(interp_values);

      uint ii = 0;
      for (VertexIterator v(mesh); !v.end(); ++v)
      {
        if (!v->is_ghost())
        {
          for (uint i = 0; i < value_dim; ++i)
          {
            values[ii++] = interp_values[v->index() + i * mesh.numVertices()];
          }
        }
      }
      delete[] interp_values;

      // Compute new vertex based offset
      uint num_values = value_dim * mesh.distdata().num_owned(0);
      offset = 0;
#if ( MPI_VERSION > 1 )
      MPI_Exscan(&num_values, &offset, 1, MPI_UNSIGNED, MPI_SUM,
                 MPI::DOLFIN_COMM);
#else
      MPI_Scan(&num_values, &offset, 1, MPI_UNSIGNED, MPI_SUM,
          MPI::DOLFIN_COMM);
      offset -= num_values;
#endif

    }
    else
    {
      u->vector().get(values);
    }

    size = value_dim * mesh.distdata().num_owned(0);

    BinaryFunctionHeader f_hdr;
    f_hdr.dim = value_dim;
    f_hdr.size = value_dim * mesh.global_numVertices();
    if (name.length() > FNAME_LENGTH) error("Function name too long.");
    strcpy(&f_hdr.name[0], name.c_str());
    if (t_)
    {
      f_hdr.t = *t_;
    }
    else f_hdr.t = counter;
    MPI_File_write_at_all(fh, byte_offset, &f_hdr, sizeof(BinaryFunctionHeader),
                          MPI_BYTE, MPI_STATUS_IGNORE);
    byte_offset += sizeof(BinaryFunctionHeader);

    MPI_File_write_at_all(fh, byte_offset + offset * sizeof(real), values, size,
                          MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += value_dim * (mesh.global_numVertices() * sizeof(real));

    delete[] values;
  }

  MPI_File_close(&fh);

  counter++;
#else
  error("MPI I/O is required to save functions in Binary.");
#endif
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(Mesh& mesh)
{

  bool byteswap;
  BinaryFileHeader hdr;
  uint pe_size = MPI::numProcesses();
  uint pe_rank = MPI::processNumber();

  if (MPI::numProcesses() == 1 || dolfin_get("Mesh read in serial"))
  {
    if(MPI::processNumber() > 0)
    {
      error("Reading serial mesh in parallel not implemented");
    }

    std::ifstream fp(filename.c_str(), std::ifstream::binary);

    uint type = 0;
    uint gdim = 0;
    fp.read((char *) &hdr, sizeof(BinaryFileHeader));
    byteswap = hdr_check(hdr, BINARY_MESH_DATA, pe_size);
    fp.read((char *) &gdim, sizeof(uint));
    fp.read((char *) &type, sizeof(uint));
    
    if (byteswap)
    {
      gdim = bswap(gdim);
      type = bswap(type);
    }

    // Create cell type to get topological dimension and number of vertices
    CellType::Type ctype = BinaryFile::cell_type(type);
    CellType* cell_type = CellType::create(ctype);
    uint const num_cellvertices = cell_type->num_entities(0);
    delete cell_type;

    // Open mesh for editing
    MeshEditor editor(mesh, ctype, gdim);

    // Read vertex data
    uint num_vertices = 0;
    fp.read((char *) &num_vertices, sizeof(uint));
    if (byteswap) num_vertices = bswap(num_vertices);
    editor.init_vertices(num_vertices);
    uint const vertex_data_size = num_vertices * gdim;
    real * vertex_data = new real[vertex_data_size];
    fp.read((char *) vertex_data, vertex_data_size * sizeof(real));
    for (uint v = 0; v < num_vertices; ++v)
    {
      if (byteswap)
      {
        for (uint gi = 0; gi < gdim; ++gi)
        {
          vertex_data[(v * gdim) + gi] = bswap(vertex_data[(v * gdim) + gi]);
        }
      }
      editor.add_vertex(v, &vertex_data[v * gdim]);
    }
    delete[] vertex_data;

    // Read cell data
    uint num_cells = 0;
    fp.read((char *) &num_cells, sizeof(uint));
    if(byteswap) num_cells = bswap(num_cells);
    editor.init_cells(num_cells);
    uint const cell_data_size = num_cells * num_cellvertices;
    uint * cell_data = new uint[cell_data_size];
    fp.read((char *) cell_data, cell_data_size * sizeof(uint));
    for (uint c = 0; c < num_cells; ++c)
    {
      if (byteswap)
      {
        for (uint vi = 0; vi < num_cellvertices; ++vi)
        {
          cell_data[(c * num_cellvertices) + vi] =
            bswap(cell_data[(c * num_cellvertices) + vi]);
        }
      }
      editor.add_cell(c, &cell_data[c * num_cellvertices]);
    }
    delete[] cell_data;
    editor.close();
    fp.close();
  }
  else
  {
#ifdef ENABLE_MPIIO

    MPI_File fh;
    MPI_Offset byte_offset;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    // Read binary header
    uint gdim, type, num_vertices;
    MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                      MPI_STATUS_IGNORE);
    byteswap = hdr_check(hdr, BINARY_MESH_DATA, pe_size);
    MPI_File_read_all(fh, &gdim, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    MPI_File_read_all(fh, &type, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    MPI_File_read_all(fh, &num_vertices, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    
    if (byteswap) 
    {
      gdim = bswap(gdim);
      type = bswap(type);
      num_vertices = bswap(num_vertices);
    }

    // Update offset: header + (gdim + type + num_vertices)
    byte_offset = sizeof(BinaryFileHeader) + 3 * sizeof(uint);

    // Create cell type to get topological dimension and number of vertices
    CellType::Type ctype = BinaryFile::cell_type(type);
    CellType * cell_type = CellType::create(ctype);
    uint const num_cellvertices = cell_type->num_entities(0);
    delete cell_type;

    //
    uint L = floor((real) num_vertices / (real) pe_size);
    uint R = num_vertices % pe_size;
    uint local_vertices = (num_vertices + pe_size - pe_rank - 1) / pe_size;

    uint vertex_offset[2] = { 0, 0 };
    uint vertex_data[2] = { local_vertices, gdim * local_vertices };
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&vertex_data[0], &vertex_offset[0], 2, MPI_UNSIGNED, MPI_SUM,
               MPI::DOLFIN_COMM);
#else
    MPI_Scan(&vertex_data[0], &vertex_offset[0], 2, MPI_UNSIGNED, MPI_SUM,
             MPI::DOLFIN_COMM);
    vertex_offset[0] -= vertex_data[0];
    vertex_offset[1] -= vertex_data[1];
#endif

    real * vertex_buffer = new real[vertex_data[1]];
    MPI_File_read_at_all(fh, byte_offset + vertex_offset[1] * sizeof(real),
                         vertex_buffer, vertex_data[1], MPI_DOUBLE,
                         MPI_STATUS_IGNORE);
    byte_offset += gdim * num_vertices * sizeof(real);

    if (byteswap)
    {
      for (uint v = 0; v < vertex_data[1]; ++v)
      {
        vertex_buffer[v] = bswap(vertex_buffer[v]);
      }
    }

    uint num_cells;
    MPI_File_read_at_all(fh, byte_offset, &num_cells, 1, MPI_UNSIGNED,
                         MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);
    if (byteswap) num_cells = bswap(num_cells);

    uint const num_local_cells = (num_cells + pe_size - pe_rank - 1) / pe_size;

    uint cell_offset = 0;
    uint cell_data = num_cellvertices * num_local_cells;
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&cell_data, &cell_offset, 1, MPI_UNSIGNED, MPI_SUM,
               MPI::DOLFIN_COMM);
#else
    MPI_Scan(&cell_data, &cell_offset, 1, MPI_UNSIGNED, MPI_SUM,
             MPI::DOLFIN_COMM);
    cell_offset -= cell_data;
#endif

    uint * cell_buffer = new uint[cell_data];
    MPI_File_read_at_all(fh, byte_offset + cell_offset * sizeof(uint),
                         cell_buffer, cell_data, MPI_UNSIGNED,
                         MPI_STATUS_IGNORE);

    if (byteswap)
    {
      for (uint c = 0; c < cell_data; ++c)
      {
        cell_buffer[c] = bswap(cell_buffer[c]);
      }
    }

    MPI_File_close(&fh);

    _set<uint> ghosted_entities;

    // Parse cells
    std::vector<atomic_cell> cells;
    std::vector<uint> * non_local_cells = new std::vector<uint>[pe_size];
    std::set<uint> owned_vertices;
    atomic_cell cell(num_cellvertices);
    for (uint i = 0; i < cell_data; i += num_cellvertices)
    {
      uint const owner = vertex_owner(L, R, cell_buffer[i]);
      if (owner == pe_rank)
      {
        cell.v[0] = cell_buffer[i];
        owned_vertices.insert(cell.v[0]);
        for (uint n = 1; n < num_cellvertices; ++n)
        {
          cell.v[n] = cell_buffer[i + n];
          if (vertex_owner(L, R, cell_buffer[i + n]) != pe_rank)
          {
            ghosted_entities.insert(cell.v[n]);
          }
          else
          {
            owned_vertices.insert(cell.v[n]);
          }
        }
        cells.push_back(cell);
      }
      else
      {
        for (uint n = 0; n < num_cellvertices; ++n)
        {
          non_local_cells[owner].push_back(cell_buffer[i + n]);
        }
      }
    }
    delete[] cell_buffer;

    /*
     * FIXME
     * Reduce communication in this section
     */
    uint local_max = 0;
    for (uint i = 0; i < pe_size; ++i)
    {
      local_max = std::max(local_max, (uint) non_local_cells[i].size());
    }

    uint buff_size = 0;
    MPI_Allreduce(&local_max, &buff_size, 1, MPI_UNSIGNED, MPI_MAX,
                  MPI::DOLFIN_COMM);
    uint *recv_buffer = new uint[buff_size];

    // Exchange data
    MPI_Status status;
    int num_recv;
    int src;
    int dest;
    for (uint i = 1; i < pe_size; ++i)
    {
      src = (pe_rank - i + pe_size) % pe_size;
      dest = (pe_rank + i) % pe_size;

      MPI_Sendrecv(&non_local_cells[dest][0], non_local_cells[dest].size(),
                   MPI_UNSIGNED, dest, 1, recv_buffer, buff_size, MPI_UNSIGNED,
                   src, 1, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &num_recv);

      // Add received cells
      for (int j = 0; j < num_recv; j += num_cellvertices)
      {
        std::memcpy(&cell.v[0], &recv_buffer[j],
                    num_cellvertices * sizeof(uint));
        owned_vertices.insert(cell.v[0]);
        cells.push_back(cell);
        for (uint n = 1; n < num_cellvertices; ++n)
        {
          if (vertex_owner(L, R, cell.v[n]) != pe_rank)
          {
            ghosted_entities.insert(cell.v[n]);
          }
          else
          {
            owned_vertices.insert(cell.v[n]);
          }
        }
      }
    }

    // Number of vertices in mesh, owned + ghosts
    std::set<uint> all_vertices;
    std::set<uint> orphaned_vertices;
    for (uint i = 0; i < local_vertices; ++i)
    {
      all_vertices.insert(vertex_offset[0] + i);
    }
    std::set_difference(
        all_vertices.begin(), all_vertices.end(), owned_vertices.begin(),
        owned_vertices.end(),
        std::inserter(orphaned_vertices, orphaned_vertices.end()));
    uint const num_local_vertices = all_vertices.size()
        + ghosted_entities.size();

    // Open mesh for editing
    MeshEditor editor(mesh, ctype, gdim);
    editor.init_vertices(num_local_vertices);
    editor.init_cells(cells.size());

    MeshDistributedData distdata(mesh.topology().dim());

    /* Add local indices
     * Use start_index for global number
     *
     * Exchange ghost vertices
     * Use ghosts to know global number
     *
     * Add cells, remember to map global->local
     */

    uint local_vertex_index = 0;
    uint v_index;
    for (std::set<uint>::iterator it = owned_vertices.begin();
        it != owned_vertices.end(); ++local_vertex_index, ++it)
    {
      v_index = *it - vertex_offset[0];
      distdata.set_map(local_vertex_index, *it, 0);
      editor.add_vertex(local_vertex_index, &vertex_buffer[(gdim * v_index)]);
    }

    _map<uint, uint> new_owner;
    // Exchange ghost points
    std::vector<uint> * ghosts = new std::vector<uint>[pe_size];
    for (_set<uint>::iterator it = ghosted_entities.begin();
    it != ghosted_entities.end(); it++)
    {
      ghosts[vertex_owner(L, R, *it)].push_back(*it);
    }

    local_max = 0;
    for (uint i = 0; i < pe_size; ++i)
    {
      local_max = std::max(local_max, (uint) ghosts[i].size());
    }

    buff_size = 0;
    MPI_Allreduce(&local_max, &buff_size, 1, MPI_UNSIGNED, MPI_MAX,
                  MPI::DOLFIN_COMM);

    delete[] recv_buffer;
    recv_buffer = new uint[buff_size];

    uint *send_new_owner = new uint[buff_size];
    uint *recv_new_owner = new uint[buff_size];

    real *send_buffer_coords = new real[buff_size * gdim];
    real *recv_buffer_coords = new real[buff_size * gdim];

    for (uint i = 1; i < pe_size; ++i)
    {
      src = (pe_rank - i + pe_size) % pe_size;
      dest = (pe_rank + i) % pe_size;

      MPI_Sendrecv(&ghosts[dest][0], ghosts[dest].size(), MPI_UNSIGNED, dest, 1,
                   recv_buffer, buff_size, MPI_UNSIGNED, src, 1,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &num_recv);

      /*
       * Check if orphaned
       * Send back global number and orphaned info
       */

      real *sp = &send_buffer_coords[0];
      uint *np = &send_new_owner[0];
      for (int k = 0; k < num_recv; ++k)
      {
        v_index = recv_buffer[k] - vertex_offset[0];

        if (orphaned_vertices.find(recv_buffer[k]) != orphaned_vertices.end())
        {
          if (new_owner.find(recv_buffer[k]) == new_owner.end())
          {
            new_owner[recv_buffer[k]] = src;
          }
          *(np++) = new_owner[recv_buffer[k]];
        }
        else
        {
          *(np++) = pe_rank;
        }

        for (uint l = 0; l < gdim; ++l)
        {
          *(sp++) = vertex_buffer[(gdim * v_index) + l];
        }
      }
      MPI_Sendrecv(send_new_owner, num_recv, MPI_UNSIGNED, src, 1,
                   recv_new_owner, buff_size, MPI_UNSIGNED, dest, 1,
                   MPI::DOLFIN_COMM, &status);

      MPI_Sendrecv(send_buffer_coords, (num_recv * gdim), MPI_DOUBLE, src, 1,
                   recv_buffer_coords, (buff_size * gdim), MPI_DOUBLE, dest, 1,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &num_recv);

      int g_i = 0;
      for (int k = 0; k < num_recv; local_vertex_index++, ++g_i, k += gdim)
      {
        distdata.set_map(local_vertex_index, ghosts[dest][g_i], 0);
        if (recv_new_owner[g_i] != pe_rank)
        {
          distdata.set_ghost(local_vertex_index, 0);
          distdata.set_ghost_owner(local_vertex_index, recv_new_owner[g_i], 0);
        }
        editor.add_vertex(local_vertex_index, &recv_buffer_coords[k]);
      }
    }

    uint local_cell_index = 0;
    uint * connectivity = new uint[num_cellvertices];
    for (std::vector<atomic_cell>::iterator it = cells.begin();
        it != cells.end(); ++local_cell_index, ++it)
    {
      for (uint n = 0; n < it->size; ++n)
      {
        connectivity[n] = distdata.get_vertex_local(it->v[n]);
      }
      editor.add_cell(local_cell_index, &connectivity[0]);
    }
    delete[] connectivity;

    editor.close();
    mesh.distdata() = distdata;

    delete[] recv_buffer;
    delete[] vertex_buffer;
    delete[] recv_buffer_coords;
    delete[] send_buffer_coords;
    delete[] recv_new_owner;
    delete[] send_new_owner;

    for (uint i = 0; i < pe_size; ++i)
    {
      non_local_cells[i].clear();
      ghosts[i].clear();
    }
    delete[] non_local_cells;
    delete[] ghosts;

#else
    error("MPI and MPI I/O is required for reading binary meshes in parallel");
#endif

  }
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(Mesh& mesh)
{
  uint const gdim = mesh.geometry().dim();
  uint const type = BinaryFile::cell_type(mesh.type().cellType());
  uint const num_vertices = mesh.global_numVertices();
  uint const num_cells = mesh.global_numCells();
  uint const num_cellvertices = mesh.type().num_entities(0);

  BinaryFileHeader hdr;
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = MPI::numProcesses();
  hdr.type = BINARY_MESH_DATA;
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  if (MPI::numProcesses() == 1 || ! mesh.is_distributed())
  {
    if(MPI::processNumber() > 0)
    {
      error("Writing serial mesh in parallel not implemented");
    }

    std::ofstream fp(filename.c_str(), std::ofstream::binary);

    // Write Header
    fp.write((char *) &hdr, sizeof(BinaryFileHeader));
    fp.write((char *) &gdim, sizeof(uint));
    fp.write((char *) &type, sizeof(uint));

    // Write vertices
    fp.write((char *) &num_vertices, sizeof(uint));
    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      fp.write((char *) v->x(), gdim * sizeof(real));
    }

    // Write cells
    fp.write((char *) &num_cells, sizeof(uint));
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      fp.write((char *) c->entities(0), num_cellvertices * sizeof(uint));
    }

    fp.close();
  }
  else
  {

#ifdef ENABLE_MPIIO

    MPI_File fh;
    MPI_Offset byte_offset;

    /* FIXME:
     * Add MPI_Info data
     * Split and cleanup implementation
     */

    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    // Write Header
    MPI_File_write_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                       MPI_STATUS_IGNORE);
    MPI_File_write_all(fh, (void *) &gdim, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    MPI_File_write_all(fh, (void *) &type, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset = sizeof(BinaryFileHeader) + 2 * sizeof(uint);

    // Write vertices
    uint vertex_offset = 0;
    uint vertex_buffer_size = gdim * mesh.distdata().num_owned(0);
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&vertex_buffer_size, &vertex_offset, 1, MPI_UNSIGNED, MPI_SUM,
               MPI::DOLFIN_COMM);
#else
    MPI_Scan(&vertex_buffer_size, &vertex_offset, 1, MPI_UNSIGNED, MPI_SUM,
             MPI::DOLFIN_COMM);
    vertex_offset -= vertex_buffer_size;
#endif
    real * vertex_buffer = new real[vertex_buffer_size];
    real * vptr = &vertex_buffer[0];
    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      if (!v->is_ghost())
      {
        std::memcpy(vptr, &v->x()[0], gdim * sizeof(real));
        vptr += gdim;
      }
    }
    MPI_File_write_all(fh, (void *) &num_vertices, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);
    MPI_File_write_at_all(fh, byte_offset + vertex_offset * sizeof(real),
                          vertex_buffer, vertex_buffer_size, MPI_DOUBLE,
                          MPI_STATUS_IGNORE);
    byte_offset += gdim * num_vertices * sizeof(real);
    delete[] vertex_buffer;

    // Write Cells
    uint cell_offset = 0;
    uint cell_buffer_size = num_cellvertices * mesh.numCells();
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&cell_buffer_size, &cell_offset, 1, MPI_UNSIGNED, MPI_SUM,
               MPI::DOLFIN_COMM);
#else
    MPI_Scan(&cell_buffer_size, &cell_offset, 1, MPI_UNSIGNED, MPI_SUM,
             MPI::DOLFIN_COMM);
    cell_offset -= cell_buffer_size;
#endif
    uint * cell_buffer = new uint[cell_buffer_size];
    uint * cp = &cell_buffer[0];
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      for (uint i = 0; i < c->num_entities(0); ++i)
      {
        *(cp++) = mesh.distdata().get_vertex_global(c->entities(0)[i]);
      }
    }
    MPI_File_write_at_all(fh, byte_offset, (void *) &num_cells, 1, MPI_UNSIGNED,
                          MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);
    MPI_File_write_at_all(fh, byte_offset + cell_offset * sizeof(uint),
                          cell_buffer, cell_buffer_size, MPI_UNSIGNED,
                          MPI_STATUS_IGNORE);

    delete[] cell_buffer;

    MPI_File_close(&fh);

#else
    error("MPI and MPI I/O is required for writing a mesh in parallel");
#endif

  }

  message(1, "Saved mesh to file %s in binary format.", filename.c_str());
}
//----------------------------------------------------------------------------
void BinaryFile::nameUpdate(const int counter)
{
  std::string filestart, extension;
  std::ostringstream fileid, newfilename;

  fileid.fill('0');
  fileid.width(6);

  filestart.assign(filename, 0, filename.find("."));
  extension.assign(filename, filename.find("."), filename.size());

  fileid << counter;
  newfilename << filestart << fileid.str() << ".bin";

  bin_filename_ = newfilename.str();
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(MeshFunction<bool>& meshfunction)
{
  write_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(MeshFunction<int>& meshfunction)
{
  write_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(MeshFunction<unsigned int>& meshfunction)
{
  write_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(MeshFunction<double>& meshfunction)
{
  write_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
template<typename T>
  void BinaryFile::write_meshfunction(MeshFunction<T>& meshfunction)
  {

#ifdef ENABLE_MPIIO
    nameUpdate(counter);

    Mesh& mesh = meshfunction.mesh();
    real *values = new real[meshfunction.size()];
    real *vp = &values[0];

    BinaryFileHeader hdr;
    hdr.magic = BINARY_MAGIC;
    hdr.pe_size = MPI::numProcesses();
    hdr.type = BINARY_MESH_FUNCTION_DATA;
#ifdef HAVE_BIG_ENDIAN
    hdr.bendian = 1;
#else
    hdr.bendian = 0;
#endif

    MPI_File fh;
    MPI_Offset byte_offset;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) bin_filename_.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    MPI_File_write_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                       MPI_STATUS_IGNORE);
    byte_offset = sizeof(BinaryFileHeader);

    uint local_size = 0;
    int mfunc_type = 0;
    if (meshfunction.dim() == mesh.topology().dim())
    {
      MPI_File_write_at_all(fh, byte_offset, &mfunc_type, 1, MPI_UNSIGNED,
                            MPI_STATUS_IGNORE);
      byte_offset += sizeof(uint);

      for (CellIterator c(mesh); !c.end(); ++c)
      {
        *(vp++) = (real) meshfunction.get(c->index());
      }

      local_size = mesh.numCells();
    }
    else if (meshfunction.dim() == 0)
    {
      mfunc_type = 1;
      MPI_File_write_at_all(fh, byte_offset, &mfunc_type, 1, MPI_UNSIGNED,
                            MPI_STATUS_IGNORE);
      byte_offset += sizeof(uint);

      for (VertexIterator v(mesh); !v.end(); ++v)
      {
        if (!mesh.distdata().is_ghost(v->index(), 0))
        {
          *(vp++) = (real) meshfunction.get(v->index());
        }
      }

      local_size = mesh.distdata().num_owned(0);
    }
    else
    {
      error("Binary output of mesh functions is implemented "
            "for cell/vertex-based functions only.");
    }

    uint offset = 0;
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED, MPI_SUM,
               MPI::DOLFIN_COMM);
#else
    MPI_Scan(&local_size, &offset, 1, MPI_UNSIGNED, MPI_SUM,
             MPI::DOLFIN_COMM);
    offset -= local_size;
#endif

    MPI_File_write_at_all(fh, byte_offset + offset * sizeof(real), values,
                          local_size * sizeof(real), MPI_BYTE,
                          MPI_STATUS_IGNORE);

    MPI_File_close(&fh);

    counter++;

    delete[] values;

#else
    error("MPI I/O required for writing mesh functions to binary files");
#endif

  }
//----------------------------------------------------------------------------
void BinaryFile::operator>>(MeshFunction<bool>& meshfunction)
{
  read_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(MeshFunction<int>& meshfunction)
{
  read_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(MeshFunction<unsigned int>& meshfunction)
{
  read_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(MeshFunction<double>& meshfunction)
{
  read_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
template<typename T>
  void BinaryFile::read_meshfunction(MeshFunction<T>& meshfunction)
  {

#ifdef ENABLE_MPIIO

    Mesh& mesh = meshfunction.mesh();
    real *values = new real[meshfunction.size()];

    uint pe_rank = MPI::processNumber();
    uint pe_size = MPI::numProcesses();
    message("load func %s", filename.c_str());
    MPI_File fh;
    MPI_Offset byte_offset;
    BinaryFileHeader hdr;
    bool byteswap;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
                      MPI_STATUS_IGNORE);

    byteswap = hdr_check(hdr, BINARY_MESH_FUNCTION_DATA, pe_size);

    byte_offset = sizeof(BinaryFileHeader);

    int mfunc_type = 0;
    MPI_File_read_at_all(fh, byte_offset, &mfunc_type, 1, MPI_UNSIGNED,
                         MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);
    if (byteswap) mfunc_type = bswap(mfunc_type);


    if ((mfunc_type == 0 && meshfunction.dim() != mesh.topology().dim())
        || (mfunc_type == 1 && meshfunction.dim() != 0))
    {
      error("Meshfunction does not match data in file");
    }

    uint local_size = (
        mfunc_type > 0 ?
            mesh.numVertices() - mesh.distdata().num_ghost(0) : mesh.numCells());

    uint offset = 0;
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED, MPI_SUM,
               MPI::DOLFIN_COMM);
#else
    MPI_Scan(&local_size, &offset, 1, MPI_UNSIGNED, MPI_SUM,
             MPI::DOLFIN_COMM);
    offset -= local_size;
#endif

    MPI_File_read_at_all(fh, byte_offset + offset * sizeof(real), values,
                         local_size * sizeof(real), MPI_BYTE,
                         MPI_STATUS_IGNORE);
    if(byteswap)
    {
      for (uint i = 0; i < (local_size * sizeof(real)); ++i)
      {
        values[i] = bswap(values[i]);
      }
    }
    if (mfunc_type == 0)
    {
      for (uint i = 0; i < meshfunction.size(); ++i)
      {
        meshfunction.set(i, static_cast<T>(values[i]));
      }
    }
    if (mfunc_type == 1)
    {
      for (uint i = 0; i < meshfunction.size(); ++i)
      {
        meshfunction.set(i, static_cast<T>(values[i]));
      }

      std::vector<uint> *ghost_buff = new std::vector<uint>[pe_size];
      for (MeshGhostIterator it(mesh.distdata(), 0); !it.end(); ++it)
        ghost_buff[it.owner()].push_back(
            mesh.distdata().get_vertex_global(it.index()));

      MPI_Status status;
      std::vector<real> send_buff;
      uint src, dest;
      uint recv_size = mesh.distdata().num_ghost(0);
      int recv_count, recv_size_gh, send_size;

      for (uint i = 0; i < pe_size; ++i)
      {
        send_size = ghost_buff[i].size();
        MPI_Reduce(&send_size, &recv_size_gh, 1, MPI_UNSIGNED, MPI_SUM, i,
                   MPI::DOLFIN_COMM);
      }
      uint *recv_ghost = new uint[recv_size_gh];
      real *recv_buff = new real[recv_size];

      for (uint j = 1; j < pe_size; j++)
      {
        src = (pe_rank - j + pe_size) % pe_size;
        dest = (pe_rank + j) % pe_size;

        MPI_Sendrecv(&ghost_buff[dest][0], ghost_buff[dest].size(),
                     MPI_UNSIGNED, dest, 1, recv_ghost, recv_size_gh,
                     MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

        for (int k = 0; k < recv_count; ++k)
          send_buff.push_back(
              meshfunction.get(
                  mesh.distdata().get_vertex_local(recv_ghost[k])));

        MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_DOUBLE, src, 2,
                     recv_buff, recv_size, MPI_DOUBLE, dest, 2,
                     MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recv_count);

        for (int j = 0; j < recv_count; j++)
        {
          meshfunction.set(
              mesh.distdata().get_vertex_local(ghost_buff[dest][j]),
              static_cast<T>(recv_buff[j]));
        }

        send_buff.clear();
      }

      delete[] ghost_buff;
      delete[] recv_buff;
      delete[] recv_ghost;

    }

    MPI_File_close(&fh);

    delete[] values;

#else
    error("MPI I/O required for reading mesh function from binary files");
#endif
  }
//----------------------------------------------------------------------------

}
