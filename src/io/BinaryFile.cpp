// Copyright (C) 2009-2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/io/BinaryFile.h>

#include <dolfin/common/byteswap.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/Vector.h>
#include <dolfin/main/MPI.h>
#include <dolfin/main/PE.h>
#include <dolfin/math/LinearDistribution.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>
#include <dolfin/parameter/parameters.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <list>
#include <stdexcept>
#include <string>
#include <typeinfo>

namespace dolfin
{

//-----------------------------------------------------------------------------

typedef struct atomic_cell
{
  size_t const   size;
  size_t * const v;

  //-----------------------------------
  atomic_cell( size_t d )
    : size( d )
    , v( new size_t[size] )
  {
  }
  //-----------------------------------
  atomic_cell( atomic_cell const & other )
    : size( other.size )
    , v( new size_t[size] )
  {
    std::copy( other.v, other.v + size, v );
  }
  //-----------------------------------
  auto operator=( atomic_cell const & other ) -> atomic_cell &
  {
    if ( &other == this )
    {
      return *this;
    }
    if ( size != other.size )
    {
      error( "Size of atomic_cells in assignment do not match" );
    }
    std::copy( other.v, other.v + size, v );
    return *this;
  }
  //-----------------------------------
  ~atomic_cell()
  {
    delete[] v;
  }
} atomic_cell;

//----------------------------------------------------------------------------

BinaryFile::BinaryFile(const std::string filename) :
    GenericFile("Binary", filename),
#if defined( HAVE_MPI )
    t_( nullptr ),
#endif
    version_(BINARY_VERSION)
{
}
//----------------------------------------------------------------------------
BinaryFile::BinaryFile(const std::string filename, real const& t) :
    GenericFile("Binary", filename),
#if defined( HAVE_MPI )
    t_(&t),
#endif
    version_(BINARY_VERSION)
{
#if !defined( HAVE_MPI )
  MAYBE_UNUSED(t);
#endif
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
  uint pe_rank = MPI::rank();
  uint pe_size = MPI::size();

  MPI_File fh;
  MPI_Offset byte_offset;
  BinaryFileHeader hdr;
  bool byteswap;

  MPI::file_open( fh, filename, MPI_MODE_RDONLY );
  byte_offset = MPI::file_read_all( fh, hdr, sizeof(BinaryFileHeader) );
  byteswap = hdr_check(hdr, BINARY_VECTOR_DATA, pe_size);
  byte_offset += MPI::file_read_at_all( fh, &offset[0], 2,
                                        byte_offset + pe_rank * 2 * sizeof(uint),
                                        pe_size * 2 );
  if (byteswap)
  {
    offset[0] = bswap(offset[0]);
    offset[1] = bswap(offset[1]);
  }
  size = offset[1];

#else
  std::ifstream fp(filename.c_str(), std::ifstream::binary);

  if ( not fp.good() )
    throw std::runtime_error( "Failed to open file \"" + filename + "\"" );

  fp.read((char *)&size, sizeof(uint));
#endif

  real *values = new real[size];

#ifdef ENABLE_MPIIO
  MPI::file_read_at_all( fh, values, offset[1],
                         byte_offset + offset[0] * sizeof(real) );
  if (byteswap)
  {
    for(uint i = 0; i < offset[1]; ++i)
    {
      values[i] = bswap(values[i]);
    }
  }
  MPI::file_close( fh );
#else
  fp.read((char *)values, size * sizeof(real));
  fp.close();
#endif

  x.init(size);
  x.set(values);
  delete[] values;

  message(1, "BinaryFile: Read vector from file %s in binary format.",
          filename.c_str());
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
  uint pe_rank = MPI::rank();
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = MPI::size();
  hdr.type = BINARY_VECTOR_DATA;
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  MPI_File fh;
  MPI_Offset byte_offset;

  MPI::file_open( fh, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE );
  byte_offset = MPI::file_write_all( fh, hdr, sizeof( BinaryFileHeader ) );
  byte_offset += MPI::file_write_at_all( fh, &offset[0], 2,
                                         byte_offset + pe_rank * 2 * sizeof(uint),
                                         hdr.pe_size * 2 );
  MPI::file_write_at_all( fh, values, offset[1],
                          byte_offset + offset[0] * sizeof(real) );

  MPI::file_close( fh );
#else
  std::ofstream fp(filename.c_str(), std::ofstream::binary);

  if ( not fp.good() )
    throw std::runtime_error( "Failed to open file \"" + filename + "\"" );

  fp.write((char *)&size, sizeof(uint));
  fp.write((char *)values, x.local_size() * sizeof(real));
  fp.close();
#endif

  delete[] values;

  message(1, "BinaryFile: Saved vector to file %s in binary format.",
          filename.c_str());
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(Function & f)
{

#ifdef ENABLE_MPIIO
  uint pe_size = MPI::size();

  MPI_File fh;
  MPI_Offset byte_offset;
  BinaryFileHeader hdr;
  bool byteswap;

  MPI::file_open( fh, filename, MPI_MODE_RDONLY );
  byte_offset = MPI::file_read_all( fh, hdr, sizeof(BinaryFileHeader) );
  byteswap = hdr_check(hdr, BINARY_FUNCTION_DATA, pe_size);

  uint nfunc;
  byte_offset += MPI::file_read_at_all( fh, nfunc, byte_offset );
  if (byteswap) nfunc = bswap(nfunc);
  if (nfunc > 1)
  {
    warning("File contains %d functions, using first with matching dim.",
            nfunc);
  }

  BinaryFunctionHeader f_hdr;
  for (uint i = 0; i < nfunc; ++i)
  {
    byte_offset += MPI::file_read_at_all( fh, f_hdr, byte_offset );
    if (byteswap) bswap_func_hdr(f_hdr);

    /* Load function if dimension match */
    if (f_hdr.dim == f.value_size())
    {

      uint size = f.value_size() * f.mesh().topology().num_owned(0);
      real *values = new real[size];
      MPI::file_read_at_all( fh, values, size,
                             byte_offset + f.vector().offset() * sizeof(real) );
      if (byteswap)
      {
        for (uint j = 0; j < size; ++j)
        {
          values[j] = bswap(values[j]);
        }
      }
      f.vector().set(values);
      delete[] values;

      MPI::file_close( fh );


      message(1, "BinaryFile: Read function from file %s in binary format.",
             filename.c_str());

      return;
    }

    /* Otherwise continue searching*/
    byte_offset += f_hdr.dim * f.mesh().global_size(0) * sizeof(real);
  }

  MPI::file_close( fh );

  error("No matching functions found in binary file");
#else
  MAYBE_UNUSED(f);
  error("MPI I/O required for loading functions written in  binary");
#endif
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(LabelList<Function>& f)
{
#ifdef ENABLE_MPIIO
  uint pe_size = MPI::size();

  MPI_File fh;
  MPI_Offset byte_offset;
  BinaryFileHeader hdr;
  bool byteswap;

  MPI::file_open( fh, filename, MPI_MODE_RDONLY );
  byte_offset = MPI::file_read_all( fh, hdr, sizeof(BinaryFileHeader) );
  byteswap = hdr_check(hdr, BINARY_FUNCTION_DATA, pe_size);

  uint nfunc;
  byte_offset += MPI::file_read_at_all( fh, nfunc, byte_offset );
  if (byteswap) nfunc = bswap(nfunc);

  if (nfunc != f.size())
  {
    error("Number of functions mismatch between set and file");
  }

  BinaryFunctionHeader f_hdr;
  for ( Label< Function > & fun : f )
  {
    byte_offset += MPI::file_read_at_all( fh, f_hdr, byte_offset );
    if (byteswap) bswap_func_hdr(f_hdr);

    Function * u = fun.first;

    if (f_hdr.dim != u->value_size())
    {
      error("Dimension of file and function set does not match");
    }

    uint size = u->value_size() * u->mesh().topology().num_owned(0);
    Array< real > values( size );

    MPI::file_read_at_all( fh, values.data(), size,
                           byte_offset + u->vector().offset() * sizeof(real) );
    if(byteswap)
    {
      for (uint i = 0; i < size; ++i)
      {
        values[i] = bswap(values[i]);
      }
    }
    u->vector().set(values.data());

    byte_offset += f_hdr.dim * u->mesh().global_size(0) * sizeof(real);
  }


  message(1, "BinaryFile: Read functions from file %s in binary format.",
          filename.c_str());

  MPI::file_close( fh );

#else
  MAYBE_UNUSED(f);
  error("MPI I/O required for loading functions written in  binary");
#endif
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(Function & u)
{
  LabelList<Function> tmp(1, Label<Function>(u, "U"));
  write_function(tmp);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(LabelList<Function>& list)
{
  write_function(list);
}
//----------------------------------------------------------------------------
void BinaryFile::write_function(
    LabelList<Function>& f)
{

#ifdef ENABLE_MPIIO

  nameUpdate(counter);

  BinaryFileHeader hdr;
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = MPI::size();
  hdr.type = BINARY_FUNCTION_DATA;
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  MPI_File fh;
  MPI_Offset byte_offset;
  MPI::file_open( fh, bin_filename_, MPI_MODE_WRONLY | MPI_MODE_CREATE );
  byte_offset = MPI::file_write_all( fh, hdr, sizeof( BinaryFileHeader ) );
  byte_offset += MPI::file_write_all( fh, static_cast<uint>( f.size() ) );

  // Assume same mesh for all data arrays
  Mesh& mesh = f[0].first->mesh();

  for ( Label< Function > & fun : f )
  {
    Function* u = fun.first;

    std::string& name = fun.second;

    // Get number of components
    uint const value_dim = u->value_size();

    // Allocate memory for function values at vertices
    uint size = mesh.size(0) * u->value_size();
    Array< real > values( size );
    uint offset = u->vector().offset();

    if ((u->vector().local_size() / value_dim) != mesh.topology().num_owned(0))
    {
      Array< real > interp_values( size );

      // Get function values at vertices
      u->interpolate_vertex_values(interp_values.data());

      uint ii = 0;
      for (VertexIterator v(mesh); !v.end(); ++v)
      {
        if (!v->is_ghost())
        {
          for (uint i = 0; i < value_dim; ++i)
          {
            values[ii++] = interp_values[v->index() + i * mesh.size(0)];
          }
        }
      }

      // Compute new vertex based offset
      uint num_values = value_dim * mesh.topology().num_owned(0);
      offset = 0;
      MPI::exscan_sum( &num_values, &offset, 1 );
    }
    else
    {
      u->vector().get(values.data());
    }

    size = value_dim * mesh.topology().num_owned(0);

    BinaryFunctionHeader f_hdr;
    f_hdr.dim = value_dim;
    f_hdr.size = value_dim * mesh.global_size(0);
    if (name.length() > FNAME_LENGTH) error("Function name too long.");
    strcpy(&f_hdr.name[0], name.c_str());
    f_hdr.t = t_ ? *t_ : counter;

    byte_offset += MPI::file_write_at_all( fh, f_hdr, byte_offset );
    byte_offset += MPI::file_write_at_all( fh, values.data(), size,
                                           byte_offset + offset * sizeof(real),
                                           value_dim * mesh.global_size(0) );
  }

  MPI::file_close( fh );

  counter++;


  message(1, "BinaryFile: Saved function to file %s in binary format.",
          bin_filename_.c_str());
#else
  MAYBE_UNUSED(f);
  error("MPI I/O is required to save functions in Binary.");
#endif
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(Mesh& mesh)
{

  bool byteswap;
  BinaryFileHeader hdr;

  if (PE::size() == 1 || dolfin_get<bool>("Mesh read in serial"))
  {
    if(PE::rank() > 0)
    {
      error("Reading serial mesh in parallel not implemented");
    }

    std::ifstream fp(filename.c_str(), std::ifstream::binary);

    if ( not fp.good() )
      throw std::runtime_error( "Failed to open file \"" + filename + "\"" );

    uint type = 0;
    uint gdim = 0;
    fp.read((char *) &hdr, sizeof(BinaryFileHeader));
    byteswap = hdr_check(hdr, BINARY_MESH_DATA, PE::size());
    fp.read((char *) &gdim, sizeof(uint));
    fp.read((char *) &type, sizeof(uint));

    if (byteswap){ gdim = bswap(gdim); type = bswap(type); }

    // Create cell type to get topological dimension and number of vertices
    CellType::Type ctype = BinaryFile::cell_type(version_, type);
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

    if (byteswap) { bswap(vertex_data, vertex_data + vertex_data_size); }

    for (uint v = 0; v < num_vertices; ++v)
    {
      editor.add_vertex(v, &vertex_data[v * gdim]);
    }
    delete[] vertex_data;

    // Read cell data
    uint num_cells = 0;
    fp.read((char *) &num_cells, sizeof(uint));
    if(byteswap) num_cells = bswap(num_cells);
    editor.init_cells(num_cells);

    std::vector< uint > cell_data( num_cells * num_cellvertices );
    fp.read((char *) cell_data.data(), cell_data.size() * sizeof(uint));

    if ( byteswap )
    {
      bswap( cell_data.data(), cell_data.data() + cell_data.size() );
    }

    std::vector< size_t > cell_data2( cell_data.size() );
    std::copy( cell_data.begin(), cell_data.end(), cell_data2.begin() );

    for (uint c = 0; c < num_cells; ++c)
    {
      editor.add_cell(c, &cell_data2[c * num_cellvertices]);
    }

    editor.close();
    fp.close();
  }
  else
  {
#ifdef ENABLE_MPIIO

    MPI::Communicator& comm = MPI::DOLFIN_COMM;
    uint pe_size = MPI::size();
    uint pe_rank = MPI::rank();

    MPI_File fh;
    MPI_Offset byte_offset;
    MPI::file_open( fh, filename, MPI_MODE_RDONLY );

    // Read binary header
    uint gdim, type, num_vertices;
    byte_offset = MPI::file_read_all( fh, hdr, sizeof(BinaryFileHeader) );
    byteswap = hdr_check(hdr, BINARY_MESH_DATA, pe_size);
    byte_offset += MPI::file_read_all( fh, gdim );
    byte_offset += MPI::file_read_all( fh, type );
    byte_offset += MPI::file_read_all( fh, num_vertices );

    if (byteswap)
    {
      gdim = bswap(gdim);
      type = bswap(type);
      num_vertices = bswap(num_vertices);
    }

    // Create cell type to get topological dimension and number of vertices
    CellType::Type ctype = BinaryFile::cell_type(version_, type);
    CellType * cell_type = CellType::create(ctype);
    uint const num_cellvertices = cell_type->num_entities(0);
    delete cell_type;

    // Load entities from the file: cells and vertices are distributed linearly
    // so there is not need to compute the offsets.
    LinearDistribution vdist(num_vertices, pe_size, pe_rank);
    uint const vertex_offset = vdist.offset * gdim;
    uint const vertex_data   = vdist.size * gdim;
    real * vertex_buffer = new real[vertex_data];
    byte_offset += MPI::file_read_at_all( fh, vertex_buffer, vertex_data,
                                          byte_offset + vertex_offset * sizeof(real),
                                          gdim * num_vertices );
    if (byteswap) { bswap(vertex_buffer, vertex_buffer + vertex_data); }

    uint num_cells;
    byte_offset += MPI::file_read_at_all( fh, num_cells, byte_offset );
    if (byteswap) num_cells = bswap(num_cells);

    LinearDistribution cdist(num_cells, pe_size, pe_rank);
    uint const cell_offset = cdist.offset * num_cellvertices;
    uint const cell_data   = cdist.size * num_cellvertices;
    uint * cell_buffer = new uint[cell_data];
    MPI::file_read_at_all( fh, cell_buffer, cell_data,
                           byte_offset + cell_offset * sizeof(uint) );

    if (byteswap) { bswap(cell_buffer, cell_buffer + cell_data); }

    MPI::file_close( fh );

    _set<uint> ghosted_vertices;

    // Parse cells: cells are assigned to processes based on the ownership of
    // the first vertex.
    Array<atomic_cell> cells;
    Array< Array<uint> > non_local_cells( pe_size );
    _ordered_set<uint> owned_vertices;
    atomic_cell cell(num_cellvertices);
    for (uint * c = cell_buffer; c !=(cell_buffer + cell_data);
         c += num_cellvertices)
    {
      uint const owner = vdist.owner(*c);
      if (owner == pe_rank)
      {
        std::copy(c, c + num_cellvertices, cell.v);

        owned_vertices.insert(cell.v[0]);
        for (uint n = 1; n < num_cellvertices; ++n)
        {
          if (vdist.owner(cell.v[n]) != pe_rank)
          {
            ghosted_vertices.insert(cell.v[n]);
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
        append( non_local_cells[owner], c, c + num_cellvertices );
      }
    }
    delete[] cell_buffer;

    dolfin_assert(vdist.size >= owned_vertices.size());

    // Send non-local cells to their owner
    {
      uint buff_size = max_array_size( non_local_cells );
      MPI::all_reduce_in_place<MPI::max>( buff_size, comm );
      Array< uint > recv_buffer( buff_size );

      // Exchange data
      for (uint i = 1; i < pe_size; ++i)
      {
        int src = (pe_rank - i + pe_size) % pe_size;
        int dst = (pe_rank + i) % pe_size;

        int num_recv = MPI::sendrecv( non_local_cells[dst], dst,
                                      recv_buffer, src, 1, comm );

        // Add received cells
        for (int j = 0; j < num_recv; j += num_cellvertices)
        {
          std::copy(&recv_buffer[j], &recv_buffer[j + num_cellvertices], cell.v);
          owned_vertices.insert(cell.v[0]);
          cells.push_back(cell);
          for (uint n = 1; n < num_cellvertices; ++n)
          {
            if (vdist.owner(cell.v[n]) != pe_rank)
            {
              ghosted_vertices.insert(cell.v[n]);
            }
            else
            {
              owned_vertices.insert(cell.v[n]);
            }
          }
        }
      }
    }

    // The local number of vertices is now known
    uint const num_local_vertices = owned_vertices.size()
                                    + ghosted_vertices.size();

    // Compute the set of vertices which are in the process range but not
    // referenced in any local cell: all adjacent cells to an orphan vertex are
    // ordered such that their first vertex belongs to a rank different than the
    // current process. This is especially a problem with refined meshes which
    // are badly numbered.
    _ordered_set<uint> orphaned_vertices;
    {
      uint v0 = vdist.offset;
      uint const v1 = v0 + vdist.size;
      dolfin_assert(vdist.size >= owned_vertices.size());
      if (owned_vertices.size() < vdist.size)
      {
        for (_ordered_set<uint>::const_iterator it = owned_vertices.begin();
             it != owned_vertices.end(); ++it, ++v0)
        {
          for (; v0 < (*it); ++v0) { orphaned_vertices.insert(v0); }
        }
        for (; v0 < v1; ++v0) { orphaned_vertices.insert(v0); }
      }
    }
    dolfin_assert(vdist.size == owned_vertices.size()+orphaned_vertices.size());

    // Open mesh for editing
    MeshEditor editor(mesh, ctype, gdim);
    editor.init_vertices(num_local_vertices);
    editor.init_cells(cells.size());

    /* Add local indices
     * Use start_index for global number
     *
     * Exchange ghost vertices
     * Use ghosts to know global number
     *
     * Add cells, remember to map global->local
     */

    uint local_vertex_index = 0;
    for (_ordered_set<uint>::iterator it = owned_vertices.begin();
        it != owned_vertices.end(); ++local_vertex_index, ++it)
    {
      mesh.distdata()[0].set_map(local_vertex_index, *it);
      editor.add_vertex(local_vertex_index,
                        &vertex_buffer[gdim *(*it - vdist.offset)]);
    }

    // Sort ghost vertices by owner and clear the set (not needed anymore)
    Array< Array<uint> > ghosts( pe_size );
    for ( uint const & gvert : ghosted_vertices )
      ghosts[vdist.owner(gvert)].push_back(gvert);

    ghosted_vertices.clear();

    // Exchange ghost vertices coordinates
    {
      uint buff_size = max_array_size( ghosts );
      MPI::all_reduce_in_place<MPI::max>( buff_size );

      Array< uint > recv_buffer( buff_size );

      Array< uint > send_new_owner( buff_size );
      Array< uint > recv_new_owner( buff_size );

      Array< real > send_buffer_coords( buff_size * gdim );
      Array< real > recv_buffer_coords( buff_size * gdim );

      // Exchange data
      _map<uint, uint> new_owner;
      for (uint i = 1; i < pe_size; ++i)
      {
        int src = (pe_rank - i + pe_size) % pe_size;
        int dst = (pe_rank + i) % pe_size;

        int num_recv = MPI::sendrecv( ghosts[dst], dst, recv_buffer, src, 1, comm );

        /*
         * Check if orphaned
         * Send back global number and orphaned info
         */

        send_new_owner.resize( num_recv );
        send_buffer_coords.resize( num_recv * gdim );
        for (int k = 0; k < num_recv; ++k)
        {
          if (orphaned_vertices.find(recv_buffer[k]) != orphaned_vertices.end())
          {
            if (new_owner.find(recv_buffer[k]) == new_owner.end())
            {
              new_owner[recv_buffer[k]] = src;
            }
            send_new_owner[k] = new_owner[recv_buffer[k]];
          }
          else
          {
            send_new_owner[k] = pe_rank;
          }

          uint const v_index = recv_buffer[k] - vdist.offset;
          for (uint l = 0; l < gdim; ++l)
          {
            send_buffer_coords[k*gdim+l] = vertex_buffer[(gdim * v_index) + l];
          }
        }

        MPI::sendrecv( send_new_owner, src, recv_new_owner, dst, 1, comm );

        num_recv = MPI::sendrecv( send_buffer_coords, src,
                                  recv_buffer_coords, dst, 1, comm );

        int g_i = 0;
        for (int k = 0; k < num_recv; local_vertex_index++, ++g_i, k += gdim)
        {
          mesh.distdata()[0].set_map(local_vertex_index, ghosts[dst][g_i]);
          if (recv_new_owner[g_i] != pe_rank)
          {
            mesh.distdata()[0].set_ghost(local_vertex_index, recv_new_owner[g_i]);
          }
          editor.add_vertex(local_vertex_index, &recv_buffer_coords[k]);
        }
      }
    }

    delete[] vertex_buffer;

    // Add cell connectivities
    uint local_cell_index = 0;
    Array< size_t > connectivity( num_cellvertices );
    for (Array<atomic_cell>::iterator it = cells.begin();
         it != cells.end(); ++local_cell_index, ++it)
    {
      mesh.distdata()[0].get_local(it->size ,it->v, connectivity.data());
      editor.add_cell(local_cell_index, connectivity.data());
    }

    editor.close();

#else
    error("MPI and MPI I/O is required for reading binary meshes in parallel");
#endif

  }

  message(1, "BinaryFile: Read Mesh from file %s in binary format.",
          filename.c_str());
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(Mesh& mesh)
{
  uint const gdim = mesh.geometry_dimension();
  uint const type = BinaryFile::cell_type(BINARY_VERSION, mesh.type().cellType());
  uint const num_vertices = mesh.global_size(0);
  uint const num_cells = mesh.num_global_cells();
  uint const num_cellvertices = mesh.type().num_entities(0);

  BinaryFileHeader hdr;
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = PE::size();
  hdr.type = BINARY_MESH_DATA;
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  if (PE::size() == 1 || ! mesh.is_distributed())
  {
    if(PE::rank() > 0)
    {
      error("Writing serial mesh in parallel not implemented");
    }

    std::ofstream fp(filename.c_str(), std::ofstream::binary);

    if ( not fp.good() )
      throw std::runtime_error( "Failed to open file \"" + filename + "\"" );

    // Write Header
    fp.write( reinterpret_cast<const char*>( &hdr ), sizeof(BinaryFileHeader));
    fp.write( reinterpret_cast<const char*>( &gdim ), sizeof(uint));
    fp.write( reinterpret_cast<const char*>( &type ), sizeof(uint));

    // Write vertices
    fp.write( reinterpret_cast<const char*>( &num_vertices ), sizeof(uint));
    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      fp.write( reinterpret_cast<char*>( v->x() ), gdim * sizeof(real));
    }

    // Write cells
    fp.write( reinterpret_cast<const char*>( &num_cells ), sizeof(uint));
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      fp.write( reinterpret_cast<const char*>( c->entities(0).data() ),
                num_cellvertices * sizeof(uint) );
    }

    fp.close();
  }
  else
  {

#ifdef ENABLE_MPIIO

    MPI::Communicator& comm = MPI::DOLFIN_COMM;

    MPI_File fh;
    MPI_Offset byte_offset;

    MPI::file_open( fh, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE );

    // Write Header
    byte_offset = MPI::file_write_all( fh, hdr, sizeof( BinaryFileHeader ) );
    byte_offset += MPI::file_write_all( fh, gdim );
    byte_offset += MPI::file_write_all( fh, type );

    // Write vertices
    size_t vertex_offset = 0;
    size_t vertex_buffer_size = gdim * mesh.topology().num_owned(0);
    MPI::offset(vertex_buffer_size, vertex_offset, comm);

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
    byte_offset += MPI::file_write_all( fh, num_vertices );
    byte_offset += MPI::file_write_at_all( fh, vertex_buffer, vertex_buffer_size,
                                           byte_offset + vertex_offset * sizeof(real),
                                           gdim * num_vertices );
    delete[] vertex_buffer;

    // Write Cells
    size_t cell_offset = 0;
    size_t cell_buffer_size = num_cellvertices * mesh.num_cells();
    MPI::offset(cell_buffer_size, cell_offset, comm);

    uint * cell_buffer = new uint[cell_buffer_size];
    uint * cp = &cell_buffer[0];
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      for (uint i = 0; i < c->num_entities(0); ++i)
      {
        *(cp++) = mesh.distdata()[0].get_global(c->entities(0)[i]);
      }
    }
    byte_offset += MPI::file_write_at_all( fh, num_cells, byte_offset );
    MPI::file_write_at_all( fh, cell_buffer, cell_buffer_size,
                            byte_offset + cell_offset * sizeof(uint) );

    delete[] cell_buffer;

    MPI::file_close( fh );

#else
    error("MPI and MPI I/O is required for writing a mesh in parallel");
#endif

  }

  message(1, "BinaryFile: Saved mesh to file %s in binary format.",
          filename.c_str());
}
//----------------------------------------------------------------------------
void BinaryFile::nameUpdate(const int counter)
{
  std::string file_start( filename, 0, filename.rfind(".") );

  std::ostringstream file_id;
  file_id.fill('0');
  file_id.width(6);
  file_id << counter;

  bin_filename_ = file_start + file_id.str() + ".bin";
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
void BinaryFile::operator<<(MeshFunction<uint>& meshfunction)
{
  write_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(MeshFunction<real>& meshfunction)
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
  hdr.pe_size = MPI::size();
  hdr.type = BINARY_MESH_FUNCTION_DATA;
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  MPI_File fh;
  MPI_Offset byte_offset;

  MPI::file_open( fh, bin_filename_, MPI_MODE_WRONLY | MPI_MODE_CREATE );
  byte_offset = MPI::file_write_all( fh, hdr, sizeof( BinaryFileHeader ) );

  uint local_size = 0;
  uint mfunc_type = 0;
  if (meshfunction.dim() == mesh.topology_dimension())
  {
    byte_offset += MPI::file_write_at_all( fh, mfunc_type, byte_offset );

    for (CellIterator c(mesh); !c.end(); ++c)
    {
      *(vp++) = (real) meshfunction(c->index());
    }

    local_size = mesh.num_cells();
  }
  else if (meshfunction.dim() == 0)
  {
    mfunc_type = 1;
    byte_offset += MPI::file_write_at_all( fh, mfunc_type, byte_offset );

    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      if (!v->is_ghost())
      {
        *(vp++) = (real) meshfunction(v->index());
      }
    }

    local_size = mesh.topology().num_owned(0);
  }
  else
  {
    error("Binary output of mesh functions is implemented "
          "for cell/vertex-based functions only.");
  }

  uint offset = 0;
  MPI::exscan_sum( &local_size, &offset, 1 );
  MPI::file_write_at_all( fh, values, local_size,
                          byte_offset + offset * sizeof(real) );
  MPI::file_close( fh );

  counter++;

  delete[] values;

  message(1, "BinaryFile: Saved MeshFunction to file %s in binary format.",
          bin_filename_.c_str());

#else
  MAYBE_UNUSED(meshfunction);
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
void BinaryFile::operator>>(MeshFunction<uint>& meshfunction)
{
  read_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(MeshFunction<real>& meshfunction)
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

  uint pe_rank = MPI::rank();
  uint pe_size = MPI::size();
  MPI_File fh;
  MPI_Offset byte_offset;
  BinaryFileHeader hdr;
  bool byteswap;
  MPI::file_open( fh, filename, MPI_MODE_RDONLY );
  byte_offset = MPI::file_read_all( fh, hdr, sizeof(BinaryFileHeader) );
  byteswap = hdr_check(hdr, BINARY_MESH_FUNCTION_DATA, pe_size);

  int mfunc_type = 0;
  byte_offset += MPI::file_read_at_all( fh, mfunc_type, byte_offset );
  if (byteswap) mfunc_type = bswap(mfunc_type);


  if ((mfunc_type == 0 && meshfunction.dim() != mesh.topology_dimension())
      || (mfunc_type == 1 && meshfunction.dim() != 0))
  {
    error("Meshfunction does not match data in file");
  }

  uint local_size = (mfunc_type > 0 ? mesh.topology().num_owned(0)
                                    : mesh.num_cells());

  uint offset = 0;
  MPI::exscan_sum( &local_size, &offset, 1 );
  MPI::file_read_at_all( fh, values, local_size,
                         byte_offset + offset * sizeof(real) );
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
      meshfunction(i) = static_cast<T>(values[i]);
    }
  }
  if (mfunc_type == 1)
  {
    for (uint i = 0; i < meshfunction.size(); ++i)
    {
      meshfunction(i) = static_cast<T>(values[i]);
    }

    std::vector<uint> *ghost_buff = new std::vector<uint>[pe_size];
    for (GhostIterator it(mesh.distdata()[0]); it.valid(); ++it)
    {
      ghost_buff[it.owner()].push_back(it.global_index());
    }

    std::vector<real> send_buff;
    uint recv_size = mesh.distdata()[0].num_ghost();
    int recv_size_gh, send_size;

    for (uint i = 0; i < pe_size; ++i)
    {
      send_size = ghost_buff[i].size();
      MPI::check_error( MPI_Reduce(&send_size, &recv_size_gh, 1, MPI_UNSIGNED,
                                   MPI_SUM, i, MPI::DOLFIN_COMM) );
    }
    uint *recv_ghost = new uint[recv_size_gh];
    real *recv_buff = new real[recv_size];

    for (uint j = 1; j < pe_size; j++)
    {
      int src = (pe_rank - j + pe_size) % pe_size;
      int dest = (pe_rank + j) % pe_size;

      int recv_count = MPI::sendrecv( &ghost_buff[dest][0],
                                      ghost_buff[dest].size(), dest,
                                      recv_ghost, recv_size_gh, src, 1 );

      for (int k = 0; k < recv_count; ++k)
        send_buff.push_back(meshfunction(mesh.distdata()[0].get_local(recv_ghost[k])));

      recv_count = MPI::sendrecv( &send_buff[0], send_buff.size(), src,
                                  recv_buff, recv_size, dest, 2 );

      for (int j = 0; j < recv_count; j++)
      {
        meshfunction(mesh.distdata()[0].get_local(ghost_buff[dest][j]))
            = static_cast<T>(recv_buff[j]);
      }

      send_buff.clear();
    }

    delete[] ghost_buff;
    delete[] recv_buff;
    delete[] recv_ghost;

  }

  MPI::file_close( fh );

  delete[] values;

  message(1, "BinaryFile: Read MesHFunction from file %s in binary format.",
          filename.c_str());

#else
  MAYBE_UNUSED(meshfunction);
  error("MPI I/O required for reading mesh function from binary files");
#endif
}

//-----------------------------------------------------------------------------
auto BinaryFile::hdr_check(BinaryFileHeader& hdr, Binary_data_t type,
                           uint pe_size) -> bool
{

  bool byteswap = false;

  if (hdr.magic == BINARY_MAGIC_V2)
  {
    message(1, "BinaryFile: Loading Binary File format version 2");
    version_ = 2;
  }
  else if (hdr.magic == BINARY_MAGIC_V1)
  {
    message(1, "BinaryFile: Loading Binary File format version 1");
    version_ = 1;
  }
  else if (bswap(hdr.magic) == BINARY_MAGIC_V2)
  {
    message(1, "BinaryFile: Loading Binary File format version 2 (endian conversion)");
    version_ = 2;
    byteswap = true;
  }
  else if (bswap(hdr.magic) == BINARY_MAGIC_V1)
  {
    message(1, "BinaryFile: Loading Binary File format version 1 (endian conversion)");
    version_ = 1;
    byteswap = true;
  }
  else
  {
    error("Corrupt header: invalid magic number (%0x)", hdr.magic);
  }

  if (byteswap)
  {
    hdr.magic = bswap(hdr.magic);
    hdr.bendian = bswap(hdr.bendian);
    hdr.pe_size = bswap(hdr.pe_size);
    hdr.type = static_cast<Binary_data_t>(bswap(hdr.type));
  }

  if (hdr.type != type)
  {
    error("Invalid data type in file");
  }

  if ((hdr.type == BINARY_FUNCTION_DATA || hdr.type == BINARY_VECTOR_DATA)
      && (hdr.pe_size != pe_size))
  {
    error("File stored on %d PEs, currently running on %d PEs", hdr.pe_size,
          pe_size);
  }

  return byteswap;
}

//-----------------------------------------------------------------------------
#ifdef ENABLE_MPIIO
void BinaryFile::bswap_func_hdr(BinaryFunctionHeader& hdr)
{
  hdr.dim = bswap(hdr.dim);
  hdr.size = bswap(hdr.size);
  hdr.t = bswap(hdr.t);
}
#endif

//-----------------------------------------------------------------------------
auto BinaryFile::cell_type(uint version, CellType::Type const type) -> uint
{
  switch (version)
    {
    case 2:
      switch (type)
        {
        case CellType::point:
          return 0;
          break;
        case CellType::interval:
          return 1;
          break;
        case CellType::triangle:
          return 2;
          break;
        case CellType::tetrahedron:
          return 3;
          break;
        case CellType::quadrilateral:
          return 4;
          break;
        case CellType::hexahedron:
          return 6;
          break;
        default:
          error("Unsupported mesh cell type in BinaryFile V2.");
          break;
        }
      break;
    case 1:
      switch (type)
        {
        case CellType::triangle:
          return 0;
          break;
        case CellType::tetrahedron:
          return 1;
          break;
        default:
          error("Unsupported mesh cell type in BinaryFile V1.");
          break;
        }
      break;
    default:
      error("Invalid version compatibility number for cell type detection.");
      break;
    }
  return 0;
}

//-----------------------------------------------------------------------------
auto BinaryFile::cell_type(uint version, uint const type) -> CellType::Type
{
  switch (version)
    {
    case 2:
      switch (type)
        {
        case 0:
          return CellType::point;
          break;
        case 1:
          return CellType::interval;
          break;
        case 2:
          return CellType::triangle;
          break;
        case 3:
          return CellType::tetrahedron;
          break;
        case 4:
          return CellType::quadrilateral;
          break;
        case 6:
          return CellType::hexahedron;
          break;
        default:
          error("Unsupported binary cell type in BinaryFile V2.");
          break;
        }
      break;
    case 1:
      switch (type)
        {
        case 0:
          return CellType::triangle;
          break;
        case 1:
          return CellType::tetrahedron;
          break;
        default:
          error("Unsupported binary cell type in BinaryFile V1.");
          break;
        }
      break;
    default:
      error("Invalid version compatibility number for cell type detection.");
      break;
    }
  return CellType::point;
}

//----------------------------------------------------------------------------

} /* namespace dolfin */
