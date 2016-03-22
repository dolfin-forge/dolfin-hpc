// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2003-10-21
// Last changed: 2008-07-01

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_XML
#include <cstring>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/timing.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/io/XMLMesh.h>
#include <dolfin/parameter/parameters.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
XMLMesh::XMLMesh(Mesh& mesh) :
    XMLObject(),
    mesh_(mesh),
    state_(ROOT),
    editor_(NULL),
    parallel_(false),
    cell_count_(0),
    vertex_offset_(0),
    vertex_range_end_(0),
    cell_offset_(0),
    cell_range_end_(0),
    local_vertices_(NULL),
    shared_vertices_(NULL),
    num_local_vertices_(0),
    num_local_cells_(0),
    vertex_owner_(NULL)
{
}
//-----------------------------------------------------------------------------
XMLMesh::~XMLMesh()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void XMLMesh::startElement(const xmlChar *name, const xmlChar **attrs)
{
  switch (state_)
  {
  case ROOT:
    if (xmlStrcasecmp(name, (xmlChar *) "mesh") == 0)
    {
      beginMesh(name, attrs);
      state_ = IN_MESH;
    }
    break;
  case IN_MESH:
    if (xmlStrcasecmp(name, (xmlChar *) "vertices") == 0)
    {
      readVertices(name, attrs);
      state_ = IN_VERTICES;
    }
    else if (xmlStrcasecmp(name, (xmlChar *) "cells") == 0)
    {
      readCells(name, attrs);
      state_ = IN_CELLS;
    }
    break;
  case IN_VERTICES:
    if (xmlStrcasecmp(name, (xmlChar *) "vertex") == 0)
    {
      readVertex(name, attrs);
    }
    break;
  case IN_CELLS:
    readCell(name, attrs);
    break;
  default:
    break;
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::endElement(const xmlChar *name)
{
  switch (state_)
  {
  case IN_MESH:
    if (xmlStrcasecmp(name, (xmlChar *) "mesh") == 0)
    {
      endMesh();
      state_ = ROOT;
    }
    break;
  case IN_VERTICES:
    if (xmlStrcasecmp(name, (xmlChar *) "vertices") == 0)
    {
      state_ = IN_MESH;
    }
    break;
  case IN_CELLS:
    if (xmlStrcasecmp(name, (xmlChar *) "cells") == 0)
    {
      state_ = IN_MESH;
    }
    break;
  default:
    break;
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::open(std::string filename)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool XMLMesh::close()
{
  return state_ == ROOT;
}
//-----------------------------------------------------------------------------
void XMLMesh::beginMesh(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  std::string type = parseString(name, attrs, "celltype");
  uint gdim = parseUnsignedInt(name, attrs, "dim");
  if (editor_ != NULL)
  {
    error("XMLMesh : mesh editor is already created.");
  }
  CellType * cell_type = CellType::create(type);
  editor_ = new MeshEditor(mesh_, cell_type->cellType(), gdim);
  delete cell_type;
  parallel_ = (MPI::numProcesses() > 1 && !dolfin_get("Mesh read in serial"));
  if(parallel_)
  {
    warning("Reading DOLFIN xml meshes in parallel is deprecated.\n"
            "For better I/O performance, consider converting to flat binary");
  }
  tic();
}
//-----------------------------------------------------------------------------
void XMLMesh::readVertices(const xmlChar *name, const xmlChar **attrs)
{
  uint num_vertices = parseUnsignedInt(name, attrs, "size");
  //
  cell_count_ = 0;
  uint pe_size = MPI::numProcesses();
  uint rank = MPI::processNumber();
  // Calculate a linear data distribution
  uint L = std::floor((real) num_vertices / (real) pe_size);
  uint R = num_vertices % pe_size;
  num_local_vertices_ = (num_vertices + pe_size - rank - 1) / pe_size;
  vertex_offset_ = rank * L + std::min(rank, R);
  vertex_range_end_ = vertex_offset_ + (num_local_vertices_ - 1);
  // Set number of vertices
  editor_->initVertices(num_local_vertices_);

  if (parallel_)
  {
    mesh_.distdata().set_num_global(0, num_vertices);
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::readCells(const xmlChar *name, const xmlChar **attrs)
{
  uint num_cells = parseUnsignedInt(name, attrs, "size");
  //
  cell_count_ = 0;
  if(parallel_)
  {
    editor_->initCells(1);
    editor_->close();
//    MeshFunction<uint> pre_partition;
//    mesh_.partition_geom(pre_partition);
//    mesh_.distribute(pre_partition);
    vertex_owner_ = new uint[mesh_.size(0)];
    std::fill_n(&vertex_owner_[0], mesh_.size(0), MPI::numProcesses());
  }
  else
  {
    editor_->initCells(num_cells);
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::readVertex(const xmlChar *name, const xmlChar **attrs)
{
  // Read index
  uint v = parseUnsignedInt(name, attrs, "index");

  if (v < vertex_offset_ || v > vertex_range_end_) return;

  // Handle differently depending on geometric dimension
  real x[Point::MAX_SIZE];
  switch (mesh_.geometry().dim())
  {
  case 3:
    x[2] = parseReal(name, attrs, "z");
  case 2:
    x[1] = parseReal(name, attrs, "y");
  case 1:
    x[0] = parseReal(name, attrs, "x");
    break;
  default:
    error("Dimension of mesh must be 1, 2 or 3: provided %d.",
          mesh_.geometry().dim());
    break;
  }
  editor_->addVertex(v - vertex_offset_, &x[0]);
  if (parallel_)
  {
    mesh_.distdata().set_map(v - vertex_offset_, v, 0);
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::readCell(const xmlChar *name, const xmlChar **attrs)
{
  static const char * const vertex_attr[8] =
      { "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7" };
  static uint v[8] = { 0 };
  //
  uint c = parseUnsignedInt(name, attrs, "index");

  if (parallel_)
  {
    v[0] = parseUnsignedInt(name, attrs, vertex_attr[0]);
    if (!mesh_.distdata().has_global(v[0], 0))
    {
     return;
    }
    uint const rank = MPI::processNumber();
    vertex_owner_[mesh_.distdata().get_vertex_local(v[0])] = rank;
    cell_buffer_.push_back(v[0]);
    for (uint i = 1; i < mesh_.type().numEntities(0); ++i)
    {
      v[i] = parseUnsignedInt(name, attrs, vertex_attr[i]);
      if (mesh_.distdata().has_global(v[i], 0))
      {
        vertex_owner_[mesh_.distdata().get_vertex_local(v[i])] = rank;
      }
      else
      {
        nonlocal_vertices_.insert(v[i]);
      }
      cell_buffer_.push_back(v[i]);
    }
  }
  else
  {
    for (uint i = 0; i < mesh_.type().numEntities(0); ++i)
    {
      v[i] = parseUnsignedInt(name, attrs, vertex_attr[i]);
    }
    editor_->addCell(c, &v[0]);
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::endMesh()
{
  if (parallel_)
  {
    Mesh new_mesh;
    delete editor_;
    uint const gdim = mesh_.geometry().dim();
    editor_ = new MeshEditor(new_mesh, mesh_.type().cellType(), gdim);

    uint const rank = MPI::processNumber();
    uint const pe_size = MPI::numProcesses();

    Array<uint> sendbuf(nonlocal_vertices_.size());
    sendbuf.assign(nonlocal_vertices_.begin(), nonlocal_vertices_.end());
    uint const shared = nonlocal_vertices_.size();
    uint orphan = 0;

    //
    MPI_Status status;
    uint src, dst;
    int recvmax;
    MPI_Allreduce(&shared, &recvmax, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
    uint *recvbuf = new uint[recvmax];
    int recvcount;
    uint * sendbuf_idx = new uint[2*recvmax];
    real * sendbuf_crd = new real[gdim*recvmax];
    uint idxbuf_size = 2*shared;
    uint *idxbuf = new uint[idxbuf_size];
    uint *idxptr = &idxbuf[0];
    uint crdbuf_size = gdim*shared;
    real *crdbuf = new real[crdbuf_size];
    real *crdptr = &crdbuf[0];
    // Exchange ghost points
    for (uint j = 1; j < pe_size; ++j)
    {

      src = (rank - j + pe_size) % pe_size;
      dst = (rank + j) % pe_size;

      MPI_Sendrecv(&sendbuf[0], sendbuf.size(), MPI_UNSIGNED, dst, 1,
                   recvbuf, recvmax, MPI_UNSIGNED, src, 1,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);

      uint count = 0;
      for (int k = 0; k < recvcount; ++k)
      {
        if (mesh_.distdata().has_global(recvbuf[k], 0))
        {
          uint const index = mesh_.distdata().get_vertex_local(recvbuf[k]);
          sendbuf_idx[2*count] = recvbuf[k];
          if (vertex_owner_[index] == pe_size)
          {
            ++orphan;
            vertex_owner_[index] = src;
            sendbuf_idx[2*count+1] = src;
          }
          else
          {
            sendbuf_idx[2*count+1] = vertex_owner_[index];
          }
          for (uint d = 0; d < gdim; ++d)
          {
            sendbuf_crd[gdim*count+d] = mesh_.geometry().x(index)[d];
          }
          ++count;
        }
      }

      MPI_Sendrecv(&sendbuf_idx[0], 2*count, MPI_UNSIGNED, src, 1, idxptr,
                   idxbuf_size, MPI_UNSIGNED, dst, 1, MPI::DOLFIN_COMM,
                   &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);
      idxbuf_size -= recvcount;
      idxptr += recvcount;

      MPI_Sendrecv(&sendbuf_crd[0], gdim*count, MPI_DOUBLE, src, 2,
                   crdptr, crdbuf_size, MPI_DOUBLE, dst, 2, MPI::DOLFIN_COMM,
                   &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
      crdbuf_size -= recvcount;
      crdptr += recvcount;
    }
    delete[] sendbuf_crd;
    delete[] sendbuf_idx;

    // Init new mesh
    editor_->initVertices(mesh_.numVertices() + shared - orphan);
    new_mesh.distdata().set_num_global(0, mesh_.global_numVertices());

    uint vertex_count = 0;
    for (VertexIterator vertex(mesh_); !vertex.end(); ++vertex)
    {
      if (vertex_owner_[vertex->index()] == rank)
      {
        new_mesh.distdata().set_map(vertex_count,
                                    mesh_.distdata().get_global(*vertex), 0);
        editor_->addVertex(vertex_count, vertex->point());
        ++vertex_count;
      }
    }

    //Add shared ghost vertices
    uint ii = 0;
    uint ci = 0;
    for (uint i = 0; i < shared; ++i, ii+=2, ci += gdim)
    {
      new_mesh.distdata().set_map(vertex_count, idxbuf[ii], 0);
      if (idxbuf[ii + 1] != rank)
      {
        new_mesh.distdata().set_ghost(vertex_count, 0);
        new_mesh.distdata().set_ghost_owner(vertex_count, idxbuf[ii + 1], 0);
      }
      editor_->addVertex(vertex_count, &crdbuf[ci]);
      ++vertex_count;
    }

    uint ndims = mesh_.type().numVertices(mesh_.topology().dim());
    editor_->initCells(cell_buffer_.size() / ndims);
    uint c = 0;
    uint * connectivity = new uint[ndims];
    for (uint i = 0; i < cell_buffer_.size(); i += ndims)
    {
      for (uint n = 0; n < ndims; ++n)
      {
        connectivity[n] =
            new_mesh.distdata().get_vertex_local(cell_buffer_[i + n]);
      }
      editor_->addCell(c++, &connectivity[0]);
    }
    delete[] connectivity;
    editor_->close();
    mesh_ = new_mesh;

    nonlocal_vertices_.clear();
    sendbuf.clear();
    cell_buffer_.clear();
    delete[] idxbuf;
    delete[] crdbuf;
    delete[] recvbuf;
    delete[] vertex_owner_;
    delete editor_;
    editor_ = NULL;
  }
  else
  {
    editor_->close();
    delete editor_;
    editor_ = NULL;
  }
  message("XMLMesh : loading took %g s.", toc());
}
//-----------------------------------------------------------------------------

}

#endif /* HAVE_XML */
