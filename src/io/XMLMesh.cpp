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
    pre_partitioning_(true),
    cell_type_(NULL),
    vertex_dist_(NULL),
    cell_dist_(NULL),
    cell_count_(0),
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
void XMLMesh::open(std::string const& filename)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool XMLMesh::close()
{
  return state_ == ROOT;
}
//-----------------------------------------------------------------------------
void XMLMesh::clear()
{
  state_ = ROOT;
  delete editor_;
  editor_ = NULL;
  parallel_ = false;
  pre_partitioning_ = false;
  delete cell_type_;
  cell_type_ = NULL;
  delete vertex_dist_;
  vertex_dist_ = NULL;
  delete cell_dist_;
  cell_dist_ = NULL;
  cell_count_ = 0;
  delete [] vertex_owner_;
  vertex_owner_ = NULL;
  nonlocal_vertices_.clear();
  cell_buffer_.clear();
}
//-----------------------------------------------------------------------------
void XMLMesh::beginMesh(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  std::string celltype = parse<std::string>(name, attrs, "celltype");
  if (cell_type_ != NULL)
  {
    error("XMLMesh : cell type is already created.");
  }
  cell_type_ = CellType::create(celltype);
  uint dim = parse<uint>(name, attrs, "dim");
  if (editor_ != NULL)
  {
    error("XMLMesh : mesh editor is already created.");
  }
  dolfin_assert(cell_type_ != NULL);
  editor_ = new MeshEditor(mesh_, *cell_type_, dim);
  parallel_ = mesh_.topology().is_distributed();
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
  uint size = parse<uint>(name, attrs, "size");
  // Calculate a linear data distribution
  if (vertex_dist_ != NULL)
  {
    error("XMLMesh : vertex distribution is already created.");
  }
  vertex_dist_ = new LinearDistribution(size, MPI::numProcesses(), MPI::processNumber());
  // Set number of vertices
  editor_->init_vertices(vertex_dist_->size, size);
  if(parallel_)
  {
    // Set process range for vertices as they are just collected contiguously
    // prior to distribution.
    mesh_.topology().distdata()[0].set_range(vertex_dist_->size, size);
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::readCells(const xmlChar *name, const xmlChar **attrs)
{
  uint size = parse<uint>(name, attrs, "size");
  //
  if (cell_dist_ != NULL)
  {
    error("XMLMesh : cell distribution is already created.");
  }
  cell_count_ = 0;
  if(parallel_)
  {
    editor_->init_cells(0);
    editor_->close();
    dolfin_assert(mesh_.topology().size(0) == vertex_dist_->size);
    dolfin_assert(mesh_.distdata()[0].local_size() == vertex_dist_->size);
    // The following section requires process range and global size to be set
    if (pre_partitioning_)
    {
      MeshFunction<uint> pre_partition;
      mesh_.partition_geom(pre_partition);
      mesh_.distribute(pre_partition);
    }
    if (vertex_owner_ != NULL)
    {
      error("XMLMesh : vertex ownership array is already created.");
    }
    vertex_owner_ = new uint[mesh_.size(0)];
    std::fill_n(&vertex_owner_[0], mesh_.size(0), MPI::numProcesses());
  }
  else
  {
    editor_->init_cells(size);
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::readVertex(const xmlChar *name, const xmlChar **attrs)
{
  // Read index
  uint index = parse<uint>(name, attrs, "index");

  if (!vertex_dist_->in_range(index)) return;

  // Handle differently depending on geometric dimension
  real x[Point::MAX_SIZE];
  switch (mesh_.geometry().dim())
  {
  case 3:
    x[2] = parse<real>(name, attrs, "z");
  case 2:
    x[1] = parse<real>(name, attrs, "y");
  case 1:
    x[0] = parse<real>(name, attrs, "x");
    break;
  default:
    error("Dimension of mesh must be 1, 2 or 3: provided %d.",
          mesh_.geometry().dim());
    break;
  }
  editor_->add_vertex(index - vertex_dist_->offset, &x[0]);
}
//-----------------------------------------------------------------------------
void XMLMesh::readCell(const xmlChar *name, const xmlChar **attrs)
{
  static const char * const vertex_attr[8] =
      { "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7" };
  static uint v[8] = { 0 };
  dolfin_assert(cell_type_ != NULL);
  /*
  ///FIXME: Assume only one cell type for now
  if (xmlStrcmp(name, (xmlChar *) cell_type_->str().c_str()))
  {
    error("XMLMesh : invalid cell element <%s>, expected <%s>",
          (const char *) name, cell_type_->str().c_str());
    if (supported_celltypes_.count(current_cell_.first) == 0)
    {
      error("XMLMesh : invalid cell element <%s>", (const char *) name);
    }
    celltype_.second = checkCellElement(name);
    if (current_cell_.second < 0)
    {
      error("XMLMesh : unknown cell element <%s>", (const char *) name);
    }
  }
  */
  //
  uint index = parse<uint>(name, attrs, "index");

  if (parallel_)
  {
    v[0] = parse<uint>(name, attrs, vertex_attr[0]);
    if (!mesh_.distdata()[0].has_global(v[0]))
    {
     return;
    }
    uint const rank = MPI::processNumber();
    vertex_owner_[mesh_.distdata()[0].get_local(v[0])] = rank;
    cell_buffer_.push_back(v[0]);
    for (uint i = 1; i < cell_type_->num_entities(0); ++i)
    {
      v[i] = parse<uint>(name, attrs, vertex_attr[i]);
      if (mesh_.distdata()[0].has_global(v[i]))
      {
        vertex_owner_[mesh_.distdata()[0].get_local(v[i])] = rank;
      }
      else
      {
        nonlocal_vertices_.insert(v[i]);
      }
      cell_buffer_.push_back(v[i]);
    }
    ++cell_count_;
  }
  else
  {
    for (uint i = 0; i < cell_type_->num_entities(0); ++i)
    {
      v[i] = parse<uint>(name, attrs, vertex_attr[i]);
    }
    editor_->add_cell(index, &v[0]);
    ++cell_count_;
  }
}
//-----------------------------------------------------------------------------
void XMLMesh::endMesh()
{
  if (parallel_)
  {
    uint const rank = MPI::processNumber();
    uint const pe_size = MPI::numProcesses();
    uint const gdim = mesh_.geometry().dim();
    Array<uint> sendbuf(nonlocal_vertices_.size());
    sendbuf.assign(nonlocal_vertices_.begin(), nonlocal_vertices_.end());
    uint const shared = nonlocal_vertices_.size();
    uint orphan = 0;

    // Exchange ghost points
    MPI_Status status;
    uint src;
    uint dst;
    uint recvmax;
    MPI::numGlobalSum(shared, recvmax);
    uint *recvbuf = new uint[recvmax];
    int recvcount;
    uint * sendbuf_v = new uint[2*recvmax];
    real * sendbuf_x = new real[gdim*recvmax];
    uint recvcnt_v = 2*shared;
    uint * recvbuf_v = new uint[recvcnt_v];
    uint * recvptr_v = &recvbuf_v[0];
    uint recvcnt_x = gdim*shared;
    real * recvbuf_x = new real[recvcnt_x];
    real * recvptr_x = &recvbuf_x[0];
    DistributedData& vdata0 = mesh_.distdata()[0];
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dst = (rank + j) % pe_size;

      MPI_Sendrecv(&sendbuf[0], sendbuf.size(), MPI_UNSIGNED, dst, 1, recvbuf,
                   recvmax, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);

      uint count = 0;
      for (int k = 0; k < recvcount; ++k)
      {
        if (vdata0.has_global(recvbuf[k]))
        {
          uint const index = vdata0.get_local(recvbuf[k]);
          sendbuf_v[2*count] = recvbuf[k];
          if (vertex_owner_[index] == pe_size)
          {
            ++orphan;
            vertex_owner_[index] = src;
            sendbuf_v[2*count+1] = src;
          }
          else
          {
            sendbuf_v[2*count+1] = vertex_owner_[index];
          }
          for (uint d = 0; d < gdim; ++d)
          {
            sendbuf_x[gdim*count+d] = mesh_.geometry().x(index)[d];
          }
          ++count;
        }
      }

      MPI_Sendrecv(&sendbuf_v[0], 2 * count, MPI_UNSIGNED, src, 1, recvptr_v,
                   recvcnt_v, MPI_UNSIGNED, dst, 1, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);
      recvcnt_v -= recvcount;
      recvptr_v += recvcount;

      MPI_Sendrecv(&sendbuf_x[0], gdim * count, MPI_DOUBLE, src, 2, recvptr_x,
                   recvcnt_x, MPI_DOUBLE, dst, 2, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
      recvcnt_x -= recvcount;
      recvptr_x += recvcount;
    }
    delete[] sendbuf_x;
    delete[] sendbuf_v;

    // Create mesh editor
    Mesh new_mesh;
    delete editor_;
    editor_ = new MeshEditor(new_mesh, mesh_.type(), gdim);

    // Add vertices
    editor_->init_vertices(mesh_.size(0) + shared - orphan);
    DistributedData& new_distdata0 = new_mesh.distdata()[0];

    uint vertex_count = 0;
    for (VertexIterator vertex(mesh_); !vertex.end(); ++vertex)
    {
      if (vertex_owner_[vertex->index()] == rank)
      {
        new_distdata0.set_map(vertex_count, vertex->global_index());
        editor_->add_vertex(vertex_count, vertex->x());
        ++vertex_count;
      }
      else if (vertex_owner_[vertex->index()] == pe_size)
      {
        error("XMLMesh : vertex %u is unassigned");
      }
    }

    // Add shared ghost vertices
    uint ii = 0;
    uint ci = 0;
    for (uint i = 0; i < shared; ++i, ii+=2, ci += gdim)
    {
      new_distdata0.set_map(vertex_count, recvbuf_v[ii]);
      if (recvbuf_v[ii + 1] != rank)
      {
        new_distdata0.set_ghost(vertex_count, recvbuf_v[ii + 1]);
        dolfin_assert(new_distdata0.is_ghost(vertex_count));
        dolfin_assert(new_distdata0.get_owner(vertex_count) == recvbuf_v[ii + 1]);
      }
      editor_->add_vertex(vertex_count, &recvbuf_x[ci]);
      ++vertex_count;
    }

    // Add cells
    // This code is only valid for homogeneous meshes
    editor_->init_cells(cell_count_);
    uint c = 0;
    uint num_cell_vertices = mesh_.type().num_entities(0);
    uint * connectivity = new uint[num_cell_vertices];
    for (uint i = 0; i < cell_buffer_.size(); i += num_cell_vertices)
    {
      for (uint n = 0; n < num_cell_vertices; ++n)
      {
        connectivity[n] = new_distdata0.get_local(cell_buffer_[i + n]);
      }
      editor_->add_cell(c++, &connectivity[0]);
    }
    delete[] connectivity;
    editor_->close();
    mesh_ = new_mesh;

    sendbuf.clear();
    delete[] recvbuf_v;
    delete[] recvbuf_x;
    delete[] recvbuf;
  }
  else
  {
    editor_->close();
  }
  clear();
  message("XMLMesh : loading time %g s.", toc());
}
//-----------------------------------------------------------------------------

}

#endif /* HAVE_XML */
