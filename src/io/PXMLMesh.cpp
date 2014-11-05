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
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/io/PXMLMesh.h>
#include <dolfin/parameter/parameters.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
PXMLMesh::PXMLMesh(Mesh& mesh) :
    XMLObject(),
    mesh_(mesh),
    state_(OUTSIDE),
    editor_(NULL),
    f_(NULL),
    a_(NULL),
    numParsedVertices_(0),
    numParsedCells_(0),
    startIndex_vert_(0),
    endIndex_vert_(0),
    startIndex_cell_(0),
    endIndex_cell_(0),
    local_vertices_(NULL),
    shared_vertices_(NULL),
    numLocalVertices_(0),
    numLocalCells_(0)

{
  warning("Reading DOLFIN xml meshes in parallel is deprecated. "
          "For better I/O performance, consider converting to flat binary");
}
//-----------------------------------------------------------------------------
PXMLMesh::~PXMLMesh()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
void PXMLMesh::startElement(const xmlChar *name, const xmlChar **attrs)
{
  switch (state_)
    {
    case OUTSIDE:

      if (xmlStrcasecmp(name, (xmlChar *) "mesh") == 0)
      {
        readMesh(name, attrs);
        state_ = INSIDE_MESH;
      }

      break;

    case INSIDE_MESH:

      if (xmlStrcasecmp(name, (xmlChar *) "vertices") == 0)
      {
        readVertices(name, attrs);
        state_ = INSIDE_VERTICES;
      }
      else if (xmlStrcasecmp(name, (xmlChar *) "cells") == 0)
      {
        readCells(name, attrs);
        state_ = INSIDE_CELLS;
      }
      else if (xmlStrcasecmp(name, (xmlChar *) "data") == 0)
      {
        state_ = INSIDE_DATA;
      }

      break;

    case INSIDE_VERTICES:

      if (xmlStrcasecmp(name, (xmlChar *) "vertex") == 0)
      {
        readVertex(name, attrs);
      }
      break;

    case INSIDE_CELLS:

      if (xmlStrcasecmp(name, (xmlChar *) "interval") == 0)
      {
        readInterval(name, attrs);
      }
      else if (xmlStrcasecmp(name, (xmlChar *) "triangle") == 0)
      {
        readTriangle(name, attrs);
      }
      else if (xmlStrcasecmp(name, (xmlChar *) "tetrahedron") == 0)
      {
        readTetrahedron(name, attrs);
      }

      break;

    case INSIDE_DATA:

      if (xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0)
      {
        error("Parsing mesh function is not supported in parallel");
        readMeshFunction(name, attrs);
        state_ = INSIDE_MESH_FUNCTION;
      }
      if (xmlStrcasecmp(name, (xmlChar *) "array") == 0)
      {
        error("Parsing array is not supported in parallel");
        readArray(name, attrs);
        state_ = INSIDE_ARRAY;
      }

      break;

    case INSIDE_MESH_FUNCTION:

      if (xmlStrcasecmp(name, (xmlChar *) "entity") == 0)
      {
        readMeshEntity(name, attrs);
      }

      break;

    case INSIDE_ARRAY:

      if (xmlStrcasecmp(name, (xmlChar *) "element") == 0)
      {
        readArrayElement(name, attrs);
      }
      break;

    default:
      break;
    }
}
//-----------------------------------------------------------------------------
void PXMLMesh::endElement(const xmlChar *name)
{
  switch (state_)
    {
    case INSIDE_MESH:

      if (xmlStrcasecmp(name, (xmlChar *) "mesh") == 0)
      {
        closeMesh();
        state_ = DONE;
      }

      break;

    case INSIDE_VERTICES:

      if (xmlStrcasecmp(name, (xmlChar *) "vertices") == 0)
      {
        state_ = INSIDE_MESH;
      }

      break;

    case INSIDE_CELLS:

      if (xmlStrcasecmp(name, (xmlChar *) "cells") == 0)
      {
        state_ = INSIDE_MESH;
      }

      break;

    case INSIDE_DATA:

      if (xmlStrcasecmp(name, (xmlChar *) "data") == 0)
      {
        state_ = INSIDE_MESH;
      }

      break;

    case INSIDE_MESH_FUNCTION:

      if (xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0)
      {
        state_ = INSIDE_DATA;
      }

      break;

    case INSIDE_ARRAY:

      if (xmlStrcasecmp(name, (xmlChar *) "array") == 0)
      {
        state_ = INSIDE_DATA;
      }

      break;

    default:
      break;
    }
}
//-----------------------------------------------------------------------------
void PXMLMesh::open(std::string filename)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool PXMLMesh::close()
{
  return state_ == DONE;
}
//-----------------------------------------------------------------------------
void PXMLMesh::readMesh(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  std::string type = parseString(name, attrs, "celltype");
  uint gdim = parseUnsignedInt(name, attrs, "dim");

  // Open mesh for editing
  if (editor_ != NULL)
  {
    error("In PXMLMesh, mesh editor is already created.");
  }
  CellType * cell_type = CellType::create(type);
  editor_ = new MeshEditor(mesh_, cell_type->cellType(), gdim);
  delete cell_type;
}
//-----------------------------------------------------------------------------
void PXMLMesh::readVertices(const xmlChar *name, const xmlChar **attrs)
{

  numParsedVertices_ = 0;
  numParsedCells_ = 0;

  // Parse values
  uint num_vertices = parseUnsignedInt(name, attrs, "size");

  uint pe_size = MPI::numProcesses();
  uint rank = MPI::processNumber();

  // Calculate a linear data distribution
  uint L = std::floor((real) num_vertices / (real) pe_size);
  uint R = num_vertices % pe_size;
  numLocalVertices_ = (num_vertices + pe_size - rank - 1) / pe_size;

  startIndex_vert_ = rank * L + std::min(rank, R);
  endIndex_vert_ = startIndex_vert_ + (numLocalVertices_ - 1);
  mesh_.distdata().set_num_global(0, num_vertices);

  // Set number of vertices
  editor_->initVertices(numLocalVertices_);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readCells(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  //  uint num_cells = parseUnsignedInt(name, attrs, "size");

  numParsedCells_ = 0;

  //FIXME
  editor_->initCells(1);
  editor_->close();
  MeshFunction<uint> pre_partition;
  mesh_.partition_geom(pre_partition);
  mesh_.distribute(pre_partition);

  for (VertexIterator vertex(mesh_); !vertex.end(); ++vertex)
  {
    own_vertex_[mesh_.distdata().get_global(*vertex)] = true;
  }
}
//-----------------------------------------------------------------------------
void PXMLMesh::readVertex(const xmlChar *name, const xmlChar **attrs)
{
  // Read index
  uint v = parseUnsignedInt(name, attrs, "index");

  if (v < startIndex_vert_ || v > endIndex_vert_) return;

  mesh_.distdata().set_map(numParsedVertices_, v, 0);

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
  editor_->addVertex(numParsedVertices_, &x[0]);
  numParsedVertices_++;
}
//-----------------------------------------------------------------------------
void PXMLMesh::readInterval(const xmlChar *name, const xmlChar **attrs)
{
  // Check dimension
  if (mesh_.topology().dim() != 1)
  {
    error("Mesh entity (interval) does not match dimension of mesh (%d).",
          mesh_.topology().dim());
  }

  // Ignore index
  uint v[2];
  v[0] = parseUnsignedInt(name, attrs, "v0");
  v[1] = parseUnsignedInt(name, attrs, "v1");

  // Return if no vertices are local
  if (!(own_vertex_[v[1]] || own_vertex_[v[0]]) || !own_vertex_[v[0]])
  {
    return;
  }

  used_vertex_[mesh_.distdata().get_vertex_local(v[0])] = true;
  if (own_vertex_[v[1]])
  {
    used_vertex_[mesh_.distdata().get_vertex_local(v[1])] = true;
  }

  // Add shared vertices to shared list
  if (!(own_vertex_[v[1]] && own_vertex_[v[0]]))
  {
    if (!own_vertex_[v[1]]) shared_buffer_[v[1]] = true;
  }

  // Add cell to cell buffer
  cell_buffer_.push_back(v[0]);
  cell_buffer_.push_back(v[1]);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readTriangle(const xmlChar *name, const xmlChar **attrs)
{
  // Check dimension
  if (mesh_.topology().dim() != 2)
  {
    error("Mesh entity (triangle) does not match dimension of mesh (%d).",
          mesh_.topology().dim());
  }

  // Ignore index
  uint v[3];
  v[0] = parseUnsignedInt(name, attrs, "v0");
  v[1] = parseUnsignedInt(name, attrs, "v1");
  v[2] = parseUnsignedInt(name, attrs, "v2");

  // Return if no vertices are local
  if (!(own_vertex_[v[1]] || own_vertex_[v[2]] || own_vertex_[v[0]])
      || !own_vertex_[v[0]])
  {
    return;
  }

  used_vertex_[mesh_.distdata().get_vertex_local(v[0])] = true;
  if (own_vertex_[v[1]])
  {
    used_vertex_[mesh_.distdata().get_vertex_local(v[1])] = true;
  }
  if (own_vertex_[v[2]])
  {
    used_vertex_[mesh_.distdata().get_vertex_local(v[2])] = true;
  }

  // Add shared vertices to shared list
  if (!(own_vertex_[v[1]] && own_vertex_[v[2]] && own_vertex_[v[0]]))
  {
    if (!own_vertex_[v[1]]) shared_buffer_[v[1]] = true;
    if (!own_vertex_[v[2]]) shared_buffer_[v[2]] = true;
  }

  // Add cell to cell buffer
  cell_buffer_.push_back(v[0]);
  cell_buffer_.push_back(v[1]);
  cell_buffer_.push_back(v[2]);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readTetrahedron(const xmlChar *name, const xmlChar **attrs)
{
  // Check dimension
  if (mesh_.topology().dim() != 3)
  {
    error("Mesh entity (tetrahedron) does not match dimension of mesh (%d).",
          mesh_.topology().dim());
  }

  // Ignore index
  uint v[4];
  v[0] = parseUnsignedInt(name, attrs, "v0");
  v[1] = parseUnsignedInt(name, attrs, "v1");
  v[2] = parseUnsignedInt(name, attrs, "v2");
  v[3] = parseUnsignedInt(name, attrs, "v3");

  // Return if no vertices are local
  if (!(own_vertex_[v[1]] || own_vertex_[v[2]] || own_vertex_[v[0]]
      || own_vertex_[v[3]]) || !own_vertex_[v[0]])
  {
    return;
  }

  used_vertex_[mesh_.distdata().get_vertex_local(v[0])] = true;
  if (own_vertex_[v[1]])
  {
    used_vertex_[mesh_.distdata().get_vertex_local(v[1])] = true;
  }
  if (own_vertex_[v[2]])
  {
    used_vertex_[mesh_.distdata().get_vertex_local(v[2])] = true;
  }
  if (own_vertex_[v[3]])
  {
    used_vertex_[mesh_.distdata().get_vertex_local(v[3])] = true;
  }

  // Add shared vertices to shared list
  if (!(own_vertex_[v[1]] && own_vertex_[v[2]] && own_vertex_[v[0]]
      && own_vertex_[v[3]]))
  {
    if (!own_vertex_[v[1]]) shared_buffer_[v[1]] = true;
    if (!own_vertex_[v[2]]) shared_buffer_[v[2]] = true;
    if (!own_vertex_[v[3]]) shared_buffer_[v[3]] = true;
  }

  // Add cell to cell buffer
  cell_buffer_.push_back(v[0]);
  cell_buffer_.push_back(v[1]);
  cell_buffer_.push_back(v[2]);
  cell_buffer_.push_back(v[3]);

}
//-----------------------------------------------------------------------------
void PXMLMesh::readMeshFunction(const xmlChar* name, const xmlChar** attrs)
{
  // Parse values
  const std::string id = parseString(name, attrs, "name");
  const std::string type = parseString(name, attrs, "type");
  uint const dim = parseUnsignedInt(name, attrs, "dim");
  uint const size = parseUnsignedInt(name, attrs, "size");

  // Only uint supported at this point
  if (strcmp(type.c_str(), "uint") != 0)
  {
    error("Only uint-valued mesh data is currently supported.");
  }

  // Check size
  mesh_.init(dim);
  if (mesh_.size(dim) != size)
  {
    error("Wrong number of values for MeshFunction, expecting %d.",
          mesh_.size(dim));
  }

  // Register data
  f_ = mesh_.data().createMeshFunction(id);
  dolfin_assert(f_);
  f_->init(mesh_, dim);

  // Set all values to zero
  *f_ = 0;
}
//-----------------------------------------------------------------------------
void PXMLMesh::readArray(const xmlChar* name, const xmlChar** attrs)
{
  // Parse values
  const std::string id = parseString(name, attrs, "name");
  const std::string type = parseString(name, attrs, "type");
  uint const size = parseUnsignedInt(name, attrs, "size");

  // Only uint supported at this point
  if (strcmp(type.c_str(), "uint") != 0)
  {
    error("Only uint-valued mesh data is currently supported.");
  }

  // Register data
  a_ = mesh_.data().createArray(id, size);
  dolfin_assert(a_);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readMeshEntity(const xmlChar* name, const xmlChar** attrs)
{
  // Read index
  uint const index = parseUnsignedInt(name, attrs, "index");

  // Read and set value
  dolfin_assert(f_);
  dolfin_assert(index < f_->size());
  uint const value = parseUnsignedInt(name, attrs, "value");
  f_->set(index, value);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readArrayElement(const xmlChar* name, const xmlChar** attrs)
{
  // Read index
  uint const index = parseUnsignedInt(name, attrs, "index");

  // Read and set value
  dolfin_assert(a_);
  dolfin_assert(index < a_->size());
  uint const value = parseUnsignedInt(name, attrs, "value");
  (*a_)[index] = value;
}
//-----------------------------------------------------------------------------
void PXMLMesh::closeMesh()
{
  Mesh new_mesh;
  delete editor_;
  editor_ = new MeshEditor(new_mesh, mesh_.type().cellType(),
                           mesh_.topology().dim(), mesh_.geometry().dim());

  int rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  Array<uint> send_buff, send_buff_indices, send_buff_orphan;
  // Make room for own vertices and shared ones
  uint shared = 0;
  uint orphan = 0;
  for (uint i = 0; i < shared_buffer_.size(); i++)
    if (shared_buffer_[i])
    {
      shared++;
      send_buff.push_back(i);
    }

  _map<uint,bool> assigned_orphan, ghost_vertex;
  for (uint i = 0; i < mesh_.numVertices(); i++)
  {
    if (!used_vertex_[i])
    {
      orphan++;
    }
  }

  int num_recv, max_nsh;
  Array<double> send_buff_coords;
  MPI_Status status;
  uint gdim = mesh_.geometry().dim();
  uint num_shared = shared;
  uint num_coords = gdim * num_shared;
  uint num_orphan = num_shared;
  real *shared_coords = new real[num_coords];
  real *shp = &shared_coords[0];
  uint *shared_indices = new uint[num_shared];
  uint *shpi = &shared_indices[0];
  uint *shared_orphans = new uint[num_shared];
  uint *shpo = &shared_orphans[0];
  MPI_Allreduce(&num_shared, &max_nsh, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  uint *shvert = new uint[max_nsh];
  uint src, dest;
  _map<uint, uint> owner_map;
  // Exchange ghost points
  for (uint j = 1; j < pe_size; j++)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, dest, 1, shvert,
                 max_nsh, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &num_recv);

    for (int k = 0; k < num_recv; k++)
      if (own_vertex_[shvert[k]])
      {
        Vertex vertex(mesh_, mesh_.distdata().get_vertex_local(shvert[k]));
        send_buff_coords.push_back(vertex.point().x());
        send_buff_coords.push_back(vertex.point().y());
        if (gdim > 2) send_buff_coords.push_back(vertex.point().z());

        send_buff_indices.push_back(shvert[k]);
        if (!used_vertex_[vertex.index()] && !assigned_orphan[vertex.index()])
        {
          owner_map[vertex.index()] = status.MPI_SOURCE;
          send_buff_orphan.push_back(pe_size);
          assigned_orphan[vertex.index()] = true;
        }
        else
        {
          if (owner_map.count(vertex.index()) == 0) send_buff_orphan.push_back(
              rank);
          else send_buff_orphan.push_back(owner_map[vertex.index()]);
        }
        ghost_vertex[mesh_.distdata().get_vertex_local(shvert[k])] = true;
      }

    MPI_Sendrecv(&send_buff_indices[0], send_buff_indices.size(), MPI_UNSIGNED,
                 src, 1, shpi, num_shared, MPI_UNSIGNED, dest, 1,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &num_recv);
    num_shared -= num_recv;
    shpi += num_recv;

    MPI_Sendrecv(&send_buff_coords[0], send_buff_coords.size(), MPI_DOUBLE, src,
                 2, shp, num_coords, MPI_DOUBLE, dest, 2, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_DOUBLE, &num_recv);
    num_coords -= num_recv;
    shp += num_recv;

    MPI_Sendrecv(&send_buff_orphan[0], send_buff_orphan.size(), MPI_UNSIGNED,
                 src, 3, shpo, num_orphan, MPI_UNSIGNED, dest, 3,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &num_recv);
    num_orphan -= num_recv;
    shpo += num_recv;

    send_buff_coords.clear();
    send_buff_indices.clear();
    send_buff_orphan.clear();

  }

  // Init new mesh
  editor_->initVertices(mesh_.numVertices() + shared - orphan);
  //  new_mesh.distdata().init(_mesh.numVertices() + shared - orphan);
  new_mesh.distdata().set_num_global(0, mesh_.global_numVertices());

  uint v = 0;
  for (VertexIterator vertex(mesh_); !vertex.end(); ++vertex)
  {
    if (used_vertex_[vertex->index()])
    {
      editor_->addVertex(v, vertex->point());
      new_mesh.distdata().set_map(v, mesh_.distdata().get_global(*vertex), 0);
      if (ghost_vertex[vertex->index()])
      {
        new_mesh.distdata().set_shared(v, 0);
      }
      v++;
    }
  }

  //Add shared ghost vertices
  uint ci = 0;
  for (uint i = 0; i < shared; i++)
  {
    new_mesh.distdata().set_map(v, shared_indices[i], 0);
    if (shared_orphans[i] < pe_size)
    {  // Why...ugly hack to set ghost owner
      new_mesh.distdata().set_ghost(v, 0);
      new_mesh.distdata().set_ghost_owner(v, shared_orphans[i], 0);
    }
    editor_->addVertex(v++, &shared_coords[ci]);
    ci += gdim;
  }

  uint ndims = mesh_.type().numVertices(mesh_.topology().dim());
  editor_->initCells(cell_buffer_.size() / ndims);
  uint c = 0;
  uint * connectivity = new uint[ndims];
  for (uint i = 0; i < cell_buffer_.size(); i += ndims)
  {
    for (uint n = 0; n < ndims; ++n)
    {
      connectivity[n] = new_mesh.distdata().get_vertex_local(
          cell_buffer_[i + n]);
    }
    editor_->addCell(c++, &connectivity[0]);
  }
  delete[] connectivity;
  editor_->close();
  mesh_ = new_mesh;

  ghost_vertex.clear();
  assigned_orphan.clear();
  shared_buffer_.clear();
  send_buff.clear();
  cell_buffer_.clear();
  used_vertex_.clear();
  own_vertex_.clear();
  delete[] shared_indices;
  delete[] shared_coords;
  delete[] shared_orphans;
  delete[] shvert;
  delete editor_;
  editor_ = NULL;
}
//-----------------------------------------------------------------------------
#else
void PXMLMesh::startElement(const xmlChar* name, const xmlChar** attrs)
{
}
//-----------------------------------------------------------------------------
void PXMLMesh::endElement(const xmlChar* name)
{
}
//-----------------------------------------------------------------------------
void PXMLMesh::open(std::string filename)
{
}
//-----------------------------------------------------------------------------
bool PXMLMesh::close()
{
  return false;
}
//-----------------------------------------------------------------------------
#endif
#endif /* HAVE_XML */

}

