// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-10-21
// Last changed: 2008-05-21

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_XML

#include <cstring>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/io/XMLMesh.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
XMLMesh::XMLMesh(Mesh& mesh) :
    XMLObject(),
    mesh_(mesh),
    state_(OUTSIDE),
    editor_(NULL),
    f_(NULL),
    a_(NULL)
{
  // Do nothing
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

      if (xmlStrcasecmp(name, (xmlChar *) "vertex") == 0) readVertex(name,
                                                                     attrs);

      break;

    case INSIDE_CELLS:

      if (xmlStrcasecmp(name, (xmlChar *) "interval") == 0) readInterval(name,
                                                                         attrs);
      else if (xmlStrcasecmp(name, (xmlChar *) "triangle") == 0) readTriangle(
          name, attrs);
      else if (xmlStrcasecmp(name, (xmlChar *) "tetrahedron") == 0) readTetrahedron(
          name, attrs);

      break;

    case INSIDE_DATA:

      if (xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0)
      {
        readMeshFunction(name, attrs);
        state_ = INSIDE_MESH_FUNCTION;
      }
      if (xmlStrcasecmp(name, (xmlChar *) "array") == 0)
      {
        readArray(name, attrs);
        state_ = INSIDE_ARRAY;
      }

      break;

    case INSIDE_MESH_FUNCTION:

      if (xmlStrcasecmp(name, (xmlChar *) "entity") == 0) readMeshEntity(name,
                                                                         attrs);

      break;

    case INSIDE_ARRAY:

      if (xmlStrcasecmp(name, (xmlChar *) "element") == 0) readArrayElement(
          name, attrs);

      break;

    default:
      ;
    }
}
//-----------------------------------------------------------------------------
void XMLMesh::endElement(const xmlChar *name)
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
      ;
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
  return state_ == DONE;
}
//-----------------------------------------------------------------------------
void XMLMesh::readMesh(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  std::string type = parseString(name, attrs, "celltype");
  uint gdim = parseUnsignedInt(name, attrs, "dim");
  
  // Create cell type to get topological dimension
  CellType* cell_type = CellType::create(type);
  uint tdim = cell_type->dim();
  delete cell_type;

  // Open mesh for editing
  if(editor_ != NULL)
  {
    error("In XMLMesh, mesh editor is already open.");
  }
  editor_ = new MeshEditor(mesh_, CellType::string2type(type), tdim, gdim);
}
//-----------------------------------------------------------------------------
void XMLMesh::readVertices(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  uint num_vertices = parseUnsignedInt(name, attrs, "size");

  // Set number of vertices
  editor_->initVertices(num_vertices);
}
//-----------------------------------------------------------------------------
void XMLMesh::readCells(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  uint num_cells = parseUnsignedInt(name, attrs, "size");

  // Set number of vertices
  editor_->initCells(num_cells);
}
//-----------------------------------------------------------------------------
void XMLMesh::readVertex(const xmlChar *name, const xmlChar **attrs)
{
  // Read index
  uint v = parseUnsignedInt(name, attrs, "index");
  
  // Handle differently depending on geometric dimension
  switch (mesh_.geometry().dim())
    {
    case 1:
      {
        real x = parseReal(name, attrs, "x");
        editor_->addVertex(v, x);
      }
      break;
    case 2:
      {
        real x = parseReal(name, attrs, "x");
        real y = parseReal(name, attrs, "y");
        editor_->addVertex(v, x, y);
      }
      break;
    case 3:
      {
        real x = parseReal(name, attrs, "x");
        real y = parseReal(name, attrs, "y");
        real z = parseReal(name, attrs, "z");
        editor_->addVertex(v, x, y, z);
      }
      break;
    default:
      error("Dimension of mesh must be 1, 2 or 3.");
    }
}
//-----------------------------------------------------------------------------
void XMLMesh::readInterval(const xmlChar *name, const xmlChar **attrs)
{
  // Check dimension
  if (mesh_.topology().dim() != 1) error(
      "Mesh entity (interval) does not match dimension of mesh (%d).",
      mesh_.topology().dim());

  // Parse values
  uint c = parseUnsignedInt(name, attrs, "index");
  uint v0 = parseUnsignedInt(name, attrs, "v0");
  uint v1 = parseUnsignedInt(name, attrs, "v1");
  
  // Add cell
  editor_->addCell(c, v0, v1);
}
//-----------------------------------------------------------------------------
void XMLMesh::readTriangle(const xmlChar *name, const xmlChar **attrs)
{
  // Check dimension
  if (mesh_.topology().dim() != 2) error(
      "Mesh entity (triangle) does not match dimension of mesh (%d).",
      mesh_.topology().dim());

  // Parse values
  uint c = parseUnsignedInt(name, attrs, "index");
  uint v0 = parseUnsignedInt(name, attrs, "v0");
  uint v1 = parseUnsignedInt(name, attrs, "v1");
  uint v2 = parseUnsignedInt(name, attrs, "v2");
  
  // Add cell
  editor_->addCell(c, v0, v1, v2);
}
//-----------------------------------------------------------------------------
void XMLMesh::readTetrahedron(const xmlChar *name, const xmlChar **attrs)
{
  // Check dimension
  if (mesh_.topology().dim() != 3) error(
      "Mesh entity (tetrahedron) does not match dimension of mesh (%d).",
      mesh_.topology().dim());

  // Parse values
  uint c = parseUnsignedInt(name, attrs, "index");
  uint v0 = parseUnsignedInt(name, attrs, "v0");
  uint v1 = parseUnsignedInt(name, attrs, "v1");
  uint v2 = parseUnsignedInt(name, attrs, "v2");
  uint v3 = parseUnsignedInt(name, attrs, "v3");
  
  // Add cell
  editor_->addCell(c, v0, v1, v2, v3);
}
//-----------------------------------------------------------------------------
void XMLMesh::readMeshFunction(const xmlChar* name, const xmlChar** attrs)
{
  // Parse values
  const std::string id = parseString(name, attrs, "name");
  const std::string type = parseString(name, attrs, "type");
  const uint dim = parseUnsignedInt(name, attrs, "dim");
  const uint size = parseUnsignedInt(name, attrs, "size");

  // Only uint supported at this point
  if (strcmp(type.c_str(), "uint") != 0) error(
      "Only uint-valued mesh data is currently supported.");

  // Check size
  mesh_.init(dim);
  if (mesh_.size(dim) != size) error(
      "Wrong number of values for MeshFunction, expecting %d.",
      mesh_.size(dim));

  // Register data
  f_ = mesh_.data().createMeshFunction(id);
  dolfin_assert(f_);
  f_->init(mesh_, dim);

  // Set all values to zero
  *f_ = 0;
}
//-----------------------------------------------------------------------------
void XMLMesh::readArray(const xmlChar* name, const xmlChar** attrs)
{
  // Parse values
  const std::string id = parseString(name, attrs, "name");
  const std::string type = parseString(name, attrs, "type");
  const uint size = parseUnsignedInt(name, attrs, "size");

  // Only uint supported at this point
  if (strcmp(type.c_str(), "uint") != 0) error(
      "Only uint-valued mesh data is currently supported.");

  // Register data
  a_ = mesh_.data().createArray(id, size);
  dolfin_assert(a_);
}
//-----------------------------------------------------------------------------
void XMLMesh::readMeshEntity(const xmlChar* name, const xmlChar** attrs)
{
  // Read index
  const uint index = parseUnsignedInt(name, attrs, "index");

  // Read and set value
  dolfin_assert(f_); dolfin_assert(index < f_->size());
  const uint value = parseUnsignedInt(name, attrs, "value");
  f_->set(index, value);
}
//-----------------------------------------------------------------------------
void XMLMesh::readArrayElement(const xmlChar* name, const xmlChar** attrs)
{
  // Read index
  const uint index = parseUnsignedInt(name, attrs, "index");

  // Read and set value
  dolfin_assert(a_); dolfin_assert(index < a_->size());
  const uint value = parseUnsignedInt(name, attrs, "value");
  (*a_)[index] = value;
}
//-----------------------------------------------------------------------------
void XMLMesh::closeMesh()
{
  editor_->close();
  delete editor_;
  editor_ = NULL;
}
//-----------------------------------------------------------------------------

#endif
