// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-16
// Last changed: 2008-05-19

#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/MeshEditor.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshEditor::MeshEditor(Mesh& mesh, CellType::Type type, uint gdim) :
    mesh_(&mesh),
    tdim_(0),
    gdim_(0),
    num_vertices_(0),
    num_cells_(0),
    vertex_index_(0),
    cell_index_(0)
{
  // Do nothing
  init(mesh, type, gdim);
}

//-----------------------------------------------------------------------------
MeshEditor::MeshEditor(Mesh& mesh, CellType::Type type, uint tdim, uint gdim) :
    mesh_(&mesh),
    tdim_(tdim),
    gdim_(gdim),
    num_vertices_(0),
    num_cells_(0),
    vertex_index_(0),
    cell_index_(0)
{
  CellType * t = CellType::create(type);
  if (tdim != t->dim())
  {
    error("In MeshEditor, cell type and topological dimension do not match.");
  }
  delete t;

  // Do nothing
  init(mesh, type, gdim);
}

//-----------------------------------------------------------------------------
MeshEditor::~MeshEditor()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void MeshEditor::init(Mesh& mesh, CellType::Type type, uint gdim)
{
  // Clear old mesh data
  mesh.clear();

  // Set cell type
  mesh._cell_type = CellType::create(type);

  // Save mesh and dimension
  this->tdim_ = mesh._cell_type->dim();

  // Initialize topological dimension
  mesh._topology.init(tdim_);

  // Initialize temporary storage for local cell data
  vertices.resize(mesh.type().numVertices(tdim_), 0.0);
}
//-----------------------------------------------------------------------------
void MeshEditor::initVertices(uint num_vertices)
{
  // Check if we are currently editing a mesh
  if (!mesh_)
  {
    error("No mesh opened, unable to edit.");
  }
  
  // Initialize mesh data
  this->num_vertices_ = num_vertices;
  mesh_->_topology.init(0, num_vertices);
  mesh_->_geometry.init(gdim_, num_vertices);
}
//-----------------------------------------------------------------------------
void MeshEditor::initCells(uint num_cells)
{
  // Check if we are currently editing a mesh
  if (!mesh_)
  {
    error("No mesh opened, unable to edit.");
  }

  // Initialize mesh data
  this->num_cells_ = num_cells;
  mesh_->_topology.init(tdim_, num_cells);
  mesh_->_topology(tdim_, 0).init(num_cells, mesh_->type().numVertices(tdim_));
}
//-----------------------------------------------------------------------------
void MeshEditor::addVertex(uint v, Point const& p)
{
  // Add vertex
  addVertexCommon(v);
  
  // Set coordinate
  mesh_->_geometry.set(v, &p[0]);
}
//-----------------------------------------------------------------------------
void MeshEditor::addVertex(uint v, real const * x)
{
  // Add vertex
  addVertexCommon(v);

  // Set coordinate
  mesh_->_geometry.set(v, x);
}
//-----------------------------------------------------------------------------
void MeshEditor::addVertex(uint v, real x)
{
  // Add vertex
  addVertexCommon(v);

  // Set coordinate
  mesh_->_geometry.set(v, 0, x);
}
//-----------------------------------------------------------------------------
void MeshEditor::addVertex(uint v, real x, real y)
{
  // Add vertex
  addVertexCommon(v);

  // Set coordinate
  mesh_->_geometry.set(v, 0, x);
  mesh_->_geometry.set(v, 1, y);
}
//-----------------------------------------------------------------------------
void MeshEditor::addVertex(uint v, real x, real y, real z)
{
  // Add vertex
  addVertexCommon(v);

  // Set coordinate
  mesh_->_geometry.set(v, 0, x);
  mesh_->_geometry.set(v, 1, y);
  mesh_->_geometry.set(v, 2, z);
}
//-----------------------------------------------------------------------------
void MeshEditor::addCell(uint c, const Array<uint>& v)
{
  // Add cell
  addCellCommon(c);

  // Set data
  mesh_->_topology(tdim_, 0).set(c, v);
}
//-----------------------------------------------------------------------------
void MeshEditor::addCell(uint c, uint v0, uint v1)
{
  dolfin_assert(mesh_->_cell_type->numEntities(0) == 2);

  // Add cell
  addCellCommon(c);

  // Set data
  vertices[0] = v0;
  vertices[1] = v1;
  mesh_->_topology(tdim_, 0).set(c, vertices);
}
//-----------------------------------------------------------------------------
void MeshEditor::addCell(uint c, uint v0, uint v1, uint v2)
{
  dolfin_assert(mesh_->_cell_type->numEntities(0) == 3);

  // Add cell
  addCellCommon(c);

  // Set data
  vertices[0] = v0;
  vertices[1] = v1;
  vertices[2] = v2;
  mesh_->_topology(tdim_, 0).set(c, vertices);
}
//-----------------------------------------------------------------------------
void MeshEditor::addCell(uint c, uint v0, uint v1, uint v2, uint v3)
{
  dolfin_assert(mesh_->_cell_type->numEntities(0) == 4);

  // Add cell
  addCellCommon(c);

  // Set data
  vertices[0] = v0;
  vertices[1] = v1;
  vertices[2] = v2;
  vertices[3] = v3;
  mesh_->_topology(tdim_, 0).set(c, vertices);
}
//-----------------------------------------------------------------------------
void MeshEditor::close()
{
  // Clear data
  clear();
}
//-----------------------------------------------------------------------------
void MeshEditor::addVertexCommon(uint v)
{
  // Check value of vertex index
  if (v >= num_vertices_)
  {
    error("Vertex index (%d) out of range [0, %d].", v, num_vertices_ - 1);
  }

  // Check if there is room for more vertices
  if (vertex_index_ >= num_vertices_)
  {
    error("Vertex list is full, %d vertices already specified.", num_vertices_);
  }
  
  // Step to next vertex
  ++vertex_index_;
}
//-----------------------------------------------------------------------------
void MeshEditor::addCellCommon(uint c)
{
  // Check value of cell index
  if (c >= num_cells_)
  {
    error("Cell index (%d) out of range [0, %d].", c, num_cells_ - 1);
  }

  // Check if there is room for more cells
  if (cell_index_ >= num_cells_)
  {
    error("Cell list is full, %d cells already specified.", num_cells_);
  }

  // Step to next cell
  ++cell_index_;
}
//-----------------------------------------------------------------------------
void MeshEditor::clear()
{
  tdim_ = 0;
  gdim_ = 0;
  num_vertices_ = 0;
  num_cells_ = 0;
  vertex_index_ = 0;
  cell_index_ = 0;
  vertices.clear();
}
//-----------------------------------------------------------------------------

}

