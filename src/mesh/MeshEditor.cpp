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
    tdim_(0),
    gdim_(0),
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
  mesh.cell_type_ = CellType::create(type);

  // Save mesh and dimension
  this->tdim_ = mesh.cell_type_->dim();
  this->gdim_ = gdim;

  // Initialize topological dimension
  mesh.topology_.init(tdim_);
}
//-----------------------------------------------------------------------------
void MeshEditor::initVertices(uint num_vertices)
{
  // Initialize mesh data
  this->num_vertices_ = num_vertices;
  mesh_->topology_.init(0, num_vertices);
  mesh_->geometry_.init(gdim_, num_vertices);
}
//-----------------------------------------------------------------------------
void MeshEditor::initCells(uint num_cells)
{
  // Initialize mesh data
  this->num_cells_ = num_cells;
  mesh_->topology_.init(tdim_, num_cells);
  mesh_->topology_(tdim_, 0).init(num_cells, mesh_->type().numVertices(tdim_));
}
//-----------------------------------------------------------------------------
void MeshEditor::addVertex(uint v, Point const& p)
{
  // Add vertex
  addVertexCommon(v);

  // Set coordinate
  mesh_->geometry_.set(v, &p[0]);
}
//-----------------------------------------------------------------------------
void MeshEditor::addVertex(uint v, real const * x)
{
  // Add vertex
  addVertexCommon(v);

  // Set coordinate
  mesh_->geometry_.set(v, x);
}
//-----------------------------------------------------------------------------
void MeshEditor::addCell(uint c, const Array<uint>& v)
{
  // Add cell
  addCellCommon(c);

  // Set data
  mesh_->topology_(tdim_, 0).set(c, v);
}
//-----------------------------------------------------------------------------
void MeshEditor::addCell(uint c, uint const * v)
{
  // Add cell
  addCellCommon(c);

  // Set data
  mesh_->topology_(tdim_, 0).set(c, v);
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
}
//-----------------------------------------------------------------------------

}

