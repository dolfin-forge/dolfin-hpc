// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-16
// Last changed: 2008-05-19

#include <dolfin/mesh/MeshEditor.h>

#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/Space.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshEditor::MeshEditor(Mesh& mesh, CellType const& type, Space const& space) :
    mesh_(&mesh),
    cell_vertices_(NULL),
    tdim_(0),
    gdim_(0),
    num_vertices_(0),
    num_cells_(0),
    vertex_index_(0),
    cell_index_(0),
    open_(false)
{
  init(mesh, type, space);
}
//-----------------------------------------------------------------------------
MeshEditor::MeshEditor(Mesh& mesh, CellType::Type cell_type, uint gdim) :
    mesh_(&mesh),
    cell_vertices_(NULL),
    tdim_(0),
    gdim_(0),
    num_vertices_(0),
    num_cells_(0),
    vertex_index_(0),
    cell_index_(0),
    open_(false)
{
  CellType * type = CellType::create(cell_type);
  EuclideanSpace space(gdim);
  init(mesh, *type, space);
  delete type;
}
//-----------------------------------------------------------------------------
MeshEditor::MeshEditor(Mesh& mesh) :
    mesh_(&mesh),
    cell_vertices_(NULL),
    tdim_(0),
    gdim_(0),
    num_vertices_(0),
    num_cells_(0),
    vertex_index_(0),
    cell_index_(0),
    open_(false)
{
  if (mesh.empty()) { error("MeshEditor : provided mesh is empty"); }
  init(mesh, mesh.type(), mesh.space());
}
//-----------------------------------------------------------------------------
MeshEditor::~MeshEditor()
{
  if(open_)
  {
    error("MeshEditor : editor has not been closed before destruction");
  }
}
//-----------------------------------------------------------------------------
void MeshEditor::init(Mesh& mesh, CellType const& type, Space const& space)
{
  // Save mesh and dimension
  this->tdim_ = type.dim();
  this->gdim_ = space.dim();

  // Initialize the topology to the given cell type and space
  mesh.init(type, space);

  open_ = true;
}
//-----------------------------------------------------------------------------
void MeshEditor::init_vertices(uint num_local, uint num_global /* = 0 */)
{
  if(!open_)
  {
    error("MeshEditor : initializing vertices on empty editor");
  }
  // Initialize mesh data
  this->num_vertices_ = num_local;
  mesh_->topology_.init(0, num_local, num_global);
  mesh_->geometry_.init(mesh_->space(), num_local);
}
//-----------------------------------------------------------------------------
void MeshEditor::init_cells(uint num_local, uint num_global /* = 0 */)
{
  if(!open_)
  {
    error("MeshEditor : initializing cells on empty editor");
  }
  // Initialize mesh data
  this->num_cells_ = num_local;
  mesh_->topology_.init(tdim_, num_local, num_global);

  // Create a shortcut to cell vertices connectivity to avoid checking its
  // existence at every cell creation
  this->cell_vertices_ = &mesh_->topology_(tdim_, 0);
}
//-----------------------------------------------------------------------------
void MeshEditor::add_vertex(uint v, real const * x)
{
  if (v >= num_vertices_)
  {
    error("Vertex index (%d) out of range [0, %d].", v, num_vertices_ - 1);
  }
  if (vertex_index_ >= num_vertices_)
  {
    error("MeshEditor : vertex list full, %d vertices added.", num_vertices_);
  }
  mesh_->geometry_.set(v, x);
  ++vertex_index_;
}
//-----------------------------------------------------------------------------
void MeshEditor::add_cell(uint c, uint const * v)
{
  if (c >= num_cells_)
  {
   error("Cell index (%d) out of range [0, %d].", c, num_cells_ - 1);
  }
  if (cell_index_ >= num_cells_)
  {
   error("MeshEditor : cell list full, %d cells added.", num_cells_);
  }
  dolfin_assert(cell_vertices_ != NULL);
  cell_vertices_->set(c, v);
  ++cell_index_;
}
//-----------------------------------------------------------------------------
void MeshEditor::close()
{
  // Check consistency of number of vertices
  if( this->num_vertices_ != mesh_->topology().size(0))
  {
    error("Mismatch between number of vertices initialized and added to mesh : "
          "%d != %d", this->num_vertices_, mesh_->topology().size(0));
  }
  // Check consistency of number of cells
  if( this->num_cells_ != mesh_->topology().size(tdim_))
  {
    error("Mismatch between number of cells initialized and added to mesh : "
          "%d != %d", this->num_cells_, mesh_->topology().size(tdim_));
  }
  // Finalize topology and geometry
  mesh_->topology_.finalize();
  mesh_->geometry_.finalize();
  // Clear data
  clear();
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
  open_ = false;
}
//-----------------------------------------------------------------------------
uint MeshEditor::current_vertex() const
{
  return vertex_index_;
}
//-----------------------------------------------------------------------------
uint MeshEditor::current_cell() const
{
  return cell_index_;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
