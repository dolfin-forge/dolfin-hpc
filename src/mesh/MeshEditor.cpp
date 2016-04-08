// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-16
// Last changed: 2008-05-19

#include <dolfin/mesh/MeshEditor.h>

#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshEditor::MeshEditor(Mesh& mesh, CellType const& cell_type, uint gdim) :
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
  init(mesh, cell_type, gdim);
}
//-----------------------------------------------------------------------------
MeshEditor::MeshEditor(Mesh& mesh, CellType::Type type, uint gdim) :
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
  CellType * cell_type = CellType::create(type);
  init(mesh, *cell_type, gdim);
  delete cell_type;
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
void MeshEditor::init(Mesh& mesh, CellType const& type, uint gdim)
{
  // Clear old mesh data
  mesh.clear();

  // Set cell type
  mesh.cell_type_ = type.clone();

  // Save mesh and dimension
  this->tdim_ = mesh.cell_type_->dim();
  this->gdim_ = gdim;

  // Initialize the topology to the given topological dimension
  mesh_->topology_.init(tdim_, MPI::numProcesses() > 1);

  // Create a shortcut to cell vertices connectivity to avoid checking its
  // existence at every cell creation
  this->cell_vertices_ = &mesh_->topology_(tdim_, 0);

  open_ = true;
}
//-----------------------------------------------------------------------------
void MeshEditor::init_vertices(uint num_local, uint num_global /* = 0 */)
{
  if(cell_vertices_ == NULL)
  {
    error("MeshEditor : initializing vertices on empty editor");
  }
  // Initialize mesh data
  this->num_vertices_ = num_local;
  mesh_->topology_.init(0    , num_local, num_global);
  mesh_->geometry_.init(gdim_, num_local);
}
//-----------------------------------------------------------------------------
void MeshEditor::init_cells(uint num_local, uint num_global /* = 0 */)
{
  if(cell_vertices_ == NULL)
  {
    error("MeshEditor : initializing cells on empty editor");
  }
  // Initialize mesh data
  this->num_cells_ = num_local;
  mesh_->topology_.init(tdim_, num_local, num_global);
}
//-----------------------------------------------------------------------------
void MeshEditor::add_cells(Array<Array<uint> > const& connectivity,
                           uint num_global /* = 0 */)
{
  // Initialize mesh data
  this->num_cells_ = connectivity.size();
  cell_vertices_->set(connectivity);
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
  mesh_->topology_(tdim_, 0).set(c, v);
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
