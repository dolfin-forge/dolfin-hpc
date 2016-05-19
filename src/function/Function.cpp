// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Dag Lindbo, 2008.
// Modified by Kristen Kaasbjerg, 2008.
// Modified by Niclas Jansson, 2008-2010.
// Modified by Aurélien Larcher 2013-2014.
//
// First added:  2007-04-02
// Last changed: 2014-02-06

#include <dolfin/function/Function.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/FunctionDecomposition.h>
#include <dolfin/function/FunctionInterpolation.h>
#include <dolfin/function/SubFunction.h>
#include <dolfin/la/Vector.h>

#include <algorithm>
#include <set>

namespace dolfin
{

//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh) :
    GenericFunction(),
    mesh_(&mesh),
    discrete_space_(NULL),
    element_(NULL),
    dofmap_(NULL),
    scratch(NULL),
    X_(NULL),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL),
    cache_mapping_(NULL)
{
  // Do nothing
}

//-----------------------------------------------------------------------------
Function::Function(Form& form, uint i) :
    GenericFunction(),
    mesh_(&form.dofmaps()[i].mesh()),
    discrete_space_(new FiniteElementSpace(form, i)),
    element_(&discrete_space_->element()),
    dofmap_(&discrete_space_->dofmap()),
    scratch(new ScratchSpace(*discrete_space_)),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL),
    cache_mapping_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
Function::Function(FiniteElementSpace const& space) :
    GenericFunction(),
    mesh_(&space.mesh()),
    discrete_space_(new FiniteElementSpace(space)),
    element_(&discrete_space_->element()),
    dofmap_(&discrete_space_->dofmap()),
    scratch(new ScratchSpace(*discrete_space_)),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL),
    cache_mapping_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, ufl::FiniteElementSpace const& finite_element) :
    GenericFunction(),
    mesh_(&mesh),
    discrete_space_(new FiniteElementSpace(mesh, finite_element)),
    element_(&discrete_space_->element()),
    dofmap_(&discrete_space_->dofmap()),
    scratch(new ScratchSpace(*discrete_space_)),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL),
    cache_mapping_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
Function::Function(SubFunction const& sub_function) :
    GenericFunction(),
    mesh_(&sub_function.function().mesh()),
    discrete_space_(new FiniteElementSpace(sub_function.function().space(),
                                           sub_function.index())),
    element_(&discrete_space_->element()),
    dofmap_(&discrete_space_->dofmap()),
    scratch(new ScratchSpace(*discrete_space_)),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL),
    cache_mapping_(NULL)
{
  // Initialize vector, scratch space and ghosts
  InitializeVector();

  // Copy subvector, naive implementation
  Function& gFunc = sub_function.function();
  DofMap const& gDm = gFunc.space().dofmap();
  uint const gLocalDim = gDm.local_dimension();
  uint const gDmOffset = gDm.sub_dofmaps_offsets()[sub_function.index()];
  uint const thisLocalDim = scratch->local_dimension;

  // Sync ghosts before getting the block
  real * gblock = gFunc.create_block();
  gFunc.sync_ghosts();
  gFunc.get_block(gblock);

  // Loop baby, loop...
  uint gBlockOffset = 0;
  uint ii = 0;
  real * this_block = this->create_block();
  for (CellIterator cell(*mesh_); !cell.end(); ++cell, gBlockOffset += gLocalDim)
  {
    for (uint dof = 0; dof < thisLocalDim; ++dof)
    {
      this_block[ii++] = gblock[gBlockOffset + gDmOffset + dof];
    }
  }
  this->set_block(this_block);

  delete[] this_block;
  delete[] gblock;
}

//-----------------------------------------------------------------------------
Function::Function(Function const& other) :
    GenericFunction(),
    mesh_(&other.mesh()),
    discrete_space_(NULL),
    element_(NULL),
    dofmap_(NULL),
    scratch(NULL),
    X_(NULL),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL),
    cache_mapping_(NULL)
{
  if(!other.empty())
  {
    *this = other;
    this->sync_ghosts();
  }
}

//-----------------------------------------------------------------------------
Function::~Function()
{
  clear();
}

//-----------------------------------------------------------------------------
bool Function::empty() const
{
  return (discrete_space_ == NULL);
}

//-----------------------------------------------------------------------------
void Function::init(Form& form, uint i)
{
  if(mesh_ != &form.dofmaps()[i].mesh())
  {
    error("Function : mesh mismatch between function and coefficient %d", i);
  }
  //
  clear();
  discrete_space_ = new FiniteElementSpace(form, i);
  element_ = &discrete_space_->element();
  dofmap_ = &discrete_space_->dofmap();
  scratch = new ScratchSpace(*discrete_space_);
  X_ = new Vector();
  //
  InitializeVector();
}

//-----------------------------------------------------------------------------
void Function::init(FiniteElementSpace const& space)
{
  if(mesh_ != &space.mesh())
  {
    error("Function : mesh mismatch between function and space");
  }
  //
  clear();
  discrete_space_ = new FiniteElementSpace(space);
  element_ = &discrete_space_->element();
  dofmap_ = &discrete_space_->dofmap();
  scratch = new ScratchSpace(*discrete_space_);
  X_ = new Vector();
  //
  InitializeVector();
}

//-----------------------------------------------------------------------------
void Function::clear()
{
  delete X_;
  X_ = NULL;
  delete discrete_space_;
  discrete_space_ = NULL;
  element_ = NULL;
  dofmap_ = NULL;
  delete scratch;
  scratch = NULL;
  delete[] indices_;
  indices_ = NULL;
  delete[] data_cache_;
  data_cache_ = NULL;
  delete cache_mapping_;
  cache_mapping_ = NULL;
  renumbered_ = true;
}

//--- UFC INTERFACE -----------------------------------------------------------
void Function::evaluate(real* values, const real* x,
                        const ufc::cell& cell) const
{
  UFCCell const * ufc_cell = static_cast<UFCCell const *>(&cell);

  // Get expansion coefficients on cell
  dofmap_->tabulate_dofs(scratch->dofs, *ufc_cell);
  X_->get(scratch->coefficients, scratch->local_dimension, scratch->dofs);

  // Compute linear combination
  for (uint j = 0; j < scratch->size; ++j)
  {
    values[j] = 0.0;
  }
  for (uint i = 0; i < element_->space_dimension(); ++i)
  {
    element_->evaluate_basis(i, scratch->values, x, *ufc_cell);
    for (uint j = 0; j < scratch->size; ++j)
    {
      values[j] += scratch->coefficients[i] * scratch->values[j];
    }
  }
}

//--- GenericFunction ---------------------------------------------------------
Mesh& Function::mesh() const
{
  return (*mesh_);
}

//-----------------------------------------------------------------------------
void Function::eval(real* values, const real* x) const
{
  // Find the cell that contains x
  Point p(x, mesh_->geometry().dim());
  Array<uint> cells;
  mesh_->intersector().overlap(p, cells);
  if (cells.size() < 1)
  {
    if (!mesh_->is_distributed())
    {
      error("Unable to evaluate function at given point (not inside domain).");
    }

    for (uint j = 0; j < scratch->size; ++j)
    {
      values[j] = dolfin::DOLFIN_REAL_MAX;
    }
    return;
  }

  Cell cell(*mesh_, cells[0]);

  // Change to global numbering
  scratch->cell.update(cell);

  // Get expansion coefficients on cell
  dofmap_->tabulate_dofs(scratch->dofs, scratch->cell);
  X_->get(scratch->coefficients, scratch->local_dimension, scratch->dofs);

  // Compute linear combination
  for (uint j = 0; j < scratch->size; ++j)
  {
    values[j] = 0.0;
  }
  for (uint i = 0; i < element_->space_dimension(); ++i)
  {
    element_->evaluate_basis(i, scratch->values, x, scratch->cell);
    for (uint j = 0; j < scratch->size; ++j)
    {
      values[j] += scratch->coefficients[i] * scratch->values[j];
      // Check that values are not NaN
      dolfin_assert(values[j] == values[j]);
    }
  }
}

//-----------------------------------------------------------------------------
uint Function::rank() const
{
  dolfin_assert(element_);
  return element_->value_rank();
}

//-----------------------------------------------------------------------------
uint Function::dim(uint i) const
{
  dolfin_assert(element_);
  return element_->value_dimension(i);
}

//-----------------------------------------------------------------------------
uint Function::value_size() const
{
  dolfin_assert(scratch);
  return scratch->size;
}

//-----------------------------------------------------------------------------
void Function::interpolate_vertex_values(real* values) const
{
  // Local data for interpolation on each cell
  uint const num_verts = mesh_->size(0);

  // Make sure vector's ghost values are updated)
  X_->apply();

  // Interpolate vertex values on each cell and pick the last value
  // if two or more cells disagree on the vertex values
  //FIXME: Well... discontinuous approximations might disagree
  if (this->space().is_cellwise_defined())
  {
    error("Interpolation to vertex values is implemented incorrectly for"
          "discontinuous approximations");
  }
  else
  {
    uint const num_cell_vertices = mesh_->type().num_entities(0);
    real* vertex_values = new real[scratch->size * num_cell_vertices];
    for (CellIterator cell(*mesh_); !cell.end(); ++cell)
    {
      // Update to current cell
      scratch->cell.update(*cell);

      // Tabulate dofs
      dofmap_->tabulate_dofs(scratch->dofs, scratch->cell);

      // Pick values from global vector
      X_->get(scratch->coefficients, scratch->local_dimension, scratch->dofs);

      // Interpolate values at the vertices
      // Values are packed by vertex and not by subspace (if any)
      element_->interpolate_vertex_values(vertex_values, scratch->coefficients,
                                         scratch->cell);

      // Copy values to array of vertex values
      for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
      {
        for (uint i = 0; i < scratch->size; ++i)
        {
          values[i * num_verts + vertex->index()] = vertex_values[vertex.pos()
              * scratch->size + i];
        }
      }
    }
    // Delete local data
    delete[] vertex_values;
  }

}

//-----------------------------------------------------------------------------
void Function::interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell) const
{
  // Check dimension
  dolfin_assert(finite_element.space_dimension() == scratch->local_dimension);

  // Tabulate dofs
  dofmap_->tabulate_dofs(scratch->dofs, cell, dolfin_cell);

  // Pick values from global vector if cache mapping is not empty
#ifdef ENABLE_FUNCTION_CACHE
  if (!cache_mapping_->empty())
  {
    for (uint i = 0; i < scratch->local_dimension; ++i)
    {
      _map<uint, uint>::const_iterator it = cache_mapping_->find(scratch->dofs[i]);
      coefficients[i] = data_cache_[it->second];
    }
  }
  else
#endif
  X_->get(coefficients, scratch->local_dimension, scratch->dofs);
}

//-----------------------------------------------------------------------------
void Function::interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell, uint facet) const
{
  interpolate(coefficients, cell, finite_element, dolfin_cell);
}

//-----------------------------------------------------------------------------
GenericVector& Function::vector() const
{
  dolfin_assert(X_);
  return *X_;
}

//-----------------------------------------------------------------------------
FiniteElementSpace const& Function::space() const
{
  dolfin_assert(discrete_space_);
  return *discrete_space_;
}

//-----------------------------------------------------------------------------
void Function::interpolate(GenericFunction const& other)
{

  FunctionInterpolation I(other, *this);
  I.compute();
}

//-----------------------------------------------------------------------------
Array<Function *> Function::decompose()
{
  return FunctionDecomposition::compute(*this);
}

//-----------------------------------------------------------------------------
uint Function::num_sub_functions() const
{
  dolfin_assert(element_);
  return element_->num_sub_elements();
}

//-----------------------------------------------------------------------------
uidx Function::block_size() const
{
  dolfin_assert(dofmap_);
  return dofmap_->dofsmapping_size();
}

//-----------------------------------------------------------------------------
real * Function::create_block() const
{
  dolfin_assert(dofmap_);
  return new real[dofmap_->dofsmapping_size()];
}

//-----------------------------------------------------------------------------
void Function::get_block(real *& values) const
{
  dolfin_assert(X_);
  dolfin_assert(dofmap_);
  if (!values)
  {
    values = new real[dofmap_->dofsmapping_size()];
  }
  X_->apply();
  X_->get(values, dofmap_->dofsmapping_size(), dofmap_->dofsmapping());
}

//-----------------------------------------------------------------------------
void Function::set_block(real *& values)
{
  dolfin_assert(X_);
  dolfin_assert(dofmap_);
  X_->set(values, dofmap_->dofsmapping_size(), dofmap_->dofsmapping());
  sync_ghosts();
}

//-----------------------------------------------------------------------------
void Function::add_block(real *& values)
{
  dolfin_assert(X_);
  dolfin_assert(dofmap_);
  X_->add(values, dofmap_->dofsmapping_size(), dofmap_->dofsmapping());
  sync_ghosts();
}

//-----------------------------------------------------------------------------
void Function::InitializeVector()
{
  if (X_->size() != dofmap_->global_dimension())
  {
    // Specific case in serial local_size == global_dimension
    X_->init(dofmap_->local_size());
  }

  InitializeGhosts();

  X_->zero();
  X_->apply();

  renumbered_ = false;
}

//-----------------------------------------------------------------------------
void Function::InitializeGhosts()
{
  if(!mesh_->is_distributed()) return;

  std::set<uint> indices;

  for (CellIterator cell(*mesh_); !cell.end(); ++cell)
  {
    // Update to current cell
    scratch->cell.update(*cell);

    // Tabulate dofs
    dofmap_->tabulate_dofs(scratch->dofs, scratch->cell);

    for (uint j = 0; j < element_->space_dimension(); ++j)
    {
      indices.insert(scratch->dofs[j]);
    }

  }
  std::map<uint, uint> map = dofmap_->get_map();
  dolfin_assert(map.size() == 0);

  X_->init_ghosted(indices.size(), indices, map);

#ifdef ENABLE_FUNCTION_CACHE
  delete[] indices_;
  delete[] data_cache_;
  delete cache_mapping_;
  cache_mapping_ = new _map<uint, uint>;

  indices_ = new uint[indices.size()];
  data_cache_ = new real[indices.size()];

  uint i = 0;
  std::set<uint>::iterator it;
  for (it = indices.begin(); it != indices.end(); it++)
  {
    indices_[i] = *it;
    (*cache_mapping_)[*it] = i++;
  }

  cache_size_ = indices.size();
#endif
}

//-----------------------------------------------------------------------------
void Function::disp() const
{
  cout << "Function" << endl;
  cout << "--------" << endl;

  // Begin indentation
  begin("");
  if(this->empty())
  {
    message("Empty");
  }
  else
  {
    this->space().disp();
    message("min : %g", this->min());
    message("max : %g", this->max());
  }
  // End indentation
  end();
  skip();
}

//-----------------------------------------------------------------------------
void Function::sync_ghosts()
{

  if(!mesh_->is_distributed()) return;

  if (dofmap_->renumbered() && !renumbered_)
  {
    InitializeGhosts();
    renumbered_ = true;
  }

  X_->apply();

#ifdef ENABLE_FUNCTION_CACHE
  if (indices_)
  {
    X_->get(data_cache_, cache_size_, indices_);
  }
#endif
}

//-----------------------------------------------------------------------------
Function& Function::operator=(Function const& other)
{
  if(this == &other)
  {
    return *this;
  }

  if(this->empty())
  {
    discrete_space_ = new FiniteElementSpace(*other.discrete_space_);
    element_ = &discrete_space_->element();
    dofmap_ = &discrete_space_->dofmap();
    scratch = new ScratchSpace(*discrete_space_);
    X_ = new Vector();
    renumbered_ = false;
    cache_size_ = 0;
    indices_ = NULL;
    data_cache_ = NULL;
    //
    InitializeVector();
  }
  else if(this->space() != other.space())
  {
    error("Function : attempt to assign function with different space");
  }

  // Copy vector
  *X_ = *other.X_;

  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator+=(Function const& other)
{
  dolfin_assert(!this->empty());
  dolfin_assert(!other.empty());
  dolfin_assert(this->space() == other.space());
  this->vector() += other.vector();
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator-=(Function const& other)
{
  dolfin_assert(!this->empty());
  dolfin_assert(!other.empty());
  dolfin_assert(this->space() == other.space());
  this->vector() -= other.vector();
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator*=(Function const& other)
{
  dolfin_assert(!this->empty());
  dolfin_assert(!other.empty());
  dolfin_assert(this->space() == other.space());
  this->vector() *= other.vector();
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::axpy(real value, Function const& other)
{
  dolfin_assert(!this->empty());
  dolfin_assert(!other.empty());
  dolfin_assert(this->space() == other.space());
  this->vector().axpy(value, other.vector());
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::swap(Function& other)
{
  std::swap(const_cast<Mesh *&>(this->mesh_), const_cast<Mesh *&>(other.mesh_));
  std::swap(this->discrete_space_, other.discrete_space_);
  std::swap(this->element_, other.element_);
  std::swap(this->dofmap_, other.dofmap_);
  std::swap(this->scratch, other.scratch);
  std::swap(this->X_, other.X_);
  std::swap(this->renumbered_, other.renumbered_);
#ifdef ENABLE_FUNCTION_CACHE
  std::swap(this->cache_size_, other.cache_size_);
  std::swap(this->indices_, other.indices_);
  std::swap(this->data_cache_, other.data_cache_);
  std::swap(this->cache_mapping_, other.cache_mapping_);
#endif
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator=(real value)
{
  dolfin_assert(!this->empty());
  this->vector() = value;
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator+=(real value)
{
  dolfin_assert(!this->empty());
  error("Not implemented");
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator-=(real value)
{
  dolfin_assert(!this->empty());
  error("Not implemented");
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator*=(real value)
{
  dolfin_assert(!this->empty());
  this->vector() *= value;
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator/=(real value)
{
  dolfin_assert(!this->empty());
  this->vector() /= value;
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::zero()
{
  dolfin_assert(!this->empty());
  this->vector().zero();
  return *this;
}

//-----------------------------------------------------------------------------
real Function::min() const
{
  dolfin_assert(!this->empty());
  return this->vector().min();
}

//-----------------------------------------------------------------------------
real Function::max() const
{
  dolfin_assert(!this->empty());
  return this->vector().max();
}

//-----------------------------------------------------------------------------

} /* *namespace dolfin */

