// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Dag Lindbo, 2008.
// Modified by Kristen Kaasbjerg, 2008.
// Modified by Niclas Jansson, 2008-2010.
// Modified by Aurélien Larcher 2013
//
// First added:  2007-04-02
// Last changed: 2010-06-16

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/UFCMesh.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/function/SubFunction.h>
#include <dolfin/function/DiscreteFunction.h>

#include <cstring>
#include <iomanip>
#include <limits>
#include <set>

namespace dolfin
{

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, GenericVector& x, Form& form,
    uint i) :
        GenericFunction(mesh),
        finite_element_(NULL),
        ufc_finite_element_(NULL),
        dof_map_(NULL),
        local_dim_(0),
        dofs_coordinates_(NULL),
        cell_dof_values_(NULL),
        local_vector(false),
        x(&x),
        intersection_detector(NULL),
        scratch(NULL),
        renumbered(false),
        _cache_size(0),
        _indices(NULL),
        data_cache(NULL)
{
  // Initialise function
  __init(mesh, form, i);
}
//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, Form& form, uint i) :
        GenericFunction(mesh),
        finite_element_(NULL),
        ufc_finite_element_(NULL),
        dof_map_(NULL),
        local_dim_(0),
        dofs_coordinates_(NULL),
        cell_dof_values_(NULL),
        local_vector(true),
        x(new Vector()),
        intersection_detector(NULL),
        scratch(NULL),
        renumbered(false),
        _cache_size(0),
        _indices(NULL),
        data_cache(NULL)
{
  // Initialise function
  __init(mesh, form, i);
}
//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, GenericVector& x,
    std::string finite_element_signature, std::string dof_map_signature) :
        GenericFunction(mesh),
        finite_element_(NULL),
        ufc_finite_element_(NULL),
        dof_map_(NULL),
        local_dim_(0),
        dofs_coordinates_(NULL),
        cell_dof_values_(NULL),
        local_vector(false),
        x(&x),
        intersection_detector(NULL),
        scratch(NULL),
        renumbered(false),
        _cache_size(0),
        _indices(NULL),
        data_cache(NULL)
{
  __init(mesh, finite_element_signature, dof_map_signature);
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh,
    std::string finite_element_signature, std::string dof_map_signature) :
        GenericFunction(mesh),

        finite_element_(NULL),
        ufc_finite_element_(NULL),
        dof_map_(NULL),
        local_dim_(0),
        dofs_coordinates_(NULL),
        cell_dof_values_(NULL),
        local_vector(true),
        x(new Vector()),
        intersection_detector(NULL),
        scratch(NULL),
        renumbered(false),
        _cache_size(0),
        _indices(NULL),
        data_cache(NULL)
{
  __init(mesh, finite_element_signature, dof_map_signature);
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, GenericVector& x,
    std::string finite_element_signature) :
        GenericFunction(mesh),
        finite_element_(NULL),
        ufc_finite_element_(NULL),
        dof_map_(NULL),
        local_dim_(0),
        dofs_coordinates_(NULL),
        cell_dof_values_(NULL),
        local_vector(false),
        x(&x),
        intersection_detector(NULL),
        scratch(NULL),
        renumbered(false),
        _cache_size(0),
        _indices(NULL),
        data_cache(NULL)
{
  __init(mesh, finite_element_signature,
      DofMap::dofmap_signature(finite_element_signature));
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh,
    std::string finite_element_signature) :
        GenericFunction(mesh),

        finite_element_(NULL),
        ufc_finite_element_(NULL),
        dof_map_(NULL),
        local_dim_(0),
        dofs_coordinates_(NULL),
        cell_dof_values_(NULL),
        local_vector(true),
        x(new Vector()),
        intersection_detector(NULL),
        scratch(NULL),
        renumbered(false),
        _cache_size(0),
        _indices(NULL),
        data_cache(NULL)
{
  __init(mesh, finite_element_signature,
      DofMap::dofmap_signature(finite_element_signature));
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(SubFunction& sub_function) :
        GenericFunction(sub_function.f->mesh),
        finite_element_(NULL),
        ufc_finite_element_(NULL),
        dof_map_(NULL),
        local_dim_(0),
        dofs_coordinates_(NULL),
        cell_dof_values_(NULL),
        local_vector(true),
        x(new Vector()),
        intersection_detector(NULL),
        scratch(NULL),
        renumbered(false),
        _cache_size(0),
        _indices(NULL),
        data_cache(NULL)
{
  // Create sub system
  SubSystem sub_system(sub_function.i);
  DiscreteFunction const& global_func = *sub_function.f;

  // Extract sub element (return value is a newly created finite_element)
  ufc_finite_element_ = global_func.finite_element_->create_sub_element(
      sub_system.array());
  finite_element_ = new FiniteElement(*ufc_finite_element_, true);
  message(
      "SubFunction finite element:"
          + std::string(ufc_finite_element_->signature()));

  // Extract sub dof map and offset
  uint offset = 0;
  ufc::dof_map * ufc_dof_map = global_func.dof_map_->create_sub_dof_map(
      sub_system.array(), offset);

  // Token is requested by the standalone function
  dof_map_ = DofMapCache::instance().acquire_dofmap(mesh,
      ufc_dof_map->signature());

  delete ufc_dof_map;

  // Initialize vector, scratch space and init ghosts
  __init();

  // Copy subvector, naive implementation
  uint const global_dm_offset =
      global_func.dof_map_->sub_dof_maps_offsets()[sub_function.i];
  DofMap const& global_dm = global_func.dofmap();
  uint * global_dofs = new uint[global_dm.local_dimension()];
  uint global_block_size = global_func.finite_element().space_dimension();
  real * global_block = new real[global_block_size];

  CellIterator cell(mesh);
  UFCCell ufccell(*cell);
  for (; !cell.end(); ++cell)
  {
    ufccell.update(*cell, mesh.distdata());
    global_dm.tabulate_dofs(global_dofs, ufccell, cell->index());
    dof_map_->tabulate_dofs(scratch->dofs, ufccell, cell->index());

    global_func.vector().get(global_block, global_block_size, global_dofs);
    for (uint dof_id = 0; dof_id < local_dim_; ++dof_id)
    {
      cell_dof_values_[dof_id] = global_block[global_dm_offset + dof_id];
    }
    this->vector().set(cell_dof_values_, local_dim_, scratch->dofs);
    this->vector().apply();
  }
  sync_ghosts();

  delete[] global_block;
  delete[] global_dofs;
}
//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(const DiscreteFunction& f) :
        GenericFunction(f.mesh),
        finite_element_(NULL),
        ufc_finite_element_(NULL),
        dof_map_(NULL),
        local_dim_(0),
        dofs_coordinates_(NULL),
        cell_dof_values_(NULL),
        local_vector(true),
        x(new Vector()),
        intersection_detector(NULL),
        scratch(NULL),
        _cache_size(0),
        _indices(NULL),
        data_cache(NULL)
{
  __init(f.mesh, f.ufc_finite_element_->signature(), f.dof_map_->signature());

  // Copy vector
  *x = *f.x;

  renumbered = f.renumbered;
}
//-----------------------------------------------------------------------------
DiscreteFunction::~DiscreteFunction()
{
  DofMapCache::instance().release_dofmap(*dof_map_);

  if (dofs_coordinates_)
  {
    for (uint i = 0; i < local_dim_; ++i)
    {
      delete[] dofs_coordinates_[i];
    }
    delete[] dofs_coordinates_;
  }

  if (cell_dof_values_)
    delete[] cell_dof_values_;

  if (finite_element_)
    delete finite_element_;

  if (local_vector)
    delete x;

  if (intersection_detector)
    delete intersection_detector;

  if (scratch)
    delete scratch;

  if (_indices)
    delete[] _indices;

  if (data_cache)
    delete[] data_cache;
}
//-----------------------------------------------------------------------------
dolfin::uint
DiscreteFunction::rank() const
{
  dolfin_assert(ufc_finite_element_);
  return ufc_finite_element_->value_rank();
}
//-----------------------------------------------------------------------------
dolfin::uint
DiscreteFunction::dim(uint i) const
{
  dolfin_assert(ufc_finite_element_);
  return ufc_finite_element_->value_dimension(i);
}
//-----------------------------------------------------------------------------
uint
DiscreteFunction::degree() const
{
  return finite_element_->degree();
}
//-----------------------------------------------------------------------------
dolfin::uint
DiscreteFunction::numSubFunctions() const
{
  dolfin_assert(ufc_finite_element_);
  return ufc_finite_element_->num_sub_elements();
}
//-----------------------------------------------------------------------------
const DiscreteFunction&
DiscreteFunction::operator=(const DiscreteFunction& f)
{
  // Check that data matches
  if (strcmp(ufc_finite_element_->signature(),
      f.ufc_finite_element_->signature()) != 0
      || strcmp(dof_map_->signature(), f.dof_map_->signature()) != 0
      || x->size() != f.x->size())
  {
    error(
        "Assignment of discrete function failed. Finite element spaces or dimensions don't match.");
  }

  // Copy vector
  *x = *f.x;

  return *this;
}
//-----------------------------------------------------------------------------
void
DiscreteFunction::interpolate(real* values) const
{
  // Local data for interpolation on each cell
  CellIterator cell(mesh);
  UFCCell ufc_cell(*cell);
  const uint num_cell_vertices = mesh.type().numVertices(mesh.topology().dim());
  real* vertex_values = new real[scratch->size * num_cell_vertices];

  // Make sure vectors ghost values are updated)
  x->apply();

  // Interpolate vertex values on each cell and pick the last value
  // if two or more cells disagree on the vertex values
  // Well... discontiuous approximations might disagree
  for (; !cell.end(); ++cell)
  {
    // Update to current cell
    ufc_cell.update(*cell, mesh.distdata());

    // Tabulate dofs
    dof_map_->tabulate_dofs(scratch->dofs, ufc_cell, cell->index());

    // Pick values from global vector
    x->get(scratch->coefficients, local_dim_, scratch->dofs);

    // Interpolate values at the vertices
    ufc_finite_element_->interpolate_vertex_values(vertex_values,
        scratch->coefficients, ufc_cell);

    // Copy values to array of vertex values
    for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
    {
      for (uint i = 0; i < scratch->size; ++i)
      {
        values[i * mesh.numVertices() + vertex->index()] =
            vertex_values[vertex.pos() * scratch->size + i];
      }
    }
  }

  // Delete local data
  delete[] vertex_values;
}
//-----------------------------------------------------------------------------
void
DiscreteFunction::interpolate(real* coefficients, const ufc::cell& cell,
    const ufc::finite_element& finite_element, const Cell& dolfin_cell) const
{
  // Check dimension
  if (finite_element.space_dimension() != local_dim_)
    error(
        "Finite element does not match for interpolation of discrete function.");

  // Tabulate dofs
  dof_map_->tabulate_dofs(scratch->dofs, cell, dolfin_cell.index());

  // Pick values from global vector
#ifdef ENABLE_FUNCTION_CACHE
  if (MPI::numProcesses() > 1)
  {
    for (uint i = 0; i < local_dim_; i++)
    {
      _map<uint, uint>::const_iterator it = cache_mapping.find(scratch->dofs[i]);
      coefficients[i] = data_cache[it->second];
    }
  }
  else
#endif
  x->get(coefficients, local_dim_, scratch->dofs);
}
//-----------------------------------------------------------------------------
void
DiscreteFunction::interpolate(Function const& other_func)
{
  real * values = new real[finite_element_->value_dimension(0)];

  Array<uint> const& value_dims = finite_element_->sub_value_dimensions(0);
  Array<uint> const& value_offs = finite_element_->sub_value_offsets(0);
  Array<uint> const& dm_dims = dof_map_->sub_dof_maps_dimensions();
  Array<uint> const& dm_offs = dof_map_->sub_dof_maps_offsets();
  uint const nb_subspaces = dm_dims.size();

  CellIterator cell(mesh);
  UFCCell ufccell(*cell);
  for (; !cell.end(); ++cell)
  {
    ufccell.update(*cell, mesh.distdata());
    dof_map_->tabulate_dofs(scratch->dofs, ufccell, cell->index());
    dof_map_->tabulate_coordinates(dofs_coordinates_, ufccell);

    uint dof_id = 0;
    for (uint sub = 0; sub < nb_subspaces; ++sub)
    {
      uint sub_val_dim = value_dims[sub];
      uint nb_nodes = dm_dims[sub] / sub_val_dim;
      uint off = dm_offs[sub];
      for (uint sub_id = 0; sub_id < nb_nodes; ++sub_id)
      {
        other_func.eval(values, dofs_coordinates_[sub_id]);
        for (uint v = 0; v < sub_val_dim; ++v)
        {
          cell_dof_values_[off + v * nb_nodes + sub_id] = values[value_offs[sub] + v];
        }
        ++dof_id;
      }
    }
    this->vector().set(cell_dof_values_, local_dim_, scratch->dofs);
    this->vector().apply();
  }
  sync_ghosts();

  delete[] values;
}
//-----------------------------------------------------------------------------
void
DiscreteFunction::eval(real* values, const real* x) const
{
  dolfin_assert(scratch);
// Initialize intersection detector if not done before
  if (!intersection_detector)
  {
    intersection_detector = new IntersectionDetector(mesh);
  }

// Find the cell that contains x
  uint const gdim = mesh.geometry().dim();
  if (gdim > 3)
  {
    error(
        "Sorry, point evaluation of functions not implemented for meshes of dimension %d.",
        gdim);
  }
  Point p;
  for (uint i = 0; i < gdim; i++)
  {
    p[i] = x[i];
  }
  Array<uint> cells;
  intersection_detector->overlap(p, cells);
  if (cells.size() < 1)
  {
    if (MPI::numProcesses() == 1)
      warning(
          "Unable to evaluate function at given point (not inside domain).");

    values[0] = 1e50;
    values[1] = 1e50;
    values[2] = 1e50;
    return;
  }

  Cell cell(mesh, cells[0]);
  UFCCell ufc_cell(cell);

// Change to global numbering
  ufc_cell.update(cell, mesh.distdata());

// Get expansion coefficients on cell
  this->interpolate(scratch->coefficients, ufc_cell, *ufc_finite_element_,
      cell);

// Compute linear combination
  for (uint j = 0; j < scratch->size; j++)
  {
    values[j] = 0.0;
  }
  for (uint i = 0; i < ufc_finite_element_->space_dimension(); i++)
  {
    ufc_finite_element_->evaluate_basis(i, scratch->values, x, ufc_cell);
    for (uint j = 0; j < scratch->size; j++)
    {
      values[j] += scratch->coefficients[i] * scratch->values[j];
    }
  }
}
//-----------------------------------------------------------------------------
std::string
DiscreteFunction::signature() const
{
  if (!ufc_finite_element_)
    error("No finite element has been associated with this DiscreteFunction.");

  return ufc_finite_element_->signature();
}
//-----------------------------------------------------------------------------
GenericVector&
DiscreteFunction::vector() const
{
  if (!x)
    error("Vector associated with DiscreteFunction has not been initialised.");

  return *x;
}
//-----------------------------------------------------------------------------
DofMap const&
DiscreteFunction::dofmap() const
{
  dolfin_assert(dof_map_);
  return *dof_map_;
}
//-----------------------------------------------------------------------------
FiniteElement const&
DiscreteFunction::finite_element() const
{
  dolfin_assert(finite_element_);
  return *finite_element_;
}
//-----------------------------------------------------------------------------
void
DiscreteFunction::get(real *& values)
{
  if (!values)
  {
    values = new real[dof_map_->dofsmapping_size()];
  }
  x->get(values, dof_map_->dofsmapping_size(), dof_map_->dofsmapping());
}
//-----------------------------------------------------------------------------
void
DiscreteFunction::set(real *& values)
{
  x->set(values, dof_map_->dofsmapping_size(), dof_map_->dofsmapping());
  sync_ghosts();
}
//-----------------------------------------------------------------------------
void
DiscreteFunction::__init(Mesh& mesh, Form& form, uint i)
{
// Create finite element
  finite_element_ = new FiniteElement(mesh, form, i);
  ufc_finite_element_ = finite_element_->ufc_finite_element_;

// Get DofMap pointer from DofMapSet
// Token has been requested by the Form
  dof_map_ = DofMapCache::instance().acquire_dofmap(mesh, form.form(), i);

//Necessary to call acquire for the moment as we release in the destructor
  if (dof_map_ != &form.dofMaps()[i])
  {
    error("Different DofMap object pointed by Form and DiscreteFunction");
  }

  __init();

}

//-----------------------------------------------------------------------------
void
DiscreteFunction::__init(Mesh& mesh, std::string finite_element_signature,
    std::string dof_map_signature)
{
// Create finite element
  finite_element_ = new FiniteElement(finite_element_signature);
  ufc_finite_element_ = finite_element_->ufc_finite_element_;

// Token is requested by the standalone function
  dof_map_ = DofMapCache::instance().acquire_dofmap(mesh, dof_map_signature);

  __init();
}

//-----------------------------------------------------------------------------
void
DiscreteFunction::__init()
{
  local_dim_ = dof_map_->local_dimension();
  if (x->size() != dof_map_->global_dimension())
  {
    if (MPI::numProcesses() > 1)
    {
      x->init(dof_map_->local_size());
    }
    else
    {
      x->init(dof_map_->global_dimension());
    }
  }

  dofs_coordinates_ = new real*[local_dim_];
  for (uint i = 0; i < local_dim_; ++i)
  {
    dofs_coordinates_[i] = new real[3]; // Internally Point is implemented for d = 3
  }
  cell_dof_values_ = new real[local_dim_];

// Initialize scratch space
  if (!scratch)
  {
    scratch = new Scratch(*ufc_finite_element_);
  }
  else
  {
    error("Scratch was not created");
  }

  if (MPI::numProcesses() > 1)
  {
    __init_ghosts();
  }

  renumbered = false;
}

//-----------------------------------------------------------------------------
void
DiscreteFunction::__init_ghosts()
{
  std::set<uint> indices;
  CellIterator cell(mesh);
  UFCCell ufc_cell(*cell);

  for (; !cell.end(); ++cell)
  {
    // Update to current cell
    ufc_cell.update(*cell, mesh.distdata());

    // Tabulate dofs
    dof_map_->tabulate_dofs(scratch->dofs, ufc_cell, cell->index());

    for (uint j = 0; j < ufc_finite_element_->space_dimension(); j++)
    {
      indices.insert(scratch->dofs[j]);
    }

  }
  std::map<uint, uint> map = dof_map_->getMap();

  x->init_ghosted(indices.size(), indices, map);

#ifdef ENABLE_FUNCTION_CACHE
  if (_indices)
    delete[] _indices;
  if (data_cache)
    delete[] data_cache;

  cache_mapping.clear();

  _indices = new uint[indices.size()];
  data_cache = new real[indices.size()];

  uint i = 0;
  std::set<uint>::iterator it;
  for (it = indices.begin(); it != indices.end(); it++)
  {
    _indices[i] = *it;
    cache_mapping[*it] = i++;
  }

  _cache_size = indices.size();
#endif
}
//-----------------------------------------------------------------------------
void
DiscreteFunction::sync_ghosts()
{

  if (MPI::numProcesses() == 1)
    return;

  if (dof_map_->renumbered() && !renumbered && MPI::numProcesses() > 1)
  {
    __init_ghosts();
    renumbered = true;
  }

  x->apply();
#ifdef ENABLE_FUNCTION_CACHE
  if (_indices)
    x->get(data_cache, _cache_size, _indices);
#endif
}
//-----------------------------------------------------------------------------
DiscreteFunction::Scratch::Scratch(ufc::finite_element& finite_element) :
    size(0), dofs(NULL), coefficients(NULL), values(NULL)
{
// Compute size of value (number of entries in tensor value)
  size = 1;
  for (uint i = 0; i < finite_element.value_rank(); i++)
  {
    size *= finite_element.value_dimension(i);
  }

// Initialize local array for mapping of dofs
  dofs = new uint[finite_element.space_dimension()];
  for (uint i = 0; i < finite_element.space_dimension(); i++)
  {
    dofs[i] = 0;
  }

// Initialize local array for expansion coefficients
  coefficients = new real[finite_element.space_dimension()];
  for (uint i = 0; i < finite_element.space_dimension(); i++)
  {
    coefficients[i] = 0.0;
  }

// Initialize local array for values
  values = new real[size];
  for (uint i = 0; i < size; i++)
  {
    values[i] = 0.0;
  }
}
//-----------------------------------------------------------------------------
DiscreteFunction::Scratch::~Scratch()
{
  if (dofs)
    delete[] dofs;

  if (coefficients)
    delete[] coefficients;

  if (values)
    delete[] values;
}
//-----------------------------------------------------------------------------

}

