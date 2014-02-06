// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Dag Lindbo, 2008.
// Modified by Kristen Kaasbjerg, 2008.
// Modified by Niclas Jansson, 2008-2010.
// Modified by Aurélien Larcher 2013.
//
// First added:  2007-04-02
// Last changed: 2013-09-13

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
    discrete_space_(mesh, form, i),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(false),
    X_(&x),
    intersection_detector_(NULL),
    renumbered(false),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, Form& form, uint i) :
    GenericFunction(mesh),
    discrete_space_(mesh, form, i),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(true),
    X_(new Vector()),
    intersection_detector_(NULL),
    renumbered(false),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, GenericVector& x,
                                   std::string const& finite_element_signature,
                                   std::string const& dof_map_signature) :
    GenericFunction(mesh),
    discrete_space_(mesh, finite_element_signature, dof_map_signature),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(false),
    X_(&x),
    intersection_detector_(NULL),
    renumbered(false),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh,
                                   std::string const& finite_element_signature,
                                   std::string const& dof_map_signature) :
    GenericFunction(mesh),
    discrete_space_(mesh, finite_element_signature, dof_map_signature),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(true),
    X_(new Vector()),
    intersection_detector_(NULL),
    renumbered(false),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, GenericVector& x,
                                   std::string const& finite_element_signature) :
    GenericFunction(mesh),
    discrete_space_(mesh, finite_element_signature),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(false),
    X_(&x),
    intersection_detector_(NULL),
    renumbered(false),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh,
                                   std::string const& finite_element_signature) :
    GenericFunction(mesh),
    discrete_space_(mesh, finite_element_signature),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(true),
    X_(new Vector()),
    intersection_detector_(NULL),
    renumbered(false),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{
  // Initialise function
  InitializeVector();
}

#if ENABLE_UFL
//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh,
                                   ufl::FiniteElementBase const& finite_element) :
    GenericFunction(mesh),
    discrete_space_(mesh, finite_element),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(true),
    X_(new Vector()),
    intersection_detector_(NULL),
    renumbered(false),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{
  // Initialise function
  InitializeVector();
}
#endif

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(SubFunction& sub_function) :
    GenericFunction(sub_function.function().mesh),
    discrete_space_(sub_function.function().space(), sub_function.index()),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(true),
    X_(new Vector()),
    intersection_detector_(NULL),
    renumbered(false),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{
  // Initialize vector, scratch space and ghosts
  InitializeVector();

  // Copy subvector, naive implementation
  DiscreteFunction& global_func = sub_function.function();
  DofMap const& global_dm = global_func.dofmap();
  FiniteElement const& global_fe = global_func.finite_element();
  uint const global_dm_offset =
      global_dm.sub_dof_maps_offsets()[sub_function.index()];

  uint * global_dofs = new uint[global_dm.local_dimension()];
  uint global_block_size = global_func.finite_element().space_dimension();
  real * global_block = new real[global_block_size];

  CellIterator cell(mesh);
  UFCCell ufccell(*cell);
  for (; !cell.end(); ++cell)
  {
    ufccell.update(*cell, mesh.distdata());
    global_dm.tabulate_dofs(global_dofs, ufccell, cell->index());
    dof_map_.tabulate_dofs(scratch.dofs, ufccell, cell->index());

    global_func.vector().get(global_block, global_block_size, global_dofs);
    for (uint dof_id = 0; dof_id < local_dimension_; ++dof_id)
    {
      scratch.coefficients[dof_id] = global_block[global_dm_offset + dof_id];
    }
    this->vector().set(scratch.coefficients, local_dimension_, scratch.dofs);
    this->vector().apply();
  }
  sync_ghosts();

  delete[] global_block;
  delete[] global_dofs;
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(const DiscreteFunction& f) :
    GenericFunction(f.mesh),
    discrete_space_(f.mesh, f.finite_element_.signature()),
    finite_element_(discrete_space_.element()),
    dof_map_(discrete_space_.dofmap()),
    scratch(discrete_space_.scratch),
    local_dimension_(0),
    local_vector_(true),
    X_(new Vector()),
    intersection_detector_(NULL),
    _cache_size(0),
    _indices(NULL),
    _data_cache(NULL)
{

  // Copy vector
  *X_ = *f.X_;

  renumbered = f.renumbered;
}

//-----------------------------------------------------------------------------
DiscreteFunction::~DiscreteFunction()
{
  delete X_;
  delete intersection_detector_;
  delete[] _indices;
  delete[] _data_cache;
}

//-----------------------------------------------------------------------------
dolfin::uint DiscreteFunction::rank() const
{
  return finite_element_.value_rank();
}

//-----------------------------------------------------------------------------
dolfin::uint DiscreteFunction::dim(uint i) const
{
  return finite_element_.value_dimension(i);
}

//-----------------------------------------------------------------------------
FiniteElementSpace const& DiscreteFunction::space() const
{
  return discrete_space_;
}

//-----------------------------------------------------------------------------
dolfin::uint DiscreteFunction::numSubFunctions() const
{
  return finite_element_.num_sub_elements();
}

//-----------------------------------------------------------------------------
const DiscreteFunction& DiscreteFunction::operator=(const DiscreteFunction& f)
{
  // Check that data matches
  if (strcmp(finite_element_.signature(), f.finite_element_.signature()) != 0
      || strcmp(dof_map_.signature(), f.dof_map_.signature()) != 0
      || X_->size() != f.X_->size())
  {
    error("Assignment of discrete function failed."
          "Finite element spaces or dimensions don't match.");
  }

  // Copy vector
  *X_ = *f.X_;

  return *this;
}

//-----------------------------------------------------------------------------
void DiscreteFunction::interpolate_vertex_values(real* values) const
{
  // Local data for interpolation on each cell
  CellIterator cell(mesh);
  UFCCell ufc_cell(*cell);
  const uint num_cell_vertices = mesh.type().numVertices(mesh.topology().dim());
  real* vertex_values = new real[scratch.size * num_cell_vertices];

  // Make sure vector's ghost values are updated)
  X_->apply();

  // Interpolate vertex values on each cell and pick the last value
  // if two or more cells disagree on the vertex values
  // Well... discontiuous approximations might disagree
  for (; !cell.end(); ++cell)
  {
    // Update to current cell
    ufc_cell.update(*cell, mesh.distdata());

    // Tabulate dofs
    dof_map_.tabulate_dofs(scratch.dofs, ufc_cell, cell->index());

    // Pick values from global vector
    X_->get(scratch.coefficients, local_dimension_, scratch.dofs);

    // Interpolate values at the vertices
    finite_element_.interpolate_vertex_values(vertex_values,
                                              scratch.coefficients, ufc_cell);

    // Copy values to array of vertex values
    for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
    {
      for (uint i = 0; i < scratch.size; ++i)
      {
        values[i * mesh.numVertices() + vertex->index()] =
            vertex_values[vertex.pos() * scratch.size + i];
      }
    }
  }

  // Delete local data
  delete[] vertex_values;
}

//-----------------------------------------------------------------------------
void DiscreteFunction::interpolate(real* coefficients, const ufc::cell& cell,
                                   const ufc::finite_element& finite_element,
                                   const Cell& dolfin_cell) const
{
  // Check dimension
  dolfin_assert(finite_element.space_dimension() == local_dimension_);

  // Tabulate dofs
  dof_map_.tabulate_dofs(scratch.dofs, cell, dolfin_cell.index());

  // Pick values from global vector
#ifdef ENABLE_FUNCTION_CACHE
  if (MPI::numProcesses() > 1)
  {
    for (uint i = 0; i < local_dimension_; i++)
    {
      _map<uint, uint>::const_iterator it = _cache_mapping.find(scratch.dofs[i]);
      coefficients[i] = _data_cache[it->second];
    }
  }
  else
#endif
  X_->get(coefficients, local_dimension_, scratch.dofs);
}

//-----------------------------------------------------------------------------
void DiscreteFunction::interpolate(Function const& other_func)
{
  Array<uint> const& value_dims = finite_element_.sub_value_dimensions(0);
  Array<uint> const& value_offs = finite_element_.sub_value_offsets(0);
  Array<uint> const& dm_dims = dof_map_.sub_dof_maps_dimensions();
  Array<uint> const& dm_offs = dof_map_.sub_dof_maps_offsets();
  uint const nb_subspaces = dm_dims.size();

  // Make sure vectors ghost values are updated)
  X_->apply();

  // Cell tabulated version
  CellIterator cell(mesh);
  UFCCell ufccell(*cell);
  real * values = new real[finite_element_.value_dimension(0)];
  real * block = new real[dof_map_.dofsmapping_size()];
  uint cell_offset = 0;
  uint const local_dim = dof_map_.local_dimension();
  for (; !cell.end(); ++cell, cell_offset += local_dim)
  {
    ufccell.update(*cell, mesh.distdata());
    dof_map_.tabulate_coordinates(scratch.coordinates, ufccell);

    uint dof_id = 0;
    for (uint sub = 0; sub < nb_subspaces; ++sub)
    {
      uint sub_val_dim = value_dims[sub];
      uint nb_nodes = dm_dims[sub] / sub_val_dim;
      uint off = dm_offs[sub];
      for (uint sub_id = 0; sub_id < nb_nodes; ++sub_id)
      {
        other_func.eval(values, scratch.coordinates[sub_id]);
        for (uint v = 0; v < sub_val_dim; ++v)
        {
          block[cell_offset + off + v * nb_nodes + sub_id] =
              values[value_offs[sub] + v];
        }
        ++dof_id;
      }
    }
  }
  this->set(block);

  delete[] block;
  delete[] values;
}

//-----------------------------------------------------------------------------
void DiscreteFunction::eval(real* values, const real* x) const
{
  // Initialize intersection detector if not done before
  if (!intersection_detector_)
  {
    intersection_detector_ = new IntersectionDetector(mesh);
  }

  // Find the cell that contains x
  uint const gdim = mesh.geometry().dim();
  if (gdim > 3)
  {
    error(
        "Sorry,"
        "point evaluation of functions not implemented for meshes of dimension %d.",
        gdim);
  }
  Point p;
  for (uint i = 0; i < gdim; i++)
  {
    p[i] = x[i];
  }
  Array<uint> cells;
  intersection_detector_->overlap(p, cells);
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
  dof_map_.tabulate_dofs(scratch.dofs, ufc_cell, cell.index());
  X_->get(scratch.coefficients, local_dimension_, scratch.dofs);

  // Compute linear combination
  for (uint j = 0; j < scratch.size; j++)
  {
    values[j] = 0.0;
  }
  for (uint i = 0; i < finite_element_.space_dimension(); i++)
  {
    finite_element_.evaluate_basis(i, scratch.values, x, ufc_cell);
    for (uint j = 0; j < scratch.size; j++)
    {
      values[j] += scratch.coefficients[i] * scratch.values[j];
    }
  }
}

//-----------------------------------------------------------------------------
std::string const DiscreteFunction::signature() const
{
  return finite_element_.signature();
}

//-----------------------------------------------------------------------------
GenericVector& DiscreteFunction::vector() const
{
  dolfin_assert(X_);
  return *X_;
}

//-----------------------------------------------------------------------------
DofMap const& DiscreteFunction::dofmap() const
{
  return dof_map_;
}

//-----------------------------------------------------------------------------
FiniteElement const& DiscreteFunction::finite_element() const
{
  return finite_element_;
}

//-----------------------------------------------------------------------------
void DiscreteFunction::get(real *& values)
{
  if (!values)
  {
    values = new real[dof_map_.dofsmapping_size()];
  }
  X_->get(values, dof_map_.dofsmapping_size(), dof_map_.dofsmapping());
}

//-----------------------------------------------------------------------------
void DiscreteFunction::set(real *& values)
{
  X_->set(values, dof_map_.dofsmapping_size(), dof_map_.dofsmapping());
  sync_ghosts();
}

//-----------------------------------------------------------------------------
void DiscreteFunction::InitializeVector()
{
  if (X_->size() != dof_map_.global_dimension())
  {
    if (MPI::numProcesses() > 1)
    {
      X_->init(dof_map_.local_size());
    }
    else
    {
      X_->init(dof_map_.global_dimension());
    }
  }
  renumbered = false;
}

//-----------------------------------------------------------------------------
void DiscreteFunction::InitializeGhosts()
{
  std::set<uint> indices;
  CellIterator cell(mesh);
  UFCCell ufc_cell(*cell);

  for (; !cell.end(); ++cell)
  {
    // Update to current cell
    ufc_cell.update(*cell, mesh.distdata());

    // Tabulate dofs
    dof_map_.tabulate_dofs(scratch.dofs, ufc_cell, cell->index());

    for (uint j = 0; j < finite_element_.space_dimension(); ++j)
    {
      indices.insert(scratch.dofs[j]);
    }

  }
  std::map<uint, uint> map = dof_map_.getMap();

  X_->init_ghosted(indices.size(), indices, map);

#ifdef ENABLE_FUNCTION_CACHE
  if (_indices)
    delete[] _indices;
  if (_data_cache)
    delete[] _data_cache;

  _cache_mapping.clear();

  _indices = new uint[indices.size()];
  _data_cache = new real[indices.size()];

  uint i = 0;
  std::set<uint>::iterator it;
  for (it = indices.begin(); it != indices.end(); it++)
  {
    _indices[i] = *it;
    _cache_mapping[*it] = i++;
  }

  _cache_size = indices.size();
#endif
}

//-----------------------------------------------------------------------------
void DiscreteFunction::sync_ghosts()
{

  if (MPI::numProcesses() == 1)
    return;

  if (dof_map_.renumbered() && !renumbered && MPI::numProcesses() > 1)
  {
    InitializeGhosts();
    renumbered = true;
  }

  X_->apply();

#ifdef ENABLE_FUNCTION_CACHE
  if (_indices)
  {
    X_->get(_data_cache, _cache_size, _indices);
  }
#endif
}

}

