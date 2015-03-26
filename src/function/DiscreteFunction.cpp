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

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
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
#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(GenericVector& x, Form& form, uint i) :
    GenericFunction(),
    mesh_(form.mesh(i)),
    discrete_space_(form, i),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(false),
    X_(&x),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, GenericVector& x, Form& form,
                                   uint i) :
    GenericFunction(),
    mesh_(mesh),
    discrete_space_(mesh, form, i),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(false),
    X_(&x),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Form& form, uint i) :
    GenericFunction(),
    mesh_(form.mesh(i)),
    discrete_space_(form, i),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(true),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh, Form& form, uint i) :
    GenericFunction(),
    mesh_(mesh),
    discrete_space_(mesh, form, i),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(true),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(GenericVector& x,
                                   FiniteElementSpace const& space) :
    GenericFunction(),
    mesh_(space.mesh()),
    discrete_space_(space),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(false),
    X_(&x),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(FiniteElementSpace const& space) :
    GenericFunction(),
    mesh_(space.mesh()),
    discrete_space_(space),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(true),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(Mesh& mesh,
                                   ufl::FiniteElementBase const& finite_element) :
    GenericFunction(),
    mesh_(mesh),
    discrete_space_(mesh, finite_element),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(true),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------
DiscreteFunction::DiscreteFunction(SubFunction const& sub_function) :
    GenericFunction(),
    mesh_(sub_function.function().mesh()),
    discrete_space_(sub_function.function().space(), sub_function.index()),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(true),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  // Initialize vector, scratch space and ghosts
  InitializeVector();

  // Copy subvector, naive implementation
  DiscreteFunction& gFunc = sub_function.function();
  DofMap const& gDm = gFunc.space().dofmap();
  uint const gLocalDim = gDm.local_dimension();
  uint const gDmOffset = gDm.sub_dofmaps_offsets()[sub_function.index()];
  uint const thisLocalDim = scratch.local_dimension;

  // Sync ghosts before getting the block
  real * gblock = gFunc.create_block();
  gFunc.sync_ghosts();
  gFunc.get_block(gblock);

  // Loop baby, loop...
  uint gBlockOffset = 0;
  uint ii = 0;
  real * this_block = this->create_block();
  for (CellIterator cell(mesh_); !cell.end(); ++cell, gBlockOffset += gLocalDim)
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
DiscreteFunction::DiscreteFunction(DiscreteFunction const& f) :
    GenericFunction(),
    mesh_(f.mesh()),
    discrete_space_(f.space()),
    element_(discrete_space_.element()),
    dofmap_(discrete_space_.dofmap()),
    scratch(discrete_space_),
    local_vector_(true),
    X_(new Vector()),
    renumbered_(false),
    cache_size_(0),
    indices_(NULL),
    data_cache_(NULL)
{
  dolfin_assert(this->space() == f.space());

  // Copy vector
  *X_ = *f.X_;

  this->sync_ghosts();
}

//-----------------------------------------------------------------------------
DiscreteFunction::~DiscreteFunction()
{
  if (local_vector_)
  {
    delete X_;
  }
  delete[] indices_;
  delete[] data_cache_;
}

//-----------------------------------------------------------------------------
DiscreteFunction const& DiscreteFunction::operator=(DiscreteFunction const& f)
{
  // Check that data matches
  if ((this->space() != f.space()) || (X_->size() != f.X_->size()))
  {
    error("Assignment of discrete function failed."
          "Finite element spaces or dimensions don't match.");
  }

  // Copy vector
  *X_ = *f.X_;

  this->sync_ghosts();

  return *this;
}

//--- UFC INTERFACE -----------------------------------------------------------
//-----------------------------------------------------------------------------
void DiscreteFunction::evaluate(real* values, const real* x,
                                const ufc::cell& cell) const
{
  //FIXME: Inconsistency of interface
  UFCCell const * ufc_cell = dynamic_cast<UFCCell const *>(&cell);

  // Get expansion coefficients on cell
  dofmap_.tabulate_dofs(scratch.dofs, *ufc_cell);
  X_->get(scratch.coefficients, scratch.local_dimension, scratch.dofs);

  // Compute linear combination
  for (uint j = 0; j < scratch.size; j++)
  {
    values[j] = 0.0;
  }
  FiniteElement const& fe = discrete_space_.element();
  for (uint i = 0; i < fe.space_dimension(); i++)
  {
    fe.evaluate_basis(i, scratch.values, x, *ufc_cell);
    for (uint j = 0; j < scratch.size; j++)
    {
      values[j] += scratch.coefficients[i] * scratch.values[j];
    }
  }
}

//--- GenericFunction ---------------------------------------------------------
Mesh& DiscreteFunction::mesh() const
{
  return mesh_;
}

//-----------------------------------------------------------------------------
uint DiscreteFunction::rank() const
{
  return element_.value_rank();
}

//-----------------------------------------------------------------------------
uint DiscreteFunction::dim(uint i) const
{
  return element_.value_dimension(i);
}

//-----------------------------------------------------------------------------
void DiscreteFunction::interpolate_vertex_values(real* values) const
{
  // Local data for interpolation on each cell
  uint const tdim = mesh_.topology().dim();
  uint const num_verts = mesh_.numVertices();

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
    uint const num_cell_vertices = mesh_.type().numVertices(tdim);
    real* vertex_values = new real[scratch.size * num_cell_vertices];
    MeshDistributedData& distdata = mesh_.distdata();
    for (CellIterator cell(mesh_); !cell.end(); ++cell)
    {
      // Update to current cell
      scratch.cell.update(*cell, distdata);

      // Tabulate dofs
      dofmap_.tabulate_dofs(scratch.dofs, scratch.cell, cell->index());

      // Pick values from global vector
      X_->get(scratch.coefficients, scratch.local_dimension, scratch.dofs);

      // Interpolate values at the vertices
      // Values are packed by vertex and not by subspace (if any)
      element_.interpolate_vertex_values(vertex_values, scratch.coefficients,
                                         scratch.cell);

      // Copy values to array of vertex values
      for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
      {
        for (uint i = 0; i < scratch.size; ++i)
        {
          values[i * num_verts + vertex->index()] = vertex_values[vertex.pos()
              * scratch.size + i];
        }
      }
    }
    // Delete local data
    delete[] vertex_values;
  }

}

//-----------------------------------------------------------------------------
void DiscreteFunction::interpolate(real* coefficients, const ufc::cell& cell,
                                   const ufc::finite_element& finite_element,
                                   const Cell& dolfin_cell) const
{
  // Check dimension
  dolfin_assert(finite_element.space_dimension() == scratch.local_dimension);

  // Tabulate dofs
  dofmap_.tabulate_dofs(scratch.dofs, cell, dolfin_cell);

  // Pick values from global vector if cache mapping is not empty
#ifdef ENABLE_FUNCTION_CACHE
  if (!cache_mapping_.empty())
  {
    for (uint i = 0; i < scratch.local_dimension; i++)
    {
      _map<uint, uint>::const_iterator it = cache_mapping_.find(scratch.dofs[i]);
      coefficients[i] = data_cache_[it->second];
    }
  }
  else
#endif
  X_->get(coefficients, scratch.local_dimension, scratch.dofs);
}

//-----------------------------------------------------------------------------
void DiscreteFunction::eval(real* values, const real* x) const
{
  // Find the cell that contains x
  uint const gdim = mesh_.geometry().dim();
  if (gdim > 3)
  {
    error("Sorry, point evaluation of functions not implemented for meshes of "
          "dimension %d.",
          gdim);
  }
  Point p;
  for (uint i = 0; i < gdim; i++)
  {
    p[i] = x[i];
  }
  Array<uint> cells;
  mesh_.intersector().overlap(p, cells);
  if (cells.size() < 1)
  {
    if (!mesh_.is_distributed())
    {
      error("Unable to evaluate function at given point (not inside domain).");
    }

    for (uint j = 0; j < scratch.size; j++)
    {
      values[j] = dolfin::DOLFIN_REAL_MAX;
    }
    return;
  }

  Cell cell(mesh_, cells[0]);

  // Change to global numbering
  scratch.cell.update(cell, mesh_.distdata());

  // Get expansion coefficients on cell
  dofmap_.tabulate_dofs(scratch.dofs, scratch.cell);
  X_->get(scratch.coefficients, scratch.local_dimension, scratch.dofs);

  // Compute linear combination
  for (uint j = 0; j < scratch.size; j++)
  {
    values[j] = 0.0;
  }
  for (uint i = 0; i < element_.space_dimension(); i++)
  {
    element_.evaluate_basis(i, scratch.values, x, scratch.cell);
    for (uint j = 0; j < scratch.size; j++)
    {
      values[j] += scratch.coefficients[i] * scratch.values[j];
      // Check that values are not NaN
      dolfin_assert(values[j] == values[j]);
    }
  }
}

//-----------------------------------------------------------------------------
GenericVector& DiscreteFunction::vector() const
{
  dolfin_assert(X_);
  return *X_;
}

//-----------------------------------------------------------------------------
FiniteElementSpace const& DiscreteFunction::space() const
{
  return discrete_space_;
}

//-----------------------------------------------------------------------------
std::string const DiscreteFunction::signature() const
{
  return element_.signature();
}

//-----------------------------------------------------------------------------
uint const DiscreteFunction::num_sub_functions() const
{
  return element_.num_sub_elements();
}

//-----------------------------------------------------------------------------
void DiscreteFunction::interpolate(Function const& other_func)
{
  if ((this->rank() != other_func.rank())
      || (this->value_size() != other_func.value_size())
      || (this->dim(0) != other_func.dim(0)))
  {
    error("Attempt to interpolate between functions of different value shape");
  }

  // Make sure vectors ghost values are updated)
  X_->apply();

  // Pretabulated version
  real * block = this->create_block();

  // Can the space be flattened into scalar finite elements ?
  Array<ufc::finite_element const*> const& flt_elms =
      this->space().element().flatten();
  bool const flattenable = (flt_elms.size() == this->value_size());

  // The other function being discrete we need to interpolate.
  if ((other_func.type() == Function::discrete) || !flattenable)
  {
    ScratchSpace other_scratch(other_func.space());
    if (this->mesh() == other_func.mesh())
    {
      uint dof = 0;
      uint const local_dim = dofmap_.local_dimension();
      for (CellIterator cell(mesh_); !cell.end(); ++cell, dof += local_dim)
      {
        scratch.cell.update(*cell, mesh_.distdata());
        this->space().element().evaluate_dofs(&block[dof], other_func,
                                              scratch.cell);
      }
    }
    else
    {
      error("Interpolation on non-matching meshes is not implemented.");
    }
  }
  // Analytical expression and flattened space (naive implementation)
  else
  {
    uint dof = 0;
    for (CellIterator cell(mesh_); !cell.end(); ++cell)
    {
      scratch.cell.update(*cell, mesh_.distdata());
      dofmap_.tabulate_coordinates(scratch.coordinates, scratch.cell);

      uint celldof = 0;
      for (uint lfspace = 0; lfspace < flt_elms.size(); ++lfspace)
      {
        for (uint ii = 0; ii < flt_elms[lfspace]->space_dimension(); ++ii)
        {
          other_func.eval(scratch.values, scratch.coordinates[celldof++]);
          block[dof++] = scratch.values[lfspace];
        }
      }
      //
      dolfin_assert(celldof == scratch.local_dimension);
    }
  }

  this->set_block(block);
  delete[] block;
}

//-----------------------------------------------------------------------------
real * DiscreteFunction::create_block() const
{
  return new real[dofmap_.dofsmapping_size()];
}

//-----------------------------------------------------------------------------
void DiscreteFunction::get_block(real *& values) const
{
  if (!values)
  {
    values = new real[dofmap_.dofsmapping_size()];
  }
  X_->apply();
  X_->get(values, dofmap_.dofsmapping_size(), dofmap_.dofsmapping());
}

//-----------------------------------------------------------------------------
void DiscreteFunction::set_block(real *& values)
{
  X_->set(values, dofmap_.dofsmapping_size(), dofmap_.dofsmapping());
  sync_ghosts();
}

//-----------------------------------------------------------------------------
void DiscreteFunction::add_block(real *& values)
{
  if (!values)
  {
    values = new real[dofmap_.dofsmapping_size()];
  }
  X_->add(values, dofmap_.dofsmapping_size(), dofmap_.dofsmapping());
}

//-----------------------------------------------------------------------------
void DiscreteFunction::InitializeVector()
{
  if (X_->size() != dofmap_.global_dimension())
  {
    // Specific case in serial local_size == global_dimension
    X_->init(dofmap_.local_size());
  }

  InitializeGhosts();

  X_->zero();
  X_->apply();

  renumbered_ = false;
}

//-----------------------------------------------------------------------------
void DiscreteFunction::InitializeGhosts()
{
  if(!mesh_.is_distributed()) return;

  std::set<uint> indices;

  MeshDistributedData& distdata = mesh_.distdata();
  for (CellIterator cell(mesh_); !cell.end(); ++cell)
  {
    // Update to current cell
    scratch.cell.update(*cell, distdata);

    // Tabulate dofs
    dofmap_.tabulate_dofs(scratch.dofs, scratch.cell, cell->index());

    for (uint j = 0; j < element_.space_dimension(); ++j)
    {
      indices.insert(scratch.dofs[j]);
    }

  }
  std::map<uint, uint> map = dofmap_.getMap();
  dolfin_assert(map.size() == 0);

  X_->init_ghosted(indices.size(), indices, map);

#ifdef ENABLE_FUNCTION_CACHE
  delete[] indices_;
  delete[] data_cache_;

  cache_mapping_.clear();

  indices_ = new uint[indices.size()];
  data_cache_ = new real[indices.size()];

  uint i = 0;
  std::set<uint>::iterator it;
  for (it = indices.begin(); it != indices.end(); it++)
  {
    indices_[i] = *it;
    cache_mapping_[*it] = i++;
  }

  cache_size_ = indices.size();
#endif
}

//-----------------------------------------------------------------------------
void DiscreteFunction::disp() const
{
  cout << "DiscreteFunction" << endl;
  cout << "----------------" << endl;

  // Begin indentation
  begin("");
  GenericFunction::disp();
  this->space().disp();
  // End indentation
  end();
  skip();
}

//-----------------------------------------------------------------------------
void DiscreteFunction::sync_ghosts()
{

  if(!mesh_.is_distributed()) return;

  if (dofmap_.renumbered() && !renumbered_)
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

}

