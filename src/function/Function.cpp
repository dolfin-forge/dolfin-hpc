// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/function/Function.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/common/AdjacentMapping.h>
#include <dolfin/common/maybe_unused.h>
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

namespace dolfin
{

//-----------------------------------------------------------------------------
Function::Function() :
    GenericFunction(),
    TimeDependent(),
    mesh_(NULL),
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
Function::Function(Mesh& mesh) :
    GenericFunction(),
    TimeDependent(),
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
    TimeDependent(),
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
    TimeDependent(),
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
    TimeDependent(),
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
    TimeDependent(sub_function.function()),
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
  gFunc.sync();
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
    TimeDependent(other),
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
  if(mesh_ == NULL)
  {
    const_cast<Mesh *&>(mesh_) = &form.dofmaps()[i].mesh();
  }
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
  if(mesh_ == NULL)
  {
    const_cast<Mesh *&>(mesh_) = &space.mesh();
  }
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
  evaluate( 1, values, x, cell );
}

//--- GenericFunction ---------------------------------------------------------
Mesh& Function::mesh() const
{
  return (*mesh_);
}

//-----------------------------------------------------------------------------
void Function::evaluate(uint n, real* values, const real* x,
                        const ufc::cell& cell) const
{
  dolfin_assert( values != NULL );
  dolfin_assert( x != NULL );

  UFCCell const * ufc_cell = static_cast<UFCCell const *>(&cell);

  // Get expansion coefficients on cell
  dofmap_->tabulate_dofs(scratch->dofs, *ufc_cell);
  X_->get(scratch->coefficients, scratch->local_dimension, scratch->dofs);
  std::fill_n(values, n * scratch->size, 0.0);
  for (uint q = 0; q < n; ++q, values+=scratch->size, x+=cell.geometric_dimension)
  {
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
    dolfin_assert( scratch->size <= Space::MAX_DIMENSION );
    dolfin_assert( element_->space_dimension() <= scratch->space_dimension );

    real coord[Space::MAX_DIMENSION];

    element_->evaluate_basis_map_coordinates( coord[0], coord[1], coord[2],
                                              x, *ufc_cell );

    element_->evaluate_basis_from_coordinates( coord[0], coord[1], coord[2],
                                               scratch->all_basis_values );

    // Compute linear combination
    for (uint i = 0; i < element_->space_dimension(); ++i)
      for (uint j = 0; j < scratch->size; ++j)
        values[j] += scratch->coefficients[i] * scratch->all_basis_values[i][j];
#else
    // Compute linear combination
    for (uint i = 0; i < element_->space_dimension(); ++i)
    {
      //FIXME: Idiotic
      element_->evaluate_basis(i, scratch->basis_values, x, *ufc_cell);
      for (uint j = 0; j < scratch->size; ++j)
      {
        values[j] += scratch->coefficients[i] * scratch->basis_values[j];
      }
    }
#endif
  }
}

//-----------------------------------------------------------------------------
void Function::eval(real* values, const real* x) const
{
  // Find the cell that contains x
  Point p(mesh_->geometry_dimension(), x);
  Array<uint> cells;
  mesh_->intersector().overlap(p, cells);
  if (cells.size() == 0)
  {
    if (!mesh_->is_distributed())
    {
      error("Unable to evaluate function at given point (not inside domain).");
    }

    for (uint j = 0; j < scratch->size; ++j)
    {
      values[j] = std::numeric_limits<real>::infinity();
    }
    return;
  }
  Cell cell(*mesh_, cells[0]);
  scratch->cell.update(cell);
  evaluate(values, x, scratch->cell);
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
  if (mesh_->num_global_cells() > 1 && this->space().is_cellwise_defined())
  {
    //FIXME: imported from naive 2013ish code before mappings existed.
    uint const num_cell_vertices = mesh_->type().num_entities(0);
    std::fill(values, values + scratch->size * mesh_->num_vertices(), 0.0);
    real* vertex_sumwghts = new real[mesh_->num_vertices()];
    std::fill(vertex_sumwghts, vertex_sumwghts + mesh_->num_vertices(), 0.0);
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

      // Sum values to array of vertex values
      for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
      {
        static real const w = 1.0;  // plan for other weights
        vertex_sumwghts[vertex->index()] += w;
        for (uint i = 0; i < scratch->size; ++i)
        {
          values[i * num_verts + vertex->index()] +=
              w * vertex_values[vertex.pos() * scratch->size + i];
        }
      }
    }

    //
    delete[] vertex_values;

    if (mesh_->is_distributed())
    {
#ifdef DOLFIN_HAVE_MPI
      uint const rank = dolfin::MPI::rank();
      uint const pe_size = dolfin::MPI::size();
      DistributedData const& dist0 = mesh_->distdata()[0];
      Array<real> * sendbuf = new Array<real> [pe_size];
      // Send sum of local weights
      for (SharedIterator it(dist0); it.valid(); ++it)
      {
        _set<uint> const& adjs = it.adj();
        for (_set<uint>::const_iterator a = adjs.begin(); a != adjs.end(); ++a)
        {
          sendbuf[*a].push_back(vertex_sumwghts[it.index()]);
          for (uint i = 0; i < scratch->size; ++i)
          {
            sendbuf[*a].push_back(values[i * num_verts + it.index()]);
          }
        }
      }

      // Exchange data

      MPI_Status status;
      uint src;
      uint dst;

      //FIXME: Overallocate
      uint recvsize = dist0.num_shared();
      uint * recvbuf = (recvsize == 0 ? NULL : new uint[recvsize]);
      int recvcount;
      for (uint j = 1; j < pe_size; ++j)
      {
        src = (rank - j + pe_size) % pe_size;
        dst = (rank + j) % pe_size;

        MPI::check_error( MPI_Sendrecv(&sendbuf[dst][0], sendbuf[dst].size(),
                                       MPI_UNSIGNED, dst, 1,  &recvbuf[0],
                                       recvsize, MPI_DOUBLE, src, 1,
                                       MPI::DOLFIN_COMM, &status) );
        MPI::check_error( MPI_Get_count(&status, MPI_DOUBLE, &recvcount) );

        // Add contributions, just simplified this part with mappings
        Array<uint> recvmapping = dist0.shared_mapping().from(src);
        dolfin_assert(recvmapping.size() == (uint) recvcount);
        for (int k = 0; k < recvcount; k += (1 + scratch->size))
        {
          dolfin_assert(dist0.has_local(recvmapping[k]) > 0);
          vertex_sumwghts[recvmapping[k]] +=  recvbuf[k];
          for (uint i = 1; i <= scratch->size; ++i)
          {
            values[i * num_verts + recvmapping[k]] += recvbuf[k + i];
          }
        }
      }

      //
      delete[] recvbuf;
      delete[] sendbuf;
#endif
    }


    // Average
    for (uint vindex = 0; vindex < num_verts; ++vindex)
    {
      for (uint i = 0; i < scratch->size; ++i)
      {
        values[i * num_verts + vindex] /= vertex_sumwghts[vindex];
      }
    }

    // Delete local data
    delete[] vertex_sumwghts;

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
  MAYBE_UNUSED(finite_element);

  // Tabulate dofs
  dofmap_->tabulate_dofs(scratch->dofs, cell, dolfin_cell);

  // Pick values from global vector if cache mapping is not empty
#ifdef ENABLE_FUNCTION_CACHE
  if (cache_mapping_ != NULL)
  {
    for (uint i = 0; i < scratch->local_dimension; ++i)
    {
      _map<uint, uint>::const_iterator it = cache_mapping_->find(scratch->dofs[i]);
      coefficients[i] = data_cache_[it->second];
    }
    return;
  }
#endif

  X_->get(coefficients, scratch->local_dimension, scratch->dofs);
}

//-----------------------------------------------------------------------------
void Function::interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell, uint) const
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
void Function::operator<<(Expression const& other)
{
  FunctionInterpolation::compute(other, *this);
}

//-----------------------------------------------------------------------------
void Function::operator<<(Coefficient const& other)
{
  FunctionInterpolation::compute(other, *this);
}

//-----------------------------------------------------------------------------
void Function::operator<<(GenericFunction const& other)
{
  FunctionInterpolation::compute(other, *this);
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
  sync();
}

//-----------------------------------------------------------------------------
void Function::add_block(real *& values)
{
  dolfin_assert(X_);
  dolfin_assert(dofmap_);
  X_->add(values, dofmap_->dofsmapping_size(), dofmap_->dofsmapping());
  sync();
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

  _ordered_set<uint> indices;

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
  _ordered_map<uint, uint> map = dofmap_->get_map();
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
  _ordered_set<uint>::iterator it;
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
  cout << "Function\n";
  cout << "--------\n";

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
void Function::sync()
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
    const_cast<Mesh *&>(mesh_) = other.mesh_;
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
  this->sync();

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
  TimeDependent::swap(other);
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
Function& Function::operator+=(real)
{
  dolfin_assert(!this->empty());
  error("Not implemented");
  return *this;
}

//-----------------------------------------------------------------------------
Function& Function::operator-=(real)
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

