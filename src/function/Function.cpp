// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/function/Function.h>

#include <dolfin/common/AdjacentMapping.h>
#include <dolfin/common/maybe_unused.h>
#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/SubFunction.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------

Function::Function()
  : GenericFunction()
  , TimeDependent()
{
  // Do nothing
}

//-----------------------------------------------------------------------------

Function::Function( Mesh & mesh )
  : GenericFunction()
  , TimeDependent()
  , mesh_( &mesh )
  , fe_space_( nullptr )
  , scratch_( nullptr )
  , X_( nullptr )
  , renumbered_( false )
{
  // Do nothing
}

//-----------------------------------------------------------------------------

Function::Function( FiniteElementSpace const & space )
  : GenericFunction()
  , TimeDependent()
  , mesh_( &space.mesh() )
  , fe_space_( new FiniteElementSpace( space ) )
  , scratch_( new ScratchSpace( *fe_space_ ) )
  , X_( new Vector() )
  , renumbered_( false )
{
  // Initialise function
  InitializeVector();
}

//-----------------------------------------------------------------------------

Function::Function( SubFunction const & sub_function )
  : GenericFunction()
  , TimeDependent( sub_function.function() )
  , mesh_( &sub_function.function().mesh() )
  , fe_space_( new FiniteElementSpace( sub_function.function().space(),
                                             sub_function.index() ) )
  , scratch_( new ScratchSpace( *fe_space_ ) )
  , X_( new Vector() )
  , renumbered_( false )
{
  // Initialize vector, scratch space and ghosts
  InitializeVector();

  // Copy subvector, naive implementation
  Function& gFunc = sub_function.function();
  DofMap const& gDm = gFunc.space().dofmap();
  size_t const gLocalDim = gDm.num_element_dofs;
  size_t const gDmOffset = gDm.sub_dofmaps_offsets()[sub_function.index()];
  size_t const thisLocalDim = scratch_->local_dimension;

  // Sync ghosts before getting the block
  real * gblock = gFunc.create_block();
  gFunc.sync();
  gFunc.get_block(gblock);

  // Loop baby, loop...
  size_t gBlockOffset = 0;
  size_t ii = 0;
  real * this_block = this->create_block();
  for (CellIterator cell(*mesh_); !cell.end(); ++cell, gBlockOffset += gLocalDim)
  {
    for (size_t dof = 0; dof < thisLocalDim; ++dof)
    {
      this_block[ii++] = gblock[gBlockOffset + gDmOffset + dof];
    }
  }
  this->set_block(this_block);

  delete[] this_block;
  delete[] gblock;
}

//-----------------------------------------------------------------------------

Function::Function( Function const & other )
  : GenericFunction()
  , TimeDependent( other )
  , mesh_( &other.mesh() )
  , fe_space_( nullptr )
  , scratch_( nullptr )
  , X_( nullptr )
  , renumbered_( false )
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

void Function::init(Form& form, size_t i)
{
  if(mesh_ == nullptr)
  {
    const_cast<Mesh *&>(mesh_) = &form.dofmaps()[i]->mesh();
  }

  if(mesh_ != &form.dofmaps()[i]->mesh())
  {
    error("Function : mesh mismatch between function and coefficient %d", i);
  }

  //
  clear();
  fe_space_ = new FiniteElementSpace( form.spaces()[i] );
  scratch_ = new ScratchSpace(*fe_space_);
  X_ = new Vector();
  //
  InitializeVector();
}

//-----------------------------------------------------------------------------

void Function::init(FiniteElementSpace const& space)
{
  if(mesh_ == nullptr)
  {
    const_cast<Mesh *&>(mesh_) = &space.mesh();
  }
  if(mesh_ != &space.mesh())
  {
    error("Function : mesh mismatch between function and space");
  }
  //
  clear();
  fe_space_ = new FiniteElementSpace(space);
  scratch_ = new ScratchSpace(*fe_space_);
  X_ = new Vector();
  //
  InitializeVector();
}

//-----------------------------------------------------------------------------

void Function::clear()
{
  if ( X_ != nullptr )
  {
    delete X_;
    X_ = nullptr;
  }

  if ( fe_space_ != nullptr )
  {
    delete fe_space_;
    fe_space_ = nullptr;
  }

  if ( scratch_ != nullptr )
  {
    delete scratch_;
    scratch_ = nullptr;
  }

  renumbered_ = true;
}

//-----------------------------------------------------------------------------

void Function::evaluate(size_t n, real* values, const real* x,
                        const ufc::cell& cell) const
{
  dolfin_assert( values != nullptr );
  dolfin_assert( x != nullptr );

  UFCCell const * ufc_cell = static_cast<UFCCell const *>(&cell);

  // Get expansion coefficients on cell
  fe_space_->dofmap().tabulate_dofs(scratch_->dofs.data(), *ufc_cell);
  X_->get(scratch_->coefficients.data(), scratch_->local_dimension, scratch_->dofs.data());
  std::fill_n(values, n * scratch_->size, 0.0);
  for (size_t q = 0; q < n; ++q, values+=scratch_->size, x+=cell.geometric_dimension)
  {
// FIXME reintroduce this, once all other things are worked out
// #ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
//     dolfin_assert( scratch_->size <= Space::MAX_DIMENSION );
//     dolfin_assert( fe_space_->element().space_dimension() <= scratch_->space_dimension );

//     real coord[Space::MAX_DIMENSION];

//     fe_space_->element().evaluate_basis_map_coordinates( coord[0], coord[1], coord[2],
//                                               x, *ufc_cell );

//     fe_space_->element().evaluate_basis_from_coordinates( coord[0], coord[1], coord[2],
//                                                scratch_->all_basis_values );

//     // Compute linear combination
//     for (size_t i = 0; i < fe_space_->element().space_dimension(); ++i)
//       for (size_t j = 0; j < scratch_->size; ++j)
//         values[j] += scratch_->coefficients[i] * scratch_->all_basis_values[i][j];
// #else
    // Compute linear combination
    for (size_t i = 0; i < fe_space_->element().space_dim; ++i)
    {
      fe_space_->element().ufc().evaluate_basis( i, scratch_->basis_values.data(),
                                                 x, ufc_cell->coordinates.data(),
                                                 ufc_cell->orientation );
      for (size_t j = 0; j < scratch_->size; ++j)
      {
        values[j] += scratch_->coefficients[i] * scratch_->basis_values[j];
      }
    }
// #endif
  }
}

//-----------------------------------------------------------------------------

void Function::eval(real* values, const real* x) const
{
  // Find the cell that contains x
  Point p( x );
  std::vector<size_t> cells;
  mesh_->intersector().overlap(p, cells);
  if (cells.size() == 0)
  {
    if (!mesh_->is_distributed())
    {
      warning("Unable to evaluate function at given point (not inside domain).");
    }

    for (size_t j = 0; j < scratch_->size; ++j)
    {
      values[j] = std::numeric_limits<real>::infinity();
    }
    return;
  }
  Cell cell(*mesh_, cells[0]);
  scratch_->cell.update(cell);
  evaluate(values, x, scratch_->cell);
}

//-----------------------------------------------------------------------------

void Function::interpolate_vertex_values(real* values) const
{
  // Local data for interpolation on each cell
  size_t const num_verts = mesh_->size(0);

  // Make sure vector's ghost values are updated)
  X_->apply();

  // Interpolate vertex values on each cell and pick the last value
  // if two or more cells disagree on the vertex values
  //FIXME: Well... discontinuous approximations might disagree
  if (mesh_->num_global_cells() > 1 && this->space().is_cellwise_defined())
  {
    //FIXME: imported from naive 2013ish code before mappings existed.
    size_t const num_cell_vertices = mesh_->type().num_entities(0);
    std::fill(values, values + scratch_->size * mesh_->num_vertices(), 0.0);
    real* vertex_sumwghts = new real[mesh_->num_vertices()];
    std::fill(vertex_sumwghts, vertex_sumwghts + mesh_->num_vertices(), 0.0);
    real* vertex_values = new real[scratch_->size * num_cell_vertices];
    for (CellIterator cell(*mesh_); !cell.end(); ++cell)
    {
      // Update to current cell
      scratch_->cell.update(*cell);

      // Tabulate dofs
      fe_space_->dofmap().tabulate_dofs(scratch_->dofs.data(), scratch_->cell);

      // Pick values from global vector
      X_->get(scratch_->coefficients.data(), scratch_->local_dimension, scratch_->dofs.data());

      // Interpolate values at the vertices
      // Values are packed by vertex and not by subspace (if any)
      // FIXME find a form, where more tha vertex_values and dof_values are used,
      // to make sure the arguments here are correct
      fe_space_->element().ufc().interpolate_vertex_values(vertex_values, scratch_->coefficients.data(),
                                                           scratch_->cell.coordinates.data(), 0);

      // Sum values to array of vertex values
      for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
      {
        static real const w = 1.0;  // plan for other weights
        vertex_sumwghts[vertex->index()] += w;
        for (size_t i = 0; i < scratch_->size; ++i)
        {
          values[i * num_verts + vertex->index()] +=
              w * vertex_values[vertex.pos() * scratch_->size + i];
        }
      }
    }

    //
    delete[] vertex_values;

    if (mesh_->is_distributed())
    {
#ifdef DOLFIN_HAVE_MPI
      size_t const rank            = dolfin::MPI::rank();
      size_t const pe_size         = dolfin::MPI::size();
      DistributedData const& dist0 = mesh_->distdata()[0];

      std::vector< std::vector<real> > sendbuf( pe_size );

      // Send sum of local weights
      for (SharedIterator it(dist0); it.valid(); ++it)
      {
        for ( size_t const & adj : it.adj() )
        {
          sendbuf[adj].push_back(vertex_sumwghts[it.index()]);
          for (size_t i = 0; i < scratch_->size; ++i)
          {
            sendbuf[adj].push_back(values[i * num_verts + it.index()]);
          }
        }
      }

      // Exchange data

      MPI_Status status;
      size_t src;
      size_t dst;

      //FIXME: Overallocate
      std::vector< size_t > recvbuf( dist0.num_shared() );
      int recvcount = 0;
      for (size_t j = 1; j < pe_size; ++j)
      {
        src = (rank - j + pe_size) % pe_size;
        dst = (rank + j) % pe_size;

        MPI::check_error( MPI_Sendrecv(sendbuf[dst].data(), sendbuf[dst].size(),
                                       MPI_type< real >::value, dst, 1,
                                       recvbuf.data(), recvbuf.size(),
                                       MPI_type< size_t >::value, src, 1,
                                       MPI::DOLFIN_COMM, &status) );
        MPI::check_error( MPI_Get_count(&status, MPI_type< size_t >::value, &recvcount) );

        // Add contributions, just simplified this part with mappings
        std::vector<size_t> recvmapping = dist0.shared_mapping().from(src);
        dolfin_assert(recvmapping.size() == (size_t) recvcount);
        for (int k = 0; k < recvcount; k += (1 + scratch_->size))
        {
          dolfin_assert(dist0.has_local(recvmapping[k]) > 0);
          vertex_sumwghts[recvmapping[k]] +=  recvbuf[k];
          for (size_t i = 1; i <= scratch_->size; ++i)
          {
            values[i * num_verts + recvmapping[k]] += recvbuf[k + i];
          }
        }
      }
#endif
    }


    // Average
    for (size_t vindex = 0; vindex < num_verts; ++vindex)
    {
      for (size_t i = 0; i < scratch_->size; ++i)
      {
        values[i * num_verts + vindex] /= vertex_sumwghts[vindex];
      }
    }

    // Delete local data
    delete[] vertex_sumwghts;

  }
  else
  {
    size_t const num_cell_vertices = mesh_->type().num_entities(0);
    real* vertex_values = new real[scratch_->size * num_cell_vertices];
    for (CellIterator cell(*mesh_); !cell.end(); ++cell)
    {
      // Update to current cell
      scratch_->cell.update(*cell);

      // Tabulate dofs
      fe_space_->dofmap().tabulate_dofs(scratch_->dofs.data(), scratch_->cell);

      // Pick values from global vector
      X_->get(scratch_->coefficients.data(), scratch_->local_dimension, scratch_->dofs.data());

      // Interpolate values at the vertices
      // Values are packed by vertex and not by subspace (if any)
      // FIXME find a form, where more tha vertex_values and dof_values are used,
      // to make sure the arguments here are correct
      fe_space_->element().ufc().interpolate_vertex_values(vertex_values, scratch_->coefficients.data(),
                                                           scratch_->cell.coordinates.data(), 0);

      // Copy values to array of vertex values
      for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
      {
        for (size_t i = 0; i < scratch_->size; ++i)
        {
          values[i * num_verts + vertex->index()] = vertex_values[vertex.pos()
              * scratch_->size + i];
        }
      }
    }
    // Delete local data
    delete[] vertex_values;
  }

}

//-----------------------------------------------------------------------------

void Function::interpolate( real *                      coefficients,
                            UFCCell const &             cell,
                            ufc::finite_element const & finite_element) const
{
  // Check dimension
  dolfin_assert(finite_element.space_dimension() == scratch_->local_dimension);
  MAYBE_UNUSED(finite_element);

  // Tabulate dofs
  fe_space_->dofmap().tabulate_dofs(scratch_->dofs.data(), cell );

  // Pick values from global vector if cache mapping is not empty
#ifdef ENABLE_FUNCTION_CACHE
  if ( not cache_mapping_.empty() )
  {
    for (size_t i = 0; i < scratch_->local_dimension; ++i)
    {
      _map<size_t, size_t>::const_iterator it = cache_mapping_.find(scratch_->dofs[i]);
      coefficients[i] = data_cache_[it->second];
    }
    return;
  }
#endif

  X_->get(coefficients, scratch_->local_dimension, scratch_->dofs.data());
}

//-----------------------------------------------------------------------------

void Function::InitializeVector()
{
  dolfin_assert( fe_space_->element().tdim == mesh_->topology().dim() );
  if ( X_->size() != fe_space_->dofmap().global_dim )
  {
    // Specific case in serial local_size == global_dimension
    X_->init(fe_space_->dofmap().local_size());
  }

  InitializeGhosts();

  X_->zero();
  X_->apply();

  renumbered_ = false;
}

//-----------------------------------------------------------------------------

void Function::InitializeGhosts()
{
  if ( not mesh_->is_distributed() )
  {
    return;
  }

  _ordered_set<size_t> indices;

  for (CellIterator cell(*mesh_); !cell.end(); ++cell)
  {
    // Update to current cell
    scratch_->cell.update(*cell);

    // Tabulate dofs
    fe_space_->dofmap().tabulate_dofs(scratch_->dofs.data(), scratch_->cell);

    for (size_t j = 0; j < fe_space_->element().space_dim; ++j)
    {
      indices.insert(scratch_->dofs[j]);
    }
  }

  ghost_mapping_.clear();
  X_->init_ghosted(indices.size(), indices, ghost_mapping_);

#ifdef ENABLE_FUNCTION_CACHE
  cache_mapping_.clear();
  indices_.resize( indices.size() );
  data_cache_.resize( indices.size() );

  size_t i = 0;
  for ( size_t const & index : indices )
  {
    indices_[i]           = index;
    cache_mapping_[index] = i++;
  }
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

  if ( not renumbered_)
  {
    InitializeGhosts();
    renumbered_ = true;
  }

  X_->apply();

#ifdef ENABLE_FUNCTION_CACHE
  if ( not indices_.empty() )
  {
    X_->get( data_cache_.data(), indices_.size(), indices_.data() );
  }
#endif
}

//-----------------------------------------------------------------------------

auto Function::operator=(Function const& other) -> Function&
{
  if(this == &other)
  {
    return *this;
  }

  if(this->empty())
  {
    const_cast<Mesh *&>(mesh_) = other.mesh_;
    fe_space_ = new FiniteElementSpace(*other.fe_space_);
    scratch_ = new ScratchSpace(*fe_space_);
    X_ = new Vector();
    renumbered_ = false;
    indices_.clear();
    data_cache_.clear();
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

auto Function::swap(Function& other) -> Function&
{
  TimeDependent::swap(other);
  std::swap(const_cast<Mesh *&>(this->mesh_), const_cast<Mesh *&>(other.mesh_));
  std::swap(this->fe_space_, other.fe_space_);
  std::swap(this->scratch_, other.scratch_);
  std::swap(this->X_, other.X_);
  std::swap(this->renumbered_, other.renumbered_);
#ifdef ENABLE_FUNCTION_CACHE
  std::swap(this->indices_, other.indices_);
  std::swap(this->data_cache_, other.data_cache_);
  std::swap(this->cache_mapping_, other.cache_mapping_);
#endif
  return *this;
}

//-----------------------------------------------------------------------------

} /* *namespace dolfin */

