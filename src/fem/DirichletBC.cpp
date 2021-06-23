// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/DirichletBC.h>

#include <dolfin/common/constants.h>
#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>
#include <dolfin/parameter/parameters.h>

#include <string>

namespace dolfin
{

//-----------------------------------------------------------------------------
DirichletBC::DirichletBC( Coefficient &     g,
                          Mesh &            mesh,
                          SubDomain const & sub_domain,
                          BCMethod          method )
  : BoundaryCondition( "Dirichlet", mesh, sub_domain )
  , conditions( { { g, sub_domain } } )
  , method_( method )
{
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(
  std::vector< std::pair< Coefficient &, SubDomain const & > > conds,
  Mesh &                                                       mesh,
  BCMethod                                                     method )
  : BoundaryCondition( "Dirichlet", mesh, conds.front().second )
  , conditions( conds )
  , method_( method )
{
  dolfin_assert( not conditions.empty() );
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC( Coefficient &     g,
                          Mesh &            mesh,
                          SubDomain const & sub_domain,
                          SubSystem const & sub_system,
                          BCMethod          method )
  : BoundaryCondition( "Dirichlet", mesh, sub_domain, sub_system )
  , conditions( { { g, sub_domain } } )
  , method_( method )
{
}
//-----------------------------------------------------------------------------
void DirichletBC::apply( GenericMatrix &      A,
                         GenericVector &      b,
                         BilinearForm const & form )
{
  apply_impl( A, b, nullptr, form );
}
//-----------------------------------------------------------------------------
void DirichletBC::apply( GenericMatrix &       A,
                         GenericVector &       b,
                         GenericVector const & x,
                         BilinearForm const &  form )
{
  apply_impl( A, b, &x, form );
}
//-----------------------------------------------------------------------------
void DirichletBC::apply_impl( GenericMatrix &       A,
                              GenericVector &       b,
                              GenericVector const * x,
                              BilinearForm const &  form )
{
  if ( form.trial_space() != form.test_space() )
  {
    error(
      "DirichletBC is implemented only for identical test and trial space" );
  }

  if ( entities_.empty() || this->invalid_mesh() )
  {
    if ( ( method_ == topological ) || ( method_ == geometric ) )
    {
      entities_.clear();

      // Build set of boundary facets
      size_t const tdim = mesh().topology_dimension();
      for ( FacetIterator f( mesh() ); !f.end(); ++f )
      {
        bool const on_boundary =
          ( f->num_entities( tdim ) == 1 ) && !f->is_shared();

        // loop over all conditions in reverse order
        for ( size_t c = conditions.size(); c > 0; --c )
        {
          if ( conditions[c - 1].second.enclosed( *f, on_boundary ) )
          {
            // Get cell to which facet belongs (there may be two, but pick
            // first)
            Cell cell( mesh(), f->entities( tdim )[0] );
            entities_.push_back( cell.index() );
            entities_.push_back( cell.index( *f ) );
            entities_.push_back( c - 1 );
            break;
          }
        }
      }
    }

    this->update_mesh_dependency();
  }

  // Simple check
  form.check( A, b );

  // Check compatibility of function g and the test (sub)space
  FiniteElementSpace const & space = form.trial_space();
  ufc::finite_element *      fe =
    space.element().create_sub_element( this->sub_system() );

  for ( size_t c = 0; c < conditions.size(); ++c )
  {
    Coefficient & g_ = conditions[c].first;

    if ( ( fe->value_rank() != g_.rank() )
         || ( fe->value_dimension( 0 ) != g_.dim( 0 ) ) )
    {
      error(
        "Rank and/or value dimension mismatch between function and space.\n"
        "Function : rank = %d, dim = %d; Space : rank = %d, dim = %d.",
        g_.rank(),
        g_.dim( 0 ),
        fe->value_rank(),
        fe->value_dimension( 0 ) );
    }
  }

  delete fe;

  // A map to hold the mapping from boundary dofs to boundary values
  _map< size_t, real > boundary_values;

  // Compute dofs and values
  switch ( method_ )
  {
    case topological:
      computeBCTopological( boundary_values, space, this->sub_system() );
      break;
    case geometric:
      computeBCGeometric( boundary_values, space, this->sub_system() );
      break;
    case pointwise:
      computeBCPointwise( boundary_values, space, this->sub_system() );
      break;
    default:
      error( "Unknown method for application of boundary conditions." );
      break;
  }

  // Copy boundary value data to arrays
  std::vector< size_t > dofs( boundary_values.size() );
  std::vector< real >   values( boundary_values.size() );
  size_t                i = 0;
  for ( std::pair< size_t const, real > & bnd_value : boundary_values )
  {
    dofs[i]     = bnd_value.first;
    values[i++] = bnd_value.second;
  }

  // Modify boundary values for nonlinear problems
  if ( x != nullptr )
  {
    std::vector< real > x_values( boundary_values.size() );
    x->get( x_values.data(), boundary_values.size(), dofs.data() );
    for ( size_t i = 0; i < boundary_values.size(); ++i )
    {
      values[i] -= x_values[i];
    }
  }

  message( "Applying boundary conditions to linear system" );

  // Modify RHS vector (b[i] = value)
  b.set( values.data(), boundary_values.size(), dofs.data() );

  // Modify linear system (A_ii = 1)
  if ( not dolfin_get< bool >( "Krylov keep PC" ) )
  {
    A.ident( boundary_values.size(), dofs.data() );
  }

  // Finalise changes to A
  if ( not dolfin_get< bool >( "Krylov keep PC" ) )
    A.apply();

  // Finalise changes to b
  b.apply();
}
//-----------------------------------------------------------------------------
void DirichletBC::zero( GenericMatrix &       A,
                              BilinearForm const &  form )
{

  if ( form.trial_space() != form.test_space() )
  {
    error(
      "DirichletBC is implemented only for identical test and trial space" );
  }

  if ( entities_.empty() || this->invalid_mesh() )
  {
    if ( ( method_ == topological ) || ( method_ == geometric ) )
    {
      entities_.clear();

      // Build set of boundary facets
      size_t const tdim = mesh().topology_dimension();
      for ( FacetIterator f( mesh() ); !f.end(); ++f )
      {
        bool const on_boundary =
          ( f->num_entities( tdim ) == 1 ) && !f->is_shared();

        // loop over all conditions in reverse order
        for ( size_t c = conditions.size(); c > 0; --c )
        {
          if ( conditions[c - 1].second.enclosed( *f, on_boundary ) )
          {
            // Get cell to which facet belongs (there may be two, but pick
            // first)
            Cell cell( mesh(), f->entities( tdim )[0] );
            entities_.push_back( cell.index() );
            entities_.push_back( cell.index( *f ) );
            entities_.push_back( c - 1 );
            break;
          }
        }
      }
    }

    this->update_mesh_dependency();
  }

  // Check compatibility of function g and the test (sub)space
  FiniteElementSpace const & space = form.trial_space();
  ufc::finite_element *      fe =
    space.element().create_sub_element( this->sub_system() );

  for ( size_t c = 0; c < conditions.size(); ++c )
  {
    Coefficient & g_ = conditions[c].first;

    if ( ( fe->value_rank() != g_.rank() )
         || ( fe->value_dimension( 0 ) != g_.dim( 0 ) ) )
    {
      error(
        "Rank and/or value dimension mismatch between function and space.\n"
        "Function : rank = %d, dim = %d; Space : rank = %d, dim = %d.",
        g_.rank(),
        g_.dim( 0 ),
        fe->value_rank(),
        fe->value_dimension( 0 ) );
    }
  }

  delete fe;

  // A map to hold the mapping from boundary dofs to boundary values
  _map< size_t, real > boundary_values;

  // Compute dofs and values
  switch ( method_ )
  {
    case topological:
      computeBCTopological( boundary_values, space, this->sub_system() );
      break;
    case geometric:
      computeBCGeometric( boundary_values, space, this->sub_system() );
      break;
    case pointwise:
      computeBCPointwise( boundary_values, space, this->sub_system() );
      break;
    default:
      error( "Unknown method for application of boundary conditions." );
      break;
  }

  // Copy boundary value data to arrays
  std::vector< size_t > dofs( boundary_values.size() );
  std::vector< real >   values( boundary_values.size() );
  size_t                i = 0;
  for ( std::pair< size_t const, real > & bnd_value : boundary_values )
  {
    dofs[i]     = bnd_value.first;
    values[i++] = bnd_value.second;
  }


  message( "Applying boundary conditions to linear system" );

  // Modify linear system (A_i = 0)
  if ( not dolfin_get< bool >( "Krylov keep PC" ) )
  {
    A.zero( boundary_values.size(), dofs.data() );
  }

  // Finalise changes to A
  if ( not dolfin_get< bool >( "Krylov keep PC" ) )
    A.apply();

}
//-----------------------------------------------------------------------------
void DirichletBC::computeBCTopological( _map< size_t, real > & boundary_values,
                                        FiniteElementSpace const & space,
                                        SubSystem const &          sub_system )
{
  // Special case
  if ( entities_.empty() )
  {
    if ( !space.mesh().is_distributed() )
    {
      warning( "Found no facets matching domain for boundary condition." );
    }
    return;
  }

  // Iterate over facets
  DofMap const &        dof_map = space.dofmap();
  std::vector< size_t > cell_dofs( dof_map.num_element_dofs );
  ScratchSpace          scratch( space, sub_system );
  for ( std::vector< size_t >::const_iterator it = entities_.begin();
        it != entities_.end(); )
  {
    // Get cell number and local facet number
    size_t const c_index = *( it++ );
    size_t const f_index = *( it++ );
    size_t const c       = *( it++ );

    // Create cell
    Cell cell( mesh(), c_index );
    scratch.cell.update( cell );

    // Tabulate dofs on cell for the full space dofmap
    dof_map.tabulate_dofs( cell_dofs.data(), scratch.cell, cell );

    // Interpolate function on cell
    conditions[c].first.interpolate( scratch.coefficients.data(), scratch.cell,
                                     *scratch.finite_element, f_index );

    // Tabulate which dofs of the subdofmap are on the facet
    scratch.dof_map->tabulate_facet_dofs( scratch.facet_dofs.data(), f_index );

    // Pick values for facet
    for ( size_t i = 0; i < scratch.dof_map->num_facet_dofs(); ++i )
    {
      size_t const dof     = cell_dofs[scratch.offset + scratch.facet_dofs[i]];
      real const   value   = scratch.coefficients[scratch.facet_dofs[i]];
      boundary_values[dof] = value;
    }
  }
}
//-----------------------------------------------------------------------------
void DirichletBC::computeBCGeometric( _map< size_t, real > & boundary_values,
                                      FiniteElementSpace const & space,
                                      SubSystem const &          sub_system )
{
  // Special case
  if ( entities_.empty() )
  {
    if ( !space.mesh().is_distributed() )
    {
      warning( "Found no facets matching domain for boundary condition." );
    }
    return;
  }

  // Initialize facets, needed for geometric search
  size_t const facet_dim = mesh().type().facet_dim();

  // Iterate over facets
  DofMap const &        dof_map = space.dofmap();
  std::vector< size_t > cell_dofs( dof_map.num_element_dofs );
  ScratchSpace          scratch( space, sub_system );
  CellType * facet_type = mesh().type().create( mesh().type().facetType() );
  Point      xdof;
  boundary_values.reserve( entities_.size() );
  for ( std::vector< size_t >::const_iterator it = entities_.begin();
        it != entities_.end(); )
  {
    // Get cell number and local facet number
    size_t const c_index = *( it++ );
    size_t const f_index = *( it++ );
    size_t const c_      = *( it++ );

    // Create facet
    Cell  cell( mesh(), c_index );
    Facet facet( mesh(), cell.entities( facet_dim )[f_index] );

    // Loop the vertices associated with the facet
    for ( VertexIterator vertex( facet ); !vertex.end(); ++vertex )
    {
      // Loop the cells associated with the vertex
      for ( CellIterator c( *vertex ); !c.end(); ++c )
      {
        scratch.cell.update( *c );
        bool interpolated = false;

        // Tabulate dofs on cell for the full space dofmap
        dof_map.tabulate_dofs( cell_dofs.data(), scratch.cell, *c );

        // Tabulate coordinates of dofs on cell
        scratch.finite_element->tabulate_dof_coordinates(
          scratch.coordinates.data(), scratch.cell.coordinates.data() );

        // Loop over all dofs on cell
        for ( size_t i = 0; i < scratch.local_dimension; ++i )
        {
          // Copy coordinates to node Point
          xdof.set( &scratch.coordinates[i * Space::MAX_DIMENSION] );
          // Check if the coordinates are on current facet and thus on
          // boundary
          if ( !facet_type->intersects( facet, xdof ) )
          {
            continue;
          }

          if ( !interpolated )
          {
            interpolated = true;
            // Interpolate function on cell for the given (sub)element
            conditions[c_].first.interpolate( scratch.coefficients.data(),
                                              scratch.cell, *scratch.finite_element );
          }

          // Set boundary value
          size_t const dof     = cell_dofs[scratch.offset + i];
          real const   value   = scratch.coefficients[i];
          boundary_values[dof] = value;
        }
      }
    }
  }
  delete facet_type;
}
//-----------------------------------------------------------------------------
void DirichletBC::computeBCPointwise( _map< size_t, real > & boundary_values,
                                      FiniteElementSpace const & space,
                                      SubSystem const &          sub_system )
{
  // Iterate over cells
  DofMap const &        dof_map = space.dofmap();
  std::vector< size_t > cell_dofs( dof_map.num_element_dofs );
  ScratchSpace          scratch( space, sub_system );
  boundary_values.reserve( entities_.size() );

  for ( CellIterator cell( mesh() ); !cell.end(); ++cell )
  {
    scratch.cell.update( *cell );

    // Tabulate dofs on cell
    dof_map.tabulate_dofs( cell_dofs.data(), scratch.cell, *cell );

    // Tabulate coordinates of dofs on cell
    scratch.finite_element->tabulate_dof_coordinates(
      scratch.coordinates.data(), scratch.cell.coordinates.data() );

    // Interpolate function only once and only on cells where necessary
    bool interpolated = false;

    // Loop all dofs on cell
    for ( size_t i = 0; i < scratch.local_dimension; ++i )
    {
      // loop over all conditions in reverse order
      for ( size_t c = conditions.size(); c > 0; --c )
      {
        // Check if the coordinates are part of the sub domain
        if ( conditions[c - 1].second.inside(
               &scratch.coordinates[i * Space::MAX_DIMENSION], true ) )
        {
          if ( !interpolated )
          {
            interpolated = true;
            // Interpolate function on cell
            conditions[c - 1].first.interpolate( scratch.coefficients.data(),
                                                 scratch.cell, *scratch.finite_element );
            break;
          }
        }
      }

      // Set boundary value
      size_t const dof     = cell_dofs[scratch.offset + i];
      real const   value   = scratch.coefficients[i];
      boundary_values[dof] = value;
    }
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
