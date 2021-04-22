// Copyright (C) 2007-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/SparsityPatternBuilder.h>

#include <dolfin/common/timing.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/la/GenericSparsityPattern.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

namespace SparsityPatternBuilder
{

//-----------------------------------------------------------------------------
void build( GenericSparsityPattern & sparsity_pattern,
            Form &                   form)
{
  message( 1, "SparsityPatternBuilder: build" );
  tic();

  Mesh &     mesh  = form.mesh();
  UFCCache & cache = form.cache();

  std::vector< size_t > local_dimensions( form.rank() );
  std::vector< size_t > local_sizes( form.rank() );
  std::vector< size_t > global_dimensions( form.rank() );

  for ( size_t i = 0; i < form.rank(); ++i )
  {
    local_dimensions[i]  = form.dofmaps()[i]->num_element_dofs;
    local_sizes[i]       = form.dofmaps()[i]->local_size();
    global_dimensions[i] = form.dofmaps()[i]->global_dim;
  }

  // Initialise sparsity pattern
  sparsity_pattern.init( form.rank(),
                         global_dimensions.data(),
                         local_sizes.data() );

  // Only build for rank >= 2 (matrices and higher order tensors)
  if ( form.rank() < 2 )
  {
    tocd( 1 );
    return;
  }

  // JANPACK doesn't need any sparsity pattern information
  if ( dolfin_get< std::string >( "linear algebra backend" ) == "JANPACK" )
  {
    tocd( 1 );
    return;
  }

  // Build sparsity pattern for cell integrals
  if ( not form.cell_integrals().empty() )
  {
    for ( CellIterator cell( mesh ); !cell.end(); ++cell )
    {
      // Update to current cell
      cache.cell.update( *cell );

      // Tabulate dofs for each dimension
      for ( size_t i = 0; i < form.rank(); ++i )
      {
        form.dofmaps()[i]->tabulate_dofs( cache.dofs[i], cache.cell );
      }

      // Fill sparsity pattern.
      sparsity_pattern.insert( local_dimensions.data(), cache.dofs.data() );
    }
  }

  // Note: no need to iterate over exterior facets since those dofs
  // are included when tabulating dofs on all cells

  // Build sparsity pattern for interior facet integrals
  if ( not form.interior_facet_integrals().empty() )
  {
    size_t const tdim = mesh.topology_dimension();

    Cell c( mesh, 0 );
    UFCCell cell0( c );
    UFCCell cell1( c );

    for ( FacetIterator facet( mesh ); !facet.end(); ++facet )
    {
      // Check if we have an interior facet
      if ( facet->num_entities( tdim ) != 2 )
      {
        continue;
      }

      // Get cells incident with facet
      Cell cell0_( mesh, facet->entities( tdim )[0] );
      Cell cell1_( mesh, facet->entities( tdim )[1] );

      // Update to current pair of cells
      cell0.update( cell0_ );
      cell1.update( cell1_ );

      // Tabulate dofs for each dimension on macro element
      for ( size_t i = 0; i < form.rank(); ++i )
      {
        const size_t offset = form.dofmaps()[i]->num_element_dofs;
        form.dofmaps()[i]->tabulate_dofs( cache.macro_dofs[i], cell0 );
        form.dofmaps()[i]->tabulate_dofs( cache.macro_dofs[i] + offset,
                                          cell1 );
      }

      // Fill sparsity pattern.
      sparsity_pattern.insert( cache.macro_local_dimensions.data(),
                               cache.macro_dofs.data() );
    }
  }

  // Build sparsity pattern for periodic facets
  // Only for square systems
  if ( mesh.has_periodic_constraint() )
  {
    bool has_facet_dofs = ( form.dofmaps()[0]->num_facet_dofs > 0 );
    bool is_square      = true;
    for ( size_t i = 1; i < form.rank(); ++i )
    {
      has_facet_dofs |= ( form.dofmaps()[i]->num_facet_dofs > 0 );
      is_square &= ( *form.dofmaps()[i] == *form.dofmaps()[i - 1] );
    }

    if ( has_facet_dofs )
    {
      PeriodicDofsMapping const & pdm =
        form.dofmaps()[0]->periodic_mapping( form.spaces()[0] );
      size_t local_dim[2];
      local_dim[0]   = 1;
      local_dim[1]   = pdm.max_local_dimension();
      size_t ** dofs = new size_t*[2];
      dofs[0]        = new size_t[1];
      dofs[1]        = new size_t[pdm.max_local_dimension()];
      for ( size_t i = 0; i < pdm.num_Gdofs(); ++i )
      {
        pdm.tabulate_dofs( i, dofs[0], dofs[1], local_dim[1] );

        // Fill sparsity pattern.
        sparsity_pattern.insert( local_dim, dofs );
      }
      delete[] dofs[1];
      delete[] dofs[0];
      delete[] dofs;
    }
  }

  // Finalize sparsity pattern
  sparsity_pattern.apply();

  tocd( 1 );
}
//-----------------------------------------------------------------------------

} // end namespace SparsityPatternBuilder

} /* namespace dolfin */
