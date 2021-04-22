// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/common/timing.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/fem/SparsityPatternBuilder.h>
#include <dolfin/fem/UFCHalo.h>
#include <dolfin/la/GenericTensor.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/log/log.h>
#include <dolfin/main/OpenMP.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>

#include <vector>

namespace dolfin
{

namespace Assembler
{

//-----------------------------------------------------------------------------

/// Assemble tensor from given (UFC) form, coefficients and sub domains.
/// This is the main assembly function in DOLFIN. All other assembly functions
/// end up calling this function.
///
/// The MeshFunction arguments can be used to specify assembly over subdomains
/// of the mesh cells, exterior facets and interior facets.
/// Either a null pointer or an empty MeshFunction may be used to specify that
/// the tensor should be assembled over the entire set of cells or facets.
void assemble( GenericTensor & A, Form & form,
               MeshValues< size_t, Cell > const *  cell_domains,
               MeshValues< size_t, Facet > const * exterior_facet_domains,
               MeshValues< size_t, Facet > const * interior_facet_domains,
               bool reset_tensor = true );

// Assemble over cells
void assembleCells( GenericTensor & A, Form & form,
                    MeshValues< size_t, Cell > const * domains );

// Assemble over exterior facets
void assembleExteriorFacets( GenericTensor & A, Form & form,
                             MeshValues< size_t, Facet > const * domains );

// Assemble over interior facets
void assembleInteriorFacets( GenericTensor & A, Form & form,
                             MeshValues< size_t, Facet > const * domains );

// Bogus-assemble periodic contributions
void initializePeriodicDofs( GenericTensor & A, Form & form );

// Initialize global tensor
void initGlobalTensor( GenericTensor & A, Form & form, bool reset_tensor );

//-----------------------------------------------------------------------------

void assemble( GenericTensor & A, Form & form, bool reset_tensor )
{
  OPENMP_PRAGMA( parallel )
  assemble( A, form, nullptr, nullptr, nullptr, reset_tensor );
}

//-----------------------------------------------------------------------------

void assemble( GenericTensor & A, Form & form,
               SubDomain const & sub_domain, bool reset_tensor )
{
  Mesh & mesh = form.mesh();

  // Extract cell domains
  MeshValues< size_t, Cell > * cell_domains = nullptr;

  // Extract facet domains
  MeshValues< size_t, Facet > * facet_domains = nullptr;

  OPENMP_PRAGMA( master )
  {
    if ( not form.cell_integrals().empty() )
    {
      cell_domains = new MeshValues< size_t, Cell >( mesh, 1 );
      sub_domain.mark( *cell_domains, 0 );
    }

    if ( not form.exterior_facet_integrals().empty()
         or not form.interior_facet_integrals().empty() )
    {
      facet_domains = new MeshValues< size_t, Facet >( mesh, 1 );
      sub_domain.mark( *facet_domains, 0 );
    }
  }

  // Assemble
  assemble( A, form, cell_domains, facet_domains, facet_domains, reset_tensor );

  // Delete domains
  OPENMP_PRAGMA( master )
  {
    delete cell_domains;
    delete facet_domains;
  }
}

//-----------------------------------------------------------------------------

void assemble( GenericTensor & A, Form & form,
               MeshValues< size_t, Cell > const &  cell_domains,
               MeshValues< size_t, Facet > const & exterior_facet_domains,
               MeshValues< size_t, Facet > const & interior_facet_domains,
               bool reset_tensor )
{
  assemble( A, form, &cell_domains,
            &exterior_facet_domains, &interior_facet_domains, reset_tensor );
}

//-----------------------------------------------------------------------------

void assemble( GenericTensor & A, Form & form,
               MeshValues< size_t, Cell > const *  cell_domains,
               MeshValues< size_t, Facet > const * exterior_facet_domains,
               MeshValues< size_t, Facet > const * interior_facet_domains,
               bool reset_tensor )
{
  // Check arguments
  OPENMP_PRAGMA( master )
  {
    if ( reset_tensor )
    {
      form.check( form.coefficients() );
    }
  }

  // Initialize global tensor
  OPENMP_PRAGMA( master )
  {
    initGlobalTensor( A, form, reset_tensor );

    // Update all ghost degrees of freedom
    for ( Coefficient * coeff : form.coefficients() )
    {
      coeff->sync();
    }
  }
  OPENMP_PRAGMA( flush )
  OPENMP_PRAGMA( barrier )

  // Assemble over cells
  assembleCells( A, form, cell_domains );

  // Assemble over exterior facets
  assembleExteriorFacets( A, form, exterior_facet_domains );

  // Assemble over interior facets
  assembleInteriorFacets( A, form, interior_facet_domains );

  // Bogus-assemble periodic dofs
  initializePeriodicDofs( A, form );

  // Finalise assembly of global tensor
  OPENMP_PRAGMA( master )
  A.apply();
  OPENMP_PRAGMA( barrier )
}

//-----------------------------------------------------------------------------

void assembleCells( GenericTensor & A, Form & form,
                    MeshValues< size_t, Cell > const * domains )
{
  if ( form.cell_integrals().empty() )
  {
    return;
  }

  message( 1, "Assembler: cells" );
  tic();

  Mesh &       mesh      = form.mesh();
  size_t const N         = mesh.num_cells();
  size_t const form_rank = form.rank();
  size_t const coef_size = form.coefficients().size();

  ufc::cell_integral const * integral  = form.cell_integrals().front();

  UFCCache & cache = form.cache();

  std::vector< size_t > local_dimensions( form.rank() );
  for ( size_t i = 0; i < form.rank(); ++i )
    local_dimensions[i] = form.dofmaps()[i]->num_element_dofs;

  CellIterator it( mesh );
OPENMP_PRAGMA( for )
for ( size_t i = 0; i < N; ++i )
{
  Cell & cell = it[i];

  // Get integral for sub domain (if any)
  if ( ( domains != nullptr ) && domains->size() > 0 )
  {
    size_t const domain = ( *domains )( cell );
    if ( domain < form.cell_integrals().size() )
    {
      integral = form.cell_integrals()[domain];
    }
    else
    {
      continue;
    }
  }

  // Update to current cell
  cache.cell.update( cell );

  // Interpolate coefficients on cell
  for ( size_t c = 0; c < coef_size; ++c )
  {
    form.coefficients()[c]->interpolate( cache.w[c], cache.cell,
                                         form.elements()[form_rank + c].ufc() );
  }

  // Tabulate dofs for each dimension
  for ( size_t d = 0; d < form_rank; ++d )
  {
    form.dofmaps()[d]->tabulate_dofs( cache.dofs[d], cache.cell );
  }

  // Tabulate cell tensor
  integral->tabulate_tensor( cache.A.data(), cache.w.data(),
                             cache.cell.coordinates.data(),
                             cache.cell.orientation );

  // Add entries to global tensor
  A.add( cache.A.data(), local_dimensions.data(), cache.dofs.data() );
}

tocd( 1 );
}

//-----------------------------------------------------------------------------

void assembleExteriorFacets( GenericTensor & A, Form & form,
                             MeshValues< size_t, Facet > const * domains )
{
  if ( form.exterior_facet_integrals().empty() )
  {
    return;
  }

  message( 1, "Assembler: exterior facets" );
  tic();

  Mesh &       mesh      = form.mesh();
  size_t const tdim      = mesh.topology_dimension();
  size_t const form_rank = form.rank();
  size_t const coef_size = form.coefficients().size();

  ufc::exterior_facet_integral const * integral = form.exterior_facet_integrals().front();

  UFCCache & cache = form.cache();

  BoundaryMesh & exterior_boundary = mesh.exterior_boundary();
  size_t const N = exterior_boundary.num_cells();
  if ( N == 0 )
  {
    return;
  }

  std::vector< size_t > local_dimensions( form.rank() );
  for ( size_t i = 0; i < form.rank(); ++i )
    local_dimensions[i] = form.dofmaps()[i]->num_element_dofs;

  FacetIterator it( mesh );
  CellIterator  c0( mesh );

OPENMP_PRAGMA( for )
for ( size_t i = 0; i < N; ++i )
{
  // Get mesh facet corresponding to boundary cell
  Facet & facet = it[exterior_boundary.facet_index( i )];

  // Get integral for sub domain (if any)
  if ( ( domains != nullptr ) && domains->size() > 0 )
  {
    size_t const domain = ( *domains )( facet );
    if ( domain < form.exterior_facet_integrals().size() )
    {
      integral = form.exterior_facet_integrals()[domain];
    }
    else
    {
      continue;
    }
  }

  dolfin_assert( integral != nullptr );

  // Get mesh cell to which mesh facet belongs
  Cell & cell = c0[facet.entities( tdim )[0]];

  // Get local index of facet with respect to the cell
  size_t const local_facet = cell.index( facet );

  // Update to current cell
  cache.cell.update( cell );

  // Interpolate coefficients on cell
  for ( size_t c = 0; c < coef_size; ++c )
  {
    form.coefficients()[c]->interpolate( cache.w[c], cache.cell,
                                         form.elements()[form.rank() + c].ufc(),
                                         local_facet );
  }

  // Tabulate dofs for each dimension
  for ( size_t d = 0; d < form_rank; ++d )
  {
    form.dofmaps()[d]->tabulate_dofs( cache.dofs[d], cache.cell );
  }

  // Tabulate exterior facet tensor
  integral->tabulate_tensor( cache.A.data(), cache.w.data(),
                             cache.cell.coordinates.data(),
                             local_facet, cache.cell.orientation );

  // Add entries to global tensor
  A.add( cache.A.data(), local_dimensions.data(), cache.dofs.data() );
}

tocd( 1 );
}

//-----------------------------------------------------------------------------

void assembleInteriorFacets( GenericTensor & A, Form & form,
                             MeshValues< size_t, Facet > const * domains )
{
  if ( form.interior_facet_integrals().empty() )
  {
    return;
  }

  message( 1, "Assembler: interior facets" );
  tic();

  Mesh &       mesh      = form.mesh();
  size_t const tdim      = mesh.topology_dimension();
  size_t const N         = mesh.size( mesh.type().facet_dim() );
  size_t const form_rank = form.rank();
  size_t const coef_size = form.coefficients().size();

  ufc::interior_facet_integral const * integral = form.interior_facet_integrals().front();

  UFCCache & cache = form.cache();

  // Halo data structure caching macro element coefficients and dofs
  UFCHalo halo( form );

  FacetIterator it( mesh );
  CellIterator  c0( mesh );
  CellIterator  c1( mesh );

OPENMP_PRAGMA( for )
for ( size_t i = 0; i < N; ++i )
{
  Facet & facet = it[i];

  // Get integral for sub domain (if any)
  if ( ( domains != nullptr ) && domains->size() > 0 )
  {
    size_t const domain = ( *domains )( facet );
    if ( domain < form.interior_facet_integrals().size() )
    {
      integral = form.interior_facet_integrals()[domain];
    }
    else
    {
      continue;
    }
  }

  // Check if we have a local interior facet
  if ( facet.num_entities( tdim ) == 2 )
  {
    // Get cells incident with facet, local index of facet
    Cell & cell0 = c0[facet.entities( tdim )[0]];
    halo.facet0  = cell0.index( facet );
    halo.cell0.update( cell0 );

    Cell & cell1 = c1[facet.entities( tdim )[1]];
    halo.facet1  = cell1.index( facet );
    halo.cell1.update( cell1 );

    // Interpolate coefficients on cell1
    for ( size_t c = 0; c < coef_size; ++c )
    {
      FiniteElement const & fe = form.elements()[form.rank() + c];

      form.coefficients()[c]->interpolate( cache.macro_w[c], halo.cell0,
                                           fe.ufc(), halo.facet0 );
      size_t const offset = fe.space_dim;
      form.coefficients()[c]->interpolate( cache.macro_w[c] + offset, halo.cell1,
                                           fe.ufc(), halo.facet1 );
    }

    // Tabulate dofs for each dimension on cell1
    for ( size_t d = 0; d < form_rank; ++d )
    {
      form.dofmaps()[d]->tabulate_dofs( cache.macro_dofs[d], halo.cell0 );
      size_t const offset = form.dofmaps()[d]->num_element_dofs;
      form.dofmaps()[d]->tabulate_dofs( cache.macro_dofs[d] + offset,
                                        halo.cell1 );
    }
  }
  // Interprocess facet
  else if ( facet.is_shared() )
  {
    // Contributions from cell0 are restored from the halo data while
    // contributions from cell1 are fetched from adjacent ranks.
    // Implementation updates pointers to coordinates, but has to copy dofs
    // and coefficients until data structures are reworked.
    halo.update( facet );
  }

  integral->tabulate_tensor( cache.macro_A.data(), cache.macro_w.data(),
                             halo.cell0.coordinates.data(),
                             halo.cell1.coordinates.data(),
                             halo.facet0, halo.facet1,
                             halo.cell0.orientation, halo.cell1.orientation );

  // Add entries to global tensor
  A.add( cache.macro_A.data(),
         cache.macro_local_dimensions.data(),
         cache.macro_dofs.data() );
}

tocd( 1 );
}

//-----------------------------------------------------------------------------

void initializePeriodicDofs( GenericTensor & A, Form & form )
{
  if ( not form.mesh().has_periodic_constraint() )
  {
    return;
  }

  // Add zero at periodic dofs to allocate entries
  /// @todo This could be fixed by a modification of the assembler's behaviour
  if ( A.rank() == 2 )
  {
    PeriodicDofsMapping const & pdm =
      form.dofmaps()[0]->periodic_mapping( form.spaces()[0] );

    std::vector< real >   block( pdm.max_local_dimension() + 1, 0.0 );
    std::vector< size_t > jcols( pdm.max_local_dimension() + 1, 0 );

    Matrix & matA  = static_cast< Matrix & >( A );
    size_t   irow  = 0;
    size_t   ncols = 0;

    for ( size_t i = 0; i < pdm.num_Gdofs(); ++i )
    {
      pdm.tabulate_dofs( i, &irow, jcols.data(), ncols );
      jcols[ncols] = irow;
      matA.add( block.data(), 1, &irow, ncols, jcols.data() );
    }
  }
}

//-----------------------------------------------------------------------------

void initGlobalTensor( GenericTensor & A, Form & form, bool reset_tensor )
{
  if ( A.rank() == 0 )
  {
    if ( reset_tensor )
    {
      A.zero();
    }

    return;
  }

  if ( reset_tensor || A.size( 0 ) == 0 )
  {
    GenericSparsityPattern * sparsity_pattern = A.factory().createPattern();
    SparsityPatternBuilder::build( *sparsity_pattern, form );
    A.init( *sparsity_pattern );

    // initialize ghost values for vectors
    if ( A.rank() == 1 )
    {
      _ordered_set< size_t > indices;

      for ( CellIterator cell( form.mesh() ); !cell.end(); ++cell )
      {
        // Update to current cell
        form.cache().cell.update( *cell );

        // Tabulate dofs
        form.spaces()[0].dofmap().tabulate_dofs( form.cache().dofs[0],
                                                 form.cache().cell );

        for ( size_t j = 0; j < form.spaces()[0].element().space_dim; ++j )
        {
          indices.insert( form.cache().dofs[0][j] );
        }
      }

      _ordered_map< size_t, size_t > empty;
      A.down_cast< GenericVector >().init_ghosted( indices.size(), indices, empty );
    }

    delete sparsity_pattern;
  }
  else
  {
    if ( ( A.rank() > 0 && A.size( 0 ) != form.dofmaps()[0]->global_dim )
         || ( A.rank() > 1 && A.size( 1 ) != form.dofmaps()[1]->global_dim ) )
    {
      error( "Assembler : dimensions of linear system do not match spaces." );
    }

    A.zero();
  }
}

//-----------------------------------------------------------------------------

} // namespace Assembler

} // namespace dolfin
