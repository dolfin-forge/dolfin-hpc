// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/PeriodicBC.h>

#include <dolfin/common/constants.h>
#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/PeriodicSubDomain.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>

#include <vector>

namespace dolfin
{

// Comparison operator for hashing coordinates. Note that two
// coordinates are considered equal if equal to within round-off.
struct lt_coordinate
{
  auto operator()( const std::vector< real > & x,
                   const std::vector< real > & y ) const -> bool
  {
    if ( x.size() > ( y.size() + DOLFIN_EPS ) )
      return false;

    for ( size_t i = 0; i < x.size(); ++i )
    {
      if ( x[i] < ( y[i] - DOLFIN_EPS ) )
        return true;
      else if ( x[i] > ( y[i] + DOLFIN_EPS ) )
        return false;
    }

    return false;
  }
};

//-----------------------------------------------------------------------------

PeriodicBC::PeriodicBC( Mesh & mesh, PeriodicSubDomain const & sub_domain )
  : BoundaryCondition( "Periodic", mesh, sub_domain )
{
}

//-----------------------------------------------------------------------------

PeriodicBC::PeriodicBC( Mesh &                    mesh,
                        PeriodicSubDomain const & sub_domain,
                        SubSystem const &         sub_system )
  : BoundaryCondition( "Periodic", mesh, sub_domain, sub_system )
{
}

//-----------------------------------------------------------------------------

void PeriodicBC::apply( GenericMatrix &      A,
                        GenericVector &      b,
                        BilinearForm const & form )
{
  message( "Applying periodic boundary conditions to linear system." );

  PeriodicSubDomain const & subdomain =
    static_cast< PeriodicSubDomain const & >( this->sub_domain() );

  /// @todo Make this work for non-scalar subsystems, like vector-valued
  /// Lagrange where more than one per element is associated with
  /// each coordinate. Note that globally there may very well be
  /// more than one dof per coordinate (for conforming elements).
  size_t const gdim = mesh().geometry_dimension();
  size_t const tdim = mesh().geometry_dimension();

  // Table of mappings from coordinates to dofs
  using coord2dofs_mapping_t = _ordered_map< std::vector< real >,
                                             std::pair< size_t, size_t >,
                                             lt_coordinate >;

  using coords2dofs_it = coord2dofs_mapping_t::iterator;


  coord2dofs_mapping_t coord_dofs;

  // vectors used for mapping coordinates
  std::vector< real > xx( gdim, 0.0 );
  std::vector< real > y( gdim, 0.0 );

  // Make sure we have the facet - cell connectivity
  mesh().init( mesh().type().facet_dim(), tdim );

  // Create local data for application of boundary conditions
  FiniteElementSpace const & space = form.trial_space();
  ScratchSpace               scratch( space );

  // Iterate over the facets of the mesh
  for ( FacetIterator facet( mesh() ); !facet.end(); ++facet )
  {
    // Get cell to which facet belongs (there may be two, but pick first)
    Cell cell( mesh(), facet->entities( tdim )[0] );
    scratch.cell.update( cell );

    // Get local index of facet with respect to the cell
    size_t const local_facet = cell.index( *facet );

    // Tabulate dofs on cell
    space.dofmap().tabulate_dofs( scratch.dofs.data(), scratch.cell );

    // Tabulate coordinates on cell
    scratch.finite_element->tabulate_dof_coordinates(
      scratch.coordinates.data(), scratch.cell.coordinates.data() );

    // Tabulate which dofs are on the facet
    scratch.dof_map->tabulate_facet_dofs( scratch.facet_dofs.data(),
                                          local_facet );

    // Iterate over facet dofs
    for ( size_t i = 0; i < scratch.dof_map->num_facet_dofs(); ++i )
    {
      // Get dof and coordinate of dof
      size_t const dof        = scratch.offset + scratch.facet_dofs[i];
      size_t const global_dof = scratch.dofs[dof];
      real *       x = &scratch.coordinates[scratch.facet_dofs[i] * gdim];

      // Map coordinate from H to G
      std::copy( x, x + gdim, y.begin() );
      subdomain.map( x, y.data() );

      // Check if coordinate is inside the domain G or in H
      bool const on_boundary = facet->num_entities( tdim ) == 1;

      // G dof
      if ( subdomain.inside( x, on_boundary ) )
      {
        // Copy coordinate to std::vector
        std::copy( x, x + gdim, xx.begin() );

        // Check if coordinate exists from before
        coords2dofs_it it = coord_dofs.find( xx );
        if ( it != coord_dofs.end() )
        {
          // // Check that we don't have more than one dof per coordinate
          // if ( it->second.first != DOLFIN_SIZE_T_UNDEF )
          // {
          //   cout << "Coordinate: x =";
          //   for ( size_t j = 0; j < gdim; j++ )
          //     cout << " " << xx[j];
          //   cout << "\nDegrees of freedom: ";
          //   cout << it->second.first << " " << global_dof << "\n";
          //   error( "More than one dof associated with coordinate. Did you "
          //          "forget to specify the subsystem?" );
          // }
          it->second.first = global_dof;
        }
        else
        {
          // Doesn't exist, so create new pair (with illegal second value)
          coord_dofs[xx] = std::make_pair( global_dof, DOLFIN_SIZE_T_UNDEF );
        }
      }
      // H dof
      else if ( subdomain.inside( y.data(), on_boundary ) )
      {
        // y = F(x) is in G, so coordinate x is in H

        // Copy coordinate to std::vector
        std::copy( y.data(), y.data() + gdim, xx.begin() );

        // Check if coordinate exists from before
        coords2dofs_it it = coord_dofs.find( xx );
        if ( it != coord_dofs.end() )
        {
          // Check that we don't have more than one dof per coordinate
          // if ( it->second.second != DOLFIN_SIZE_T_UNDEF )
          // {
          //   cout << "Coordinate: x =";
          //   for ( size_t j = 0; j < gdim; j++ )
          //     cout << " " << xx[j];
          //   cout << "\nDegrees of freedom: ";
          //   cout << it->second.second << " " << global_dof << "\n";
          //   error( "More than one dof associated with coordinate. Did you "
          //          "forget to specify the subsystem?" );
          // }
          it->second.second = global_dof;
        }
        else
        {
          // Doesn't exist, so create new pair (with illegal first value)
          coord_dofs[xx] = std::make_pair( DOLFIN_SIZE_T_UNDEF, global_dof );
        }
      }
    }
  }

  // Given pairs <dofG, dofH>
  // Insert 1 at (dof0, dof0)
  std::vector< size_t > rows( coord_dofs.size() );
  size_t                i = 0;
  for ( coords2dofs_it it = coord_dofs.begin(); it != coord_dofs.end(); ++it )
  {
    rows[i++] = it->second.first;
  }
  A.ident( coord_dofs.size(), rows.data() );

  // Insert -1 at (dof0, dof1) and 0 on right-hand side
  real const vals = -1.0;
  real const zero = 0.0;

  for ( coords2dofs_it it = coord_dofs.begin(); it != coord_dofs.end(); ++it )
  {
    // Check that we got both dofs
    size_t const & dof0 = it->second.first;
    size_t const & dof1 = it->second.second;

    if ( dof0 == DOLFIN_SIZE_T_UNDEF or dof1 == DOLFIN_SIZE_T_UNDEF )
    {
      cout << "At coordinate: x =";
      for ( size_t j = 0; j < gdim; j++ )
        cout << " " << it->first[j];
      cout << "\n";
      error( "PeriodicBC : Unable to find a pair of matching dofs for "
             "periodic boundary condition." );
    }

    /// @todo Perhaps this can be done more efficiently?

    // Set x_i - x_j = 0
    A.set( &vals, 1, &dof0, 1, &dof1 );
    b.set( &zero, 1, &dof0 );
  }

  // Apply changes
  A.apply();
  b.apply();

  /*
  //For the moment just test preallocation
  if(A.rank() == 2)
  {
    //
    Matrix& matA = static_cast<Matrix&>(A);
    PeriodicDofsMapping const& pdm = dof_map_set[0].periodic_mapping();
    real * block = new real[pdm.max_local_dimension() + 1];
    std::fill_n(block, pdm.max_local_dimension() + 1, 0.0);
    size_t irow = 0;
    size_t * jcols = new size_t[pdm.max_local_dimension() + 1];
    std::fill_n(jcols, pdm.max_local_dimension() + 1, 0.0);
    size_t ncols = 0;
    for (size_t i = 0; i < pdm.num_Gdofs(); ++i)
    {
      pdm.tabulate_dofs(i, &irow, jcols, ncols);
      jcols[ncols] = irow;
      matA.add(block, 1, &irow, ncols, jcols);
    }
    delete[] jcols;
    delete[] block;
  }
  else if(A.rank() == 1)
  {
    //
    PeriodicDofsMapping const& pdm = dof_map_set[0].periodic_mapping();
    Vector& vecA = static_cast<Vector&>(A);
    size_t const num_rows = pdm.num_Gdofs();
    real * block = new real[num_rows];
    std::memset(&block[0], 0, num_rows);
    vecA.set(&block[0], num_rows, pdm.get_Gindices());
    delete [] block;
  }
  */
}

//-----------------------------------------------------------------------------

void PeriodicBC::apply( GenericMatrix &,
                        GenericVector &,
                        const GenericVector &,
                        BilinearForm const & )
{
  error( "PeriodicBC : not implemented for nonlinear systems." );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
