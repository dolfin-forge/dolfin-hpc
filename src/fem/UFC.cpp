// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/UFC.h>

#include <dolfin/fem/DofMapSet.h>
#include <dolfin/fem/Form.h>
#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFC::UFC( Form const & form )
  : form( form )
  , finite_elements( form.rank() )
  , coefficient_elements( form.num_coefficients() )
  , cell_integrals( std::max( form.max_cell_subdomain_id(), 1ul ) )
  , exterior_facet_integrals(
      std::max( form.max_exterior_facet_subdomain_id(), 1ul ) )
  , interior_facet_integrals(
      std::max( form.max_interior_facet_subdomain_id(), 1ul ) )
  , cell( Cell( form.mesh(), 0 ) )
  , cell0( Cell( form.mesh(), 0 ) )
  , cell1( Cell( form.mesh(), 0 ) )
  , facet0( 0 )
  , facet1( 0 )
  , local_dimensions( form.rank() )
  , macro_local_dimensions( form.rank() )
  , local_sizes( form.rank() )
  , global_dimensions( form.rank() )
  , dofs( form.rank() )
  , macro_dofs( form.rank() )
  , w( form.num_coefficients() )
  , macro_w( form.num_coefficients() )
{
  // Create cell integrals
  if ( form.max_cell_subdomain_id() != 0 )
  {
    for ( size_t i = 0; i < form.max_cell_subdomain_id(); ++i )
      cell_integrals[i] = form.create_cell_integral( i );
  }
  else
  {
    cell_integrals[0] = form.create_default_cell_integral();
  }

  // Create exterior facet integrals
  if ( form.max_exterior_facet_subdomain_id() != 0 )
  {
    for ( size_t i = 0; i < form.max_exterior_facet_subdomain_id(); ++i )
      exterior_facet_integrals[i] = form.create_exterior_facet_integral( i );
  }
  else
  {
    exterior_facet_integrals[0] = form.create_default_exterior_facet_integral();
  }

  // Create interior facet integrals
  if ( form.max_interior_facet_subdomain_id() != 0 )
  {
    for ( size_t i = 0; i < form.max_interior_facet_subdomain_id(); ++i )
      interior_facet_integrals[i] = form.create_interior_facet_integral( i );
  }
  else
  {
    interior_facet_integrals[0] = form.create_default_interior_facet_integral();
  }

  size_t num_entries_A       = 1;
  size_t num_entries_macro_A = 1;

  for ( size_t i = 0; i < form.rank(); ++i )
  {
    // Create finite elements
    finite_elements[i] = form.create_finite_element( i );

    // Initialize local dimensions
    local_dimensions[i] = form.dofmaps()[i].num_element_dofs();
    num_entries_A *= local_dimensions[i];

    // Initialize local dimensions for macro element
    macro_local_dimensions[i] = 2 * form.dofmaps()[i].num_element_dofs();
    num_entries_macro_A *= macro_local_dimensions[i];

    // Initialize local sizes
    local_sizes[i] = form.dofmaps()[i].local_size();

    dolfin_assert( form.mesh().num_entities().size()
                   == finite_elements[i]->topological_dimension() + 1 );

    // Initialize global dimensions
    global_dimensions[i] =
      form.dofmaps()[i].global_dimension( form.mesh().num_entities() );

    // Initialize dofs
    dofs[i] = new size_t[local_dimensions[i]];
    std::fill_n( dofs[i], local_dimensions[i], 0 );

    // Initialize dofs on macro element
    macro_dofs[i] = new size_t[macro_local_dimensions[i]];
    std::fill_n( macro_dofs[i], macro_local_dimensions[i], 0 );
  }

  // Initialize local tensor
  A.resize( num_entries_A, 0. );

  // Initialize local tensor for macro element
  macro_A.resize( num_entries_macro_A, 0. );

  // Initialize coefficients
  for ( size_t i = 0; i < form.num_coefficients(); ++i )
  {
    // Create finite elements for coefficients
    coefficient_elements[i] = form.create_finite_element( form.rank() + i );
    w[i] = new real[coefficient_elements[i]->space_dimension()];
    std::fill_n( w[i], coefficient_elements[i]->space_dimension(), 0. );

    // Initialize coefficients on macro element
    macro_w[i] = new real[2 * coefficient_elements[i]->space_dimension()];
    std::fill_n(
      macro_w[i], 2 * coefficient_elements[i]->space_dimension(), 0. );
  }
}

//-----------------------------------------------------------------------------

UFC::~UFC()
{
  destruct( finite_elements );
  destruct( coefficient_elements );
  destruct( cell_integrals );
  destruct( exterior_facet_integrals );
  destruct( interior_facet_integrals );
  destruct( dofs );
  destruct( macro_dofs );
  destruct( w );
  destruct( macro_w );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
