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
{
  init( form, form.mesh(), form.dofmaps() );
}
//-----------------------------------------------------------------------------
void UFC::init( ufc::form const & form, Mesh & mesh, DofMapSet const & dofmaps )
{
  // Create finite elements
  finite_elements = new ufc::finite_element *[form.rank()];
  for ( uint i = 0; i < form.rank(); ++i )
  {
    finite_elements[i] = form.create_finite_element( i );
  }

  // Create finite elements for coefficients
  coefficient_elements = new ufc::finite_element *[form.num_coefficients()];
  for ( uint i = 0; i < form.num_coefficients(); ++i )
  {
    coefficient_elements[i] = form.create_finite_element( form.rank() + i );
  }

  // Create cell integrals
  cell_integrals = new ufc::cell_integral *[form.max_cell_subdomain_id()];
  for ( uint i = 0; i < form.max_cell_subdomain_id(); ++i )
  {
    cell_integrals[i] = form.create_cell_integral( i );
  }

  // Create exterior facet integrals
  exterior_facet_integrals =
    new ufc::exterior_facet_integral *[form.max_exterior_facet_subdomain_id()];
  for ( uint i = 0; i < form.max_exterior_facet_subdomain_id(); ++i )
  {
    exterior_facet_integrals[i] = form.create_exterior_facet_integral( i );
  }

  // Create interior facet integrals
  interior_facet_integrals =
    new ufc::interior_facet_integral *[form.max_interior_facet_subdomain_id()];
  for ( uint i = 0; i < form.max_interior_facet_subdomain_id(); ++i )
  {
    interior_facet_integrals[i] = form.create_interior_facet_integral( i );
  }

  // Initialize cells with first cell in mesh
  if ( mesh.num_cells() > 0 )
  {
    Cell cell( mesh, 0 );
    this->cell.init( cell );
    this->cell0.init( cell );
    this->cell1.init( cell );
  }
  facet0 = 0;
  facet1 = 0;

  // Initialize local tensor
  uint num_entries = 1;
  for ( uint i = 0; i < form.rank(); ++i )
  {
    num_entries *= dofmaps[i].num_element_support_dofs();
  }
  A = new real[num_entries];
  std::fill_n( A, num_entries, 0. );

  // Initialize local tensor for macro element
  num_entries = 1;
  for ( uint i = 0; i < form.rank(); ++i )
  {
    num_entries *= 2 * dofmaps[i].num_element_support_dofs();
  }
  macro_A = new real[num_entries];
  std::fill_n( macro_A, num_entries, 0. );

  // Initialize local dimensions
  local_dimensions = new uint[form.rank()];
  for ( uint i = 0; i < form.rank(); ++i )
  {
    local_dimensions[i] = dofmaps[i].num_element_support_dofs();
  }

  // Initialize local dimensions for macro element
  macro_local_dimensions = new uint[form.rank()];
  for ( uint i = 0; i < form.rank(); ++i )
  {
    macro_local_dimensions[i] = 2 * dofmaps[i].num_element_support_dofs();
  }

  // Initialize local sizes
  local_sizes = new uint[form.rank()];
  for ( uint i = 0; i < form.rank(); ++i )
  {
    local_sizes[i] = dofmaps[i].local_size();
  }

  // Initialize global dimensions
  global_dimensions = new uint[form.rank()];
  for ( uint i = 0; i < form.rank(); ++i )
  {
    // FIXME is i the correct here?
    ufc::finite_element * fe = form.create_finite_element( i );

    // FIXME num_entities should maybe be stored somewhere else
    std::vector< size_t > num_entities( fe->topological_dimension(),
                                        0 );
    for ( uint d = 0; d <= fe->topological_dimension(); ++d )
    {
      if ( mesh.topology().connectivity( d ) )
      {
        num_entities[d] = mesh.topology().global_size( d );
      }
    }

    global_dimensions[i] = dofmaps[i].global_dimension( num_entities );

    delete fe;
  }

  // Initialize dofs
  dofs = new uint *[form.rank()];
  for ( uint i = 0; i < form.rank(); ++i )
  {
    dofs[i] = new uint[local_dimensions[i]];
    std::fill_n( dofs[i], local_dimensions[i], 0 );
  }

  // Initialize dofs on macro element
  macro_dofs = new uint *[form.rank()];
  for ( uint i = 0; i < form.rank(); ++i )
  {
    macro_dofs[i] = new uint[macro_local_dimensions[i]];
    std::fill_n( macro_dofs[i], macro_local_dimensions[i], 0 );
  }

  // Initialize coefficients
  w = new real *[form.num_coefficients()];
  for ( uint i = 0; i < form.num_coefficients(); ++i )
  {
    w[i] = new real[coefficient_elements[i]->space_dimension()];
    std::fill_n( w[i], coefficient_elements[i]->space_dimension(), 0. );
  }

  // Initialize coefficients on macro element
  macro_w = new real *[form.num_coefficients()];
  for ( uint i = 0; i < form.num_coefficients(); ++i )
  {
    macro_w[i] = new real[2 * coefficient_elements[i]->space_dimension()];
    std::fill_n(
      macro_w[i], 2 * coefficient_elements[i]->space_dimension(), 0. );
  }
}
//-----------------------------------------------------------------------------
UFC::~UFC()
{
  // Delete finite elements
  for ( uint i = 0; i < form.rank(); ++i )
  {
    delete finite_elements[i];
  }
  delete[] finite_elements;

  // Delete coefficient finite elements
  for ( uint i = 0; i < form.num_coefficients(); ++i )
  {
    delete coefficient_elements[i];
  }
  delete[] coefficient_elements;

  // Delete cell integrals
  for ( uint i = 0; i < form.max_cell_subdomain_id(); ++i )
  {
    delete cell_integrals[i];
  }
  delete[] cell_integrals;

  // Delete exterior facet integrals
  for ( uint i = 0; i < form.max_exterior_facet_subdomain_id(); ++i )
  {
    delete exterior_facet_integrals[i];
  }
  delete[] exterior_facet_integrals;

  // Delete interior facet integrals
  for ( uint i = 0; i < form.max_interior_facet_subdomain_id(); ++i )
  {
    delete interior_facet_integrals[i];
  }
  delete[] interior_facet_integrals;

  // Delete local tensor
  delete[] A;

  // Delete local tensor for macro element
  delete[] macro_A;

  // Delete local dimensions
  delete[] local_dimensions;

  // Delete local sizes
  delete[] local_sizes;

  // Delete global dimensions
  delete[] global_dimensions;

  // Delete local dimensions for macro element
  delete[] macro_local_dimensions;

  // Delete dofs
  for ( uint i = 0; i < form.rank(); ++i )
  {
    delete[] dofs[i];
  }
  delete[] dofs;

  // Delete macro dofs
  for ( uint i = 0; i < form.rank(); ++i )
  {
    delete[] macro_dofs[i];
  }
  delete[] macro_dofs;

  // Delete coefficients
  for ( uint i = 0; i < form.num_coefficients(); ++i )
  {
    delete[] w[i];
  }
  delete[] w;

  // Delete macro coefficients
  for ( uint i = 0; i < form.num_coefficients(); ++i )
  {
    delete[] macro_w[i];
  }
  delete[] macro_w;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
