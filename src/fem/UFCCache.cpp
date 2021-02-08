// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/UFCCache.h>

#include <dolfin/fem/Form.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

UFCCache::~UFCCache()
{
  for ( size_t * d : dofs )
  {
    if ( d != nullptr )
    {
      delete[] d;
      d = nullptr;
    }
  }

  for ( size_t * md : macro_dofs )
  {
    if ( md != nullptr )
    {
      delete[] md;
      md = nullptr;
    }
  }

  for ( real * w_ : w )
  {
    if ( w_ != nullptr )
    {
      delete[] w_;
      w_ = nullptr;
    }
  }

  for ( real * mw : macro_w )
  {
    if ( mw != nullptr )
    {
      delete[] mw;
      mw = nullptr;
    }
  }
}

//-----------------------------------------------------------------------------

void UFCCache::init( Form const & form )
{
  cell.init( Cell( form.mesh(), 0 ) );

  macro_local_dimensions.resize( form().rank() );
  dofs.resize( form().rank() );
  macro_dofs.resize( form().rank() );
  w.resize( form().num_coefficients() );
  macro_w.resize( form().num_coefficients() );

  size_t num_entries_A       = 1;
  size_t num_entries_macro_A = 1;

  for ( size_t i = 0; i < form().rank(); ++i )
  {
    num_entries_A *= form.dofmaps()[i]->num_element_dofs;

    // Initialize local dimensions for macro element
    macro_local_dimensions[i] = 2 * form.dofmaps()[i]->num_element_dofs;
    num_entries_macro_A *= macro_local_dimensions[i];

    // Initialize dofs
    dofs[i] = new size_t[form.dofmaps()[i]->num_element_dofs];
    std::fill_n( dofs[i], form.dofmaps()[i]->num_element_dofs, 0 );

    // Initialize dofs on macro element
    macro_dofs[i] = new size_t[macro_local_dimensions[i]];
    std::fill_n( macro_dofs[i], macro_local_dimensions[i], 0 );
  }

  // Initialize local tensor
  A.resize( num_entries_A, 0. );

  // Initialize local tensor for macro element
  macro_A.resize( num_entries_macro_A, 0. );

  // Initialize coefficients
  for ( size_t i = 0; i < form().num_coefficients(); ++i )
  {
    // Create finite elements for coefficients
    w[i] = new real[form.elements()[form.rank()+i].space_dim];
    std::fill_n( w[i], form.elements()[form.rank()+i].space_dim, 0. );

    // Initialize coefficients on macro element
    macro_w[i] = new real[2 * form.elements()[form.rank()+i].space_dim];
    std::fill_n( macro_w[i], 2 * form.elements()[form.rank()+i].space_dim, 0. );
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
