// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/Form.h>

#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/CoefficientMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Function.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Form::Form( Mesh & mesh )
  : mesh_( mesh )
  , dof_map_set_( *this, mesh )
{
}

//-----------------------------------------------------------------------------
auto Form::check( std::vector< Coefficient * > const & coefficients ) const
  -> bool
{
  // Check that we get the correct number of coefficients
  if ( coefficients.size() != this->form().num_coefficients() )
  {
    error( "Incorrect number of coefficients: %d given but %d required.",
           coefficients.size(),
           this->form().num_coefficients() );
  }

  // Check that all coefficients have valid value dimensions
  for ( size_t i = 0; i < coefficients.size(); ++i )
  {
    message( 1,
             "Form: Checking coefficient %d: %s",
             i,
             this->coefficient_name( i ).c_str() );
    if ( coefficients[i] == nullptr )
    {
      error( "Got nullptr pointer as coefficient %d labeled as '%s'.",
             i,
             this->coefficient_name( i ).c_str() );
    }

    ufc::finite_element * fe = this->form().create_finite_element( i + this->form().rank() );
    Function * fptr = dynamic_cast< Function * >( this->coefficients()[i] );
    if ( fptr != nullptr )
    {
      if ( fptr->empty() )
      {
        error( "Coefficient %i is empty", i );
      }
      else if ( strcmp( fptr->space().element()().signature(), fe->signature() )
                != 0 )
      {
        error( "Mismatch of discrete space for Coefficient %i", i );
      }
    }
    else
    {
      size_t coef_rank = coefficients[i]->rank();
      size_t fe_rank   = fe->value_rank();
      message( 1,
               "Form: Coefficient rank: expected  = %d, provided = %d, ",
               fe_rank,
               coef_rank );
      if ( fe_rank != coef_rank )
      {
        error(
          "Invalid value rank of Coefficient '%s' with index %d:\n"
          "Got %d but expecting %d.\n"
          "You may need to provide the rank of a user defined Coefficient.",
          this->coefficient_name( i ).c_str(),
          i,
          coef_rank,
          fe_rank );
      }

      for ( size_t j = 0; j < coef_rank; ++j )
      {
        size_t dim    = coefficients[i]->dim( j );
        size_t fe_dim = fe->value_dimension( j );
        if ( dim != fe_dim )
        {
          error(
            "Invalid value dimension %d of Coefficient '%s' with index %d:\n"
            "got %d but expecting %d.\n"
            "You may need to provide the dimension of a user defined "
            "Coefficient.",
            j,
            this->coefficient_name( i ).c_str(),
            i,
            dim,
            fe_dim );
        }
      }
    }
    delete fe;
  }

  // Check that the cell dimension matches the mesh dimension
  if ( this->form().rank() + this->form().num_coefficients() > 0 )
  {
    ufc::finite_element * element = this->form().create_finite_element( 0 );
    dolfin_assert( element );
    CellType::Type celltype = mesh().type().cellType();
    ufc::shape     shape    = element->cell_shape();

    if ( celltype == CellType::interval && shape != ufc::shape::interval )
    {
      error( "Mesh cell type (intervals) does not match cell type of form." );
    }
    if ( celltype == CellType::triangle && shape != ufc::shape::triangle )
    {
      error( "Mesh cell type (triangles) does not match cell type of form." );
    }
    if ( celltype == CellType::tetrahedron && shape != ufc::shape::tetrahedron )
    {
      error( "Mesh cell type (tetrahedra) does not match cell type of form." );
    }
    delete element;
  }
  return true;
}

//-----------------------------------------------------------------------------
auto Form::is_valid_index( size_t i ) const -> bool
{
  // Check argument
  size_t const num_arguments = form().rank() + form().num_coefficients();
  if ( i >= num_arguments )
  {
    error( "Illegal function index %d. Form only has %d arguments.",
           i,
           num_arguments );
  }
  return true;
}

//----------------------------------------------------------------------------
void Form::assemble( GenericTensor & T, bool reset_tensor )
{
  Assembler::assemble( T, *this, reset_tensor );
}

//-----------------------------------------------------------------------------
void Form::init( std::vector< Coefficient * > & coefficients,
                 CoefficientMap &               map )
{
  coefficients.clear();
  for ( size_t i = 0; i < this->form().num_coefficients(); ++i )
  {
    std::string   name = this->coefficient_name( i );
    Coefficient * c    = map[name];
    if ( c != nullptr )
    {
      coefficients.push_back( map[name] );
    }
    else
    {
      error( "Missing coefficient named '%s' in CoefficientMap.",
             name.c_str() );
    }
  }
  Form::init( coefficients );
}

//----------------------------------------------------------------------------
void Form::init( std::vector< Coefficient * > & coefficients )
{
  if ( coefficients.size() != this->form().num_coefficients() )
  {
    error( "Form : invalid number of coefficients" );
  }
  for ( size_t i = 0; i < this->form().num_coefficients(); ++i )
  {
    Function * fptr = dynamic_cast< Function * >( this->coefficients()[i] );
    if ( fptr != nullptr && fptr->empty() )
    {
      fptr->init( *this, this->form().rank() + i );
      dolfin_assert( !fptr->empty() );
    }
  }
  Form::check( coefficients );
}

//----------------------------------------------------------------------------

}
