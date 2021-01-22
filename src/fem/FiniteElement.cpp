// Copyright (C) 2013 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/FiniteElement.h>

#include <dolfin/fem/Elements.h>
#include <dolfin/fem/Form.h>

#include <algorithm>
#include <iomanip>

namespace dolfin
{

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( ufc::finite_element const & element,
                              bool const                  owner )
  : ufc_finite_element_( ( owner ? &element : element.create() ) )
  , sub_value_dims_( nullptr )
{
  Initialize();
}

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( ufc::finite_element const & element,
                              size_t const                i )
  : ufc_finite_element_( element.create_sub_element( i ) )
  , sub_value_dims_( nullptr )
{
  Initialize();
}

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( ufc::finite_element const & element,
                              Array< size_t > const &     sub_system )
  : ufc_finite_element_(
    FiniteElement::create_sub_element( element, sub_system ) )
  , sub_value_dims_( nullptr )
{
  Initialize();
}

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( CellType const &, Form & form, size_t const i )
  : ufc_finite_element_( nullptr )
  , sub_value_dims_( nullptr )
{
  // Check argument
  size_t const num_arguments = form.rank() + form.num_coefficients();
  if ( i >= num_arguments )
  {
    error( "Illegal function index %d. Form only has %d arguments.",
           i,
           num_arguments );
  }

  // Create finite element
  ufc_finite_element_ = form.create_finite_element( i );

  Initialize();
}

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( FiniteElement const & other )
  : ufc_finite_element_( other.create() )
  , sub_value_dims_( nullptr )
{
  Initialize();
}

//-----------------------------------------------------------------------------

FiniteElement::~FiniteElement()
{
  while ( !flattened_.empty() )
  {
    delete flattened_.back();
    flattened_.pop_back();
  }
  delete[] sub_value_dims_;
  delete[] sub_value_offs_;
  delete ufc_finite_element_;
  ufc_finite_element_ = nullptr;
}

//-----------------------------------------------------------------------------

void FiniteElement::Initialize()
{
  dolfin_assert( ufc_finite_element_ );

  // Add sub value dimensions for mixed elements, packed by axis
  size_t const max_dim = ufc_finite_element_->value_rank() + 1;
  sub_value_dims_      = new Array< size_t >[max_dim];
  sub_value_offs_      = new Array< size_t >[max_dim];
  size_t nb_subs       = this->num_sub_elements();
  if ( nb_subs > 0 )
  {
    size_t * off = new size_t[max_dim];
    std::fill_n( off, max_dim, 0 );
    for ( size_t e = 0; e < nb_subs; ++e )
    {
      ufc::finite_element * sub_fe =
        ufc_finite_element_->create_sub_element( e );
      for ( size_t a = 0; a < max_dim; ++a )
      {
        sub_value_dims_[a].push_back( sub_fe->value_dimension( a ) );
        sub_value_offs_[a].push_back( off[a] );
        off[a] += sub_fe->value_dimension( a );
      }
      delete sub_fe;
    }
    delete[] off;
  }
  else
  {
    for ( size_t a = 0; a < max_dim; ++a )
    {
      sub_value_dims_[a].push_back( value_dimension( a ) );
      sub_value_offs_[a].push_back( 0 );
    }
  }
}

//-----------------------------------------------------------------------------

auto FiniteElement::create_sub_element(
  const ufc::finite_element & finite_element,
  Array< size_t > const &     sub_system ) -> ufc::finite_element *
{
  // If the subsystem is empty return self
  if ( sub_system.size() == 0 )
  {
    // error("Unable to extract sub system (no sub system specified).");
    return finite_element.create();
  }

  // Check if there are any sub systems
  if ( finite_element.num_sub_elements() == 0 )
  {
    error( "Unable to extract sub system (there are no sub systems)." );
  }

  // Check the number of available sub systems
  if ( sub_system[0] >= finite_element.num_sub_elements() )
  {
    error( "Unable to extract sub system %d (only %d sub systems defined).",
           sub_system[0],
           finite_element.num_sub_elements() );
  }

  // Create sub system
  ufc::finite_element * sub_element =
    finite_element.create_sub_element( sub_system[0] );

  // Return sub system if sub sub system should not be extracted
  if ( sub_system.size() == 1 )
  {
    return sub_element;
  }

  // Otherwise, recursively extract the sub sub system
  Array< size_t > sub_sub_system;
  for ( size_t i = 1; i < sub_system.size(); i++ )
  {
    sub_sub_system.push_back( sub_system[i] );
  }
  ufc::finite_element * sub_sub_element =
    create_sub_element( *sub_element, sub_sub_system );
  delete sub_element;

  return sub_sub_element;
}

//-----------------------------------------------------------------------------

void FiniteElement::flatten( ufc::finite_element const *            element,
                             Array< ufc::finite_element const * > & stack,
                             size_t                                 maxlevel )
{
  // Single root element or max level is set to zero, return immediately
  if ( element->num_sub_elements() == 0 || maxlevel == 0 )
  {
    stack.push_back( element->create() );
    return;
  }
  // Go one level down
  for ( size_t s = 0; s < element->num_sub_elements(); ++s )
  {
    ufc::finite_element const * sub = element->create_sub_element( s );
    if ( sub->num_sub_elements() == 0 )
    {
      // Leaf element
      stack.push_back( sub );
    }
    else
    {
      // Branch
      FiniteElement::flatten( sub, stack, maxlevel - 1 );
    }
  }
}

//-----------------------------------------------------------------------------

void FiniteElement::flatten( ufc::finite_element const *            element,
                             Array< ufc::finite_element const * > & stack )
{
  // Single root element or max level is set to zero, return immediately
  if ( element->num_sub_elements() == 0 )
  {
    stack.push_back( element->create() );
    return;
  }
  // Go one level down
  for ( size_t s = 0; s < element->num_sub_elements(); ++s )
  {
    ufc::finite_element const * sub = element->create_sub_element( s );
    if ( sub->num_sub_elements() == 0 )
    {
      // Leaf element
      stack.push_back( sub );
    }
    else
    {
      // Branch
      FiniteElement::flatten( sub, stack );
    }
  }
}

//-----------------------------------------------------------------------------

auto FiniteElement::is_vectorizable() const -> bool
{
  bool                                         ret = true;
  Array< ufc::finite_element const * > const & flt = this->flatten();
  for ( size_t s = 1; s < flt.size(); ++s )
  {
    if ( std::strcmp( flt[0]->signature(), flt[s]->signature() ) != 0 )
    {
      ret = false;
      break;
    }
  }
  return ret;
}

//-----------------------------------------------------------------------------

void FiniteElement::disp() const
{
  std::string shape;
  switch ( ufc_finite_element_->cell_shape() )
  {
    case ufc::shape::interval:
      shape = "interval";
      break;
    case ufc::shape::triangle:
      shape = "triangle";
      break;
    case ufc::shape::quadrilateral:
      shape = "quadrilateral";
      break;
    case ufc::shape::tetrahedron:
      shape = "tetrahedron";
      break;
    case ufc::shape::hexahedron:
      shape = "hexahedron";
      break;
    case ufc::shape::vertex:
      shape = "vertex";
      break;
    default:
      shape = "unknown shape";
      break;
  }

  section( "FiniteElement" );
  begin( "ufc::finite_element info" );
  prm( "Signature", ufc_finite_element_->signature() );
  prm( "Cell shape", shape.c_str() );
  prm( "Topological dimension", ufc_finite_element_->topological_dimension() );
  prm( "Geometric dimension", ufc_finite_element_->geometric_dimension() );
  prm( "Space dimension", ufc_finite_element_->space_dimension() );
  prm( "Value rank", ufc_finite_element_->value_rank() );
  prm( "Value dimension", ufc_finite_element_->value_dimension( 0 ) );
  prm( "Nb of sub elements", ufc_finite_element_->num_sub_elements() );
  end();
  end();
}

//-----------------------------------------------------------------------------

} // namespace dolfin
