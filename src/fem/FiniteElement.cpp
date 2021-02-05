// Copyright (C) 2013 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/FiniteElement.h>

#include <dolfin/common/assert.h>
#include <dolfin/fem/Elements.h>
#include <dolfin/log/log.h>

#include <algorithm>
#include <iomanip>

namespace dolfin
{

//-----------------------------------------------------------------------------

namespace helper
{

inline auto init_vdims( ufc::finite_element const & ufc_element )
  -> std::vector< size_t >
{
  std::vector< size_t > vdims( ufc_element.value_size() );

  for ( size_t i = 0; i < vdims.size(); ++i )
    vdims[i] = ufc_element.value_dimension( i );

  return vdims;
}

//-----------------------------------------------------------------------------

inline auto init_rvdims( ufc::finite_element const & ufc_element )
  -> std::vector< size_t >
{
  std::vector< size_t > rvdims( ufc_element.reference_value_size() );

  for ( size_t i = 0; i < rvdims.size(); ++i )
    rvdims[i] = ufc_element.reference_value_dimension( i );

  return rvdims;
}

//-----------------------------------------------------------------------------

inline auto init_flattened( ufc::finite_element const & ufc_element )
  -> std::vector< ufc::finite_element const * >
{
  std::vector< ufc::finite_element const * > flt;
  FiniteElement::flatten( &ufc_element, flt );

  return flt;
}

//-----------------------------------------------------------------------------

inline auto init_sv_dims( ufc::finite_element const & ufc_element )
  -> std::vector< std::vector< size_t > >
{
  std::vector< std::vector< size_t > > svd( ufc_element.value_rank() + 1 );

  // Add sub value dimensions for mixed elements, packed by axis
  size_t const nb_subs = ufc_element.num_sub_elements();
  if ( nb_subs > 0 )
  {
    for ( size_t e = 0; e < nb_subs; ++e )
    {
      ufc::finite_element * sub_fe = ufc_element.create_sub_element( e );
      for ( size_t a = 0; a < svd.size(); ++a )
      {
        svd[a].push_back( sub_fe->value_dimension( a ) );
      }
      delete sub_fe;
    }
  }
  else
  {
    for ( size_t a = 0; a < svd.size(); ++a )
    {
      svd[a].push_back( ufc_element.value_dimension( a ) );
    }
  }

  return svd;
}

//-----------------------------------------------------------------------------

inline auto init_sv_offs( ufc::finite_element const & ufc_element )
  -> std::vector< std::vector< size_t > >
{
  std::vector< std::vector< size_t > > svo( ufc_element.value_rank() + 1 );

  // Add sub value dimensions for mixed elements, packed by axis
  size_t const nb_subs = ufc_element.num_sub_elements();
  if ( nb_subs > 0 )
  {
    std::vector< size_t > off( svo.size(), 0 );
    for ( size_t e = 0; e < nb_subs; ++e )
    {
      ufc::finite_element * sub_fe = ufc_element.create_sub_element( e );
      for ( size_t a = 0; a < svo.size(); ++a )
      {
        svo[a].push_back( off[a] );
        off[a] += sub_fe->value_dimension( a );
      }
      delete sub_fe;
    }
  }
  else
  {
    for ( size_t a = 0; a < svo.size(); ++a )
    {
      svo[a].push_back( 0 );
    }
  }

  return svo;
}

} // namespace helper

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( ufc::finite_element const & ufc_element )
  : element( ufc_element.create() )
  , signature( element->signature() )
  , family( element->family() )
  , shape( element->cell_shape() )
  , degree( element->degree() )
  , tdim( element->topological_dimension() )
  , gdim( element->geometric_dimension() )
  , space_dim( element->space_dimension() )
  , num_sub_elements( element->num_sub_elements() )
  , value_rank( element->value_rank() )
  , value_size( element->value_size() )
  , value_dims( helper::init_vdims( *element ) )
  , ref_value_rank( element->reference_value_rank() )
  , ref_value_size( element->reference_value_size() )
  , ref_value_dims( helper::init_vdims( *element ) )
  , sub_value_dims_( helper::init_sv_dims( *element ) )
  , sub_value_offs_( helper::init_sv_offs( *element ) )
  , flattened_( helper::init_flattened( *element ) )
{
}

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( ufc::finite_element const & ufc_element,
                              size_t const                i )
  : element( ufc_element.create_sub_element( i ) )
  , signature( element->signature() )
  , family( element->family() )
  , shape( element->cell_shape() )
  , degree( element->degree() )
  , tdim( element->topological_dimension() )
  , gdim( element->geometric_dimension() )
  , space_dim( element->space_dimension() )
  , num_sub_elements( element->num_sub_elements() )
  , value_rank( element->value_rank() )
  , value_size( element->value_size() )
  , value_dims( helper::init_vdims( *element ) )
  , ref_value_rank( element->reference_value_rank() )
  , ref_value_size( element->reference_value_size() )
  , ref_value_dims( helper::init_vdims( *element ) )
  , sub_value_dims_( helper::init_sv_dims( *element ) )
  , sub_value_offs_( helper::init_sv_offs( *element ) )
  , flattened_( helper::init_flattened( *element ) )
{
}

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( ufc::finite_element const &   ufc_element,
                              std::vector< size_t > const & sub_system )
  : element( FiniteElement::create_sub_element( ufc_element, sub_system ) )
  , signature( element->signature() )
  , family( element->family() )
  , shape( element->cell_shape() )
  , degree( element->degree() )
  , tdim( element->topological_dimension() )
  , gdim( element->geometric_dimension() )
  , space_dim( element->space_dimension() )
  , num_sub_elements( element->num_sub_elements() )
  , value_rank( element->value_rank() )
  , value_size( element->value_size() )
  , value_dims( helper::init_vdims( *element ) )
  , ref_value_rank( element->reference_value_rank() )
  , ref_value_size( element->reference_value_size() )
  , ref_value_dims( helper::init_vdims( *element ) )
  , sub_value_dims_( helper::init_sv_dims( *element ) )
  , sub_value_offs_( helper::init_sv_offs( *element ) )
  , flattened_( helper::init_flattened( *element ) )
{
}

//-----------------------------------------------------------------------------

FiniteElement::FiniteElement( FiniteElement const & other )
  : element( other.ufc().create() )
  , signature( element->signature() )
  , family( element->family() )
  , shape( element->cell_shape() )
  , degree( element->degree() )
  , tdim( element->topological_dimension() )
  , gdim( element->geometric_dimension() )
  , space_dim( element->space_dimension() )
  , num_sub_elements( element->num_sub_elements() )
  , value_rank( element->value_rank() )
  , value_size( element->value_size() )
  , value_dims( helper::init_vdims( *element ) )
  , ref_value_rank( element->reference_value_rank() )
  , ref_value_size( element->reference_value_size() )
  , ref_value_dims( helper::init_vdims( *element ) )
  , sub_value_dims_( helper::init_sv_dims( *element ) )
  , sub_value_offs_( helper::init_sv_offs( *element ) )
  , flattened_( helper::init_flattened( *element ) )
{
}

//-----------------------------------------------------------------------------

FiniteElement::~FiniteElement()
{
  delete element;
  destruct( flattened_ );
}

//-----------------------------------------------------------------------------

auto FiniteElement::create_sub_element(
  const ufc::finite_element &   finite_element,
  std::vector< size_t > const & sub_system ) -> ufc::finite_element *
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
  std::vector< size_t > sub_sub_system;
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

void FiniteElement::flatten(
  ufc::finite_element const *                  element,
  std::vector< ufc::finite_element const * > & stack )
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
  bool                                               ret = true;
  std::vector< ufc::finite_element const * > const & flt = this->flatten();
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
  std::string shape_;
  switch ( element->cell_shape() )
  {
    case ufc::shape::interval:
      shape_ = "interval";
      break;
    case ufc::shape::triangle:
      shape_ = "triangle";
      break;
    case ufc::shape::quadrilateral:
      shape_ = "quadrilateral";
      break;
    case ufc::shape::tetrahedron:
      shape_ = "tetrahedron";
      break;
    case ufc::shape::hexahedron:
      shape_ = "hexahedron";
      break;
    case ufc::shape::vertex:
      shape_ = "vertex";
      break;
    default:
      shape_ = "unknown shape";
      break;
  }

  section( "FiniteElement" );
  begin( "ufc::finite_element info" );
  prm( "Signature",             signature );
  prm( "Family",                family );
  prm( "Cell shape",            shape_ );
  prm( "Degree",                degree );
  prm( "Topological dimension", tdim );
  prm( "Geometric dimension",   gdim );
  prm( "Space dimension",       space_dim );
  prm( "Value rank",            value_rank );
  prm( "Value size",            value_size );
  prm( "Value dimension",       to_string( value_dims ) );
  prm( "Ref. Value rank",       ref_value_rank );
  prm( "Ref. Value size",       ref_value_size );
  prm( "Ref. Value dimension",  to_string( ref_value_dims ) );
  prm( "# of sub elements",     num_sub_elements );
  end();
  end();
}

//-----------------------------------------------------------------------------

} // namespace dolfin
