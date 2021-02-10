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
{
}

//-----------------------------------------------------------------------------

Form::~Form()
{
  for ( DofMap * dofmap : dofmaps_ )
    DofMap::release( *dofmap );

  destruct( cell_integrals_ );
  destruct( exterior_facet_integrals_ );
  destruct( interior_facet_integrals_ );
  destruct( vertex_integrals_ );
  destruct( custom_integrals_ );
  destruct( cutcell_integrals_ );
  destruct( interface_integrals_ );
  destruct( overlap_integrals_ );
}

//-----------------------------------------------------------------------------

auto Form::check( std::vector< Coefficient * > const & coefficients ) const
  -> bool
{
  // Check that we get the correct number of coefficients
  if ( coefficients.size() != this->form().num_coefficients() )
  {
    error( "Incorrect number of coefficients: %d given but %d required.",
           coefficients.size(), this->form().num_coefficients() );
  }

  // Check that all coefficients have valid value dimensions
  for ( size_t i = 0; i < coefficients.size(); ++i )
  {
    message( 1, "Form: Checking coefficient %d: %s",
             i, this->coefficient_name( i ).c_str() );
    if ( coefficients[i] == nullptr )
    {
      error( "Got nullptr pointer as coefficient %d labeled as '%s'.",
             i, this->coefficient_name( i ).c_str() );
    }

    FiniteElement const & fe = this->elements()[ i + rank_ ];
    Function * fptr = dynamic_cast< Function * >( this->coefficients()[i] );
    if ( fptr != nullptr )
    {
      if ( fptr->empty() )
      {
        error( "Coefficient %i is empty", i );
      }
      else if ( fptr->space().element().signature != fe.signature )
      {
        error( "Mismatch of discrete space for Coefficient %i", i );
      }
    }
    else
    {
      message( 1, "Form: Coefficient rank: expected  = %d, provided = %d, ",
               fe.value_rank, coefficients[i]->rank() );
      if ( fe.value_rank != coefficients[i]->rank() )
      {
        error( "Invalid value rank of Coefficient '%s' with index %d:\n"
               "Got %d but expecting %d.\n"
               "You may need to provide the rank of a user defined Coefficient.",
               this->coefficient_name( i ).c_str(), i,
               coefficients[i]->rank(), fe.value_rank );
      }

      for ( size_t j = 0; j < coefficients[i]->rank(); ++j )
      {
        if ( coefficients[i]->dim( j ) != fe.value_dims[j] )
        {
          error( "Invalid value dimension %d of Coefficient '%s' with index %d:"
                 "\ngot %d but expecting %d.\nYou may need to provide the "
                 "dimension of a user defined Coefficient.",
                 j, this->coefficient_name( i ).c_str(), i,
                 coefficients[i]->dim( j ), fe.value_dims[j] );
        }
      }
    }
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
           i, num_arguments );
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
  signature_        = std::string( this->form().signature() );
  rank_             = this->form().rank();
  num_coefficients_ = this->form().num_coefficients();

  // init FiniteElementSpaces
  if ( spaces_.empty() )
  {
    for ( size_t i = 0; i < rank_ + num_coefficients_ ; ++i )
    {
      // create FiniteElement
      ufc::finite_element * fe = this->form().create_finite_element( i );
      dolfin_assert( fe != nullptr );
      elements_.emplace_back( *fe );
      delete fe;

      // create DofMap
      ufc::dofmap * dof = this->form().create_dofmap( i );
      dolfin_assert( dof != nullptr );
      dofmaps_.emplace_back( &(DofMap::acquire( mesh_, *dof, false )) );
      delete dof;

      // create FiniteElementSpaces
      spaces_.emplace_back( elements_[i], *dofmaps_[i] );
    }

    // initialize UFCCache
    cache_.init( *this );
  }

  // Create cell integrals
  if ( this->form().max_cell_subdomain_id() > 0 )
  {
    for ( size_t i = 0; i < this->form().max_cell_subdomain_id(); ++i )
      cell_integrals_.push_back( this->form().create_cell_integral( i ) );
  }
  else if ( this->form().has_cell_integrals() )
  {
    cell_integrals_.push_back( this->form().create_default_cell_integral() );
  }

  // Create exterior facet integrals
  if ( this->form().max_exterior_facet_subdomain_id() > 0 )
  {
    for ( size_t i = 0; i < this->form().max_exterior_facet_subdomain_id(); ++i )
      exterior_facet_integrals_.push_back( this->form().create_exterior_facet_integral( i ) );
  }
  else if ( this->form().has_exterior_facet_integrals() )
  {
    exterior_facet_integrals_.push_back( this->form().create_default_exterior_facet_integral() );
  }

  // Create interior facet integrals
  if ( this->form().max_interior_facet_subdomain_id() > 0 )
  {
    for ( size_t i = 0; i < this->form().max_interior_facet_subdomain_id(); ++i )
      interior_facet_integrals_.push_back( this->form().create_interior_facet_integral( i ) );
  }
  else if ( this->form().has_interior_facet_integrals() )
  {
    interior_facet_integrals_.push_back( this->form().create_default_interior_facet_integral() );
  }

  // Create vertex integrals
  if ( this->form().max_vertex_subdomain_id() > 0 )
  {
    for ( size_t i = 0; i < this->form().max_vertex_subdomain_id(); ++i )
      vertex_integrals_.push_back( this->form().create_vertex_integral( i ) );
  }
  else if ( this->form().has_vertex_integrals() )
  {
    vertex_integrals_.push_back( this->form().create_default_vertex_integral() );
  }

  // Create custom integrals
  if ( this->form().max_custom_subdomain_id() > 0 )
  {
    for ( size_t i = 0; i < this->form().max_custom_subdomain_id(); ++i )
      custom_integrals_.push_back( this->form().create_custom_integral( i ) );
  }
  else if ( this->form().has_custom_integrals() )
  {
    custom_integrals_.push_back( this->form().create_default_custom_integral() );
  }

  // Create cutcell integrals
  if ( this->form().max_cutcell_subdomain_id() > 0 )
  {
    for ( size_t i = 0; i < this->form().max_cutcell_subdomain_id(); ++i )
      cutcell_integrals_.push_back( this->form().create_cutcell_integral( i ) );
  }
  else if ( this->form().has_cutcell_integrals() )
  {
    cutcell_integrals_.push_back( this->form().create_default_cutcell_integral() );
  }

  // Create interface integrals
  if ( this->form().max_interface_subdomain_id() > 0 )
  {
    for ( size_t i = 0; i < this->form().max_interface_subdomain_id(); ++i )
      interface_integrals_.push_back( this->form().create_interface_integral( i ) );
  }
  else if ( this->form().has_interface_integrals() )
  {
    interface_integrals_.push_back( this->form().create_default_interface_integral() );
  }

  // Create overlap integrals
  if ( this->form().max_overlap_subdomain_id() > 0 )
  {
    for ( size_t i = 0; i < this->form().max_overlap_subdomain_id(); ++i )
      overlap_integrals_.push_back( this->form().create_overlap_integral( i ) );
  }
  else if ( this->form().has_overlap_integrals() )
  {
    overlap_integrals_.push_back( this->form().create_default_overlap_integral() );
  }

  // init coefficients
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
