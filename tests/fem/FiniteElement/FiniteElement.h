#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/FiniteElement.h>

#include "../../elements/element_library.inc"

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_FiniteElement )
{
  using namespace ElementLibrary;

  for ( char const * element_str : elements )
  {
    ufc::finite_element * ufc_element = create_finite_element( element_str );
    FiniteElement element( *ufc_element, false );
    ck_assert( strcmp( ufc_element->signature(), element.signature() ) == 0 );
    ck_assert( ufc_element->cell_shape() == element.cell_shape() );
    ck_assert( ufc_element->topological_dimension() == element.topological_dimension() );
    ck_assert( ufc_element->geometric_dimension() == element.geometric_dimension() );
    ck_assert( ufc_element->space_dimension() == element.space_dimension() );
    ck_assert( ufc_element->value_rank() == element.value_rank() );
    ck_assert( ufc_element->value_dimension(0) == element.value_dimension(0) );
    ck_assert( ufc_element->value_size() == element.value_size() );
    ck_assert( ufc_element->reference_value_rank() == element.reference_value_rank() );
    ck_assert( ufc_element->reference_value_dimension(0) == element.reference_value_dimension( 0 ) );
    ck_assert( ufc_element->reference_value_size() == element.reference_value_size() );
    ck_assert( ufc_element->degree() == element.degree() );
    ck_assert( strcmp( ufc_element->family(), element.family() ) == 0 );
    element.is_vectorizable();
    for ( size_t i = 0; i < element.num_sub_elements(); ++i )
    {
      ufc::finite_element * sub = element.create_sub_element( i );
      delete sub;
    }
    delete ufc_element;
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
