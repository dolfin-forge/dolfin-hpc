#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_ElementLibrary )
{
  ufc::finite_element * element = NULL;
  ufl::ElementList const& L = ElementLibrary::elements();
  for (ufl::FiniteElementSpace const * e = L.first(); L.valid(); e = L.next())
  {
    element = ElementLibrary::create_finite_element(e->repr());
    delete element;
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
