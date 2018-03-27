#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_FiniteElement )
  {
    ufl::ElementList const& L = ElementLibrary::elements();
    for (ufl::FiniteElementSpace const * e = L.first(); L.valid(); e = L.next())
    {
      FiniteElement F(*ElementLibrary::create_finite_element(e->repr()), true);
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
