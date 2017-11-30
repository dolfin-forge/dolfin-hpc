#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_ElementLibrary )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_ElementLibrary");
  {
    ufl::ElementList const& list = ElementLibrary::elements();
    ufc::finite_element * fe = NULL;
    for (ufl::FiniteElementSpace const * it = list.first(); list.valid();
        it = list.next())
    {
      it->display();
      //
      fe = ElementLibrary::create_finite_element(it->repr());
      delete fe;
    }
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
