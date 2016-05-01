#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/Test.h>
#include <dolfin/log/log.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

#include <check.h>

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
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
