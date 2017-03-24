#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/Test.h>
#include <dolfin/log/log.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_FiniteElement )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_FiniteElement");
  {
    ufl::ElementList const& list = ElementLibrary::elements();
    for (ufl::FiniteElementSpace const * it = list.first(); list.valid();
        it = list.next())
    {
      begin("%s", it->str().c_str());
      //---

      FiniteElement e(*ElementLibrary::create_finite_element(it->repr()), true);

      //---
      end();
    }
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
