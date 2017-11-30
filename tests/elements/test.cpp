#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "ElementLibrary/ElementLibrary.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "elements")
{
  DOLFIN_TCASE_CREATE("ElementLibrary");
  DOLFIN_TCASE_ADD(test_ElementLibrary);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/elements", suite)
//-----------------------------------------------------------------------------

#endif
