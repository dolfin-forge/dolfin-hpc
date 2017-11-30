#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "Array/Array.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "common")
{
  DOLFIN_TCASE_CREATE("Array");
  DOLFIN_TCASE_ADD(test_Array);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/common", suite)
//-----------------------------------------------------------------------------

#endif
