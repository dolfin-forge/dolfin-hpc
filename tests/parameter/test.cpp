#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "Parameter/Parameter.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "main")
{
  DOLFIN_TCASE_CREATE("parameter");
  DOLFIN_TCASE_ADD(test_parameter);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/parameter", suite)
//-----------------------------------------------------------------------------

#endif
