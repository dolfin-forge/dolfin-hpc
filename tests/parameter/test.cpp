#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "ParameterValue/ParameterValue.h"
#include "Parameter/Parameter.h"
#include "Parametrized/Parametrized.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "main")
{
  DOLFIN_TCASE_CREATE("ParameterValue");
  DOLFIN_TCASE_ADD(test_ParameterValue);

  DOLFIN_TCASE_CREATE("Parameter");
  DOLFIN_TCASE_ADD(test_Parameter);

  DOLFIN_TCASE_CREATE("Parametrized");
  DOLFIN_TCASE_ADD(test_Parametrized);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/parameter", suite)
//-----------------------------------------------------------------------------

#endif
