#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include "alarm/alarm.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "main")
{
  DOLFIN_TCASE_CREATE("alarm");
  DOLFIN_TCASE_ADD(test_alarm);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/main", suite)
//-----------------------------------------------------------------------------

#endif
