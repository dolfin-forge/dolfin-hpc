#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "LogStream/LogStream.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "log")
{
  DOLFIN_TCASE_CREATE("LogStream");
  DOLFIN_TCASE_ADD(test_LogStream);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/log", suite)
//-----------------------------------------------------------------------------

#endif
