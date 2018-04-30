#include <dolfin_tests.h>

#ifdef HAVE_CHECK

// Include all test cases
#include "Cell/UFL_Cell.h"
#include "Domain/UFL_Domain.h"
#include "FiniteElement/UFL_FiniteElement.h"
#include "MixedElement/UFL_MixedElement.h"
#include "Space/UFL_Space.h"
#include "type/UFL_type.h"

DOLFIN_SUITE_BEGIN(suite, "ufl")
{
  DOLFIN_TCASE_CREATE("UFL_Cell");
  DOLFIN_TCASE_ADD(test_UFL_Cell);
  
  DOLFIN_TCASE_CREATE("UFL_Domain");
  DOLFIN_TCASE_ADD(test_UFL_Domain);

  DOLFIN_TCASE_CREATE("UFL_FiniteElement");
  DOLFIN_TCASE_ADD(test_UFL_FiniteElement);

  DOLFIN_TCASE_CREATE("UFL_MixedElement");
  DOLFIN_TCASE_ADD(test_UFL_MixedElement);

  DOLFIN_TCASE_CREATE("UFL_Space");
  DOLFIN_TCASE_ADD(test_UFL_Space);

  DOLFIN_TCASE_CREATE("UFL_type");
  DOLFIN_TCASE_ADD(test_UFL_type_int);
  DOLFIN_TCASE_ADD(test_UFL_type_real);
  DOLFIN_TCASE_ADD(test_UFL_type_string);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/ufl", suite)
//-----------------------------------------------------------------------------

#endif
