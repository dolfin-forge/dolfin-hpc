#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "Analytic/Analytic.h"
#include "Constant/Constant.h"
#include "Expression/Expression.h"
#include "Field/Field.h"
#include "Function/Function.h"
#include "Operators/Operators.h"
#include "Real/Real.h"
#include "SpaceTimeFunction/SpaceTimeFunction.h"
#include "Value/Value.h"
#include "UFCFunction/UFCFunction.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "function")
{
  DOLFIN_TCASE_CREATE("Constant");
  DOLFIN_TCASE_ADD(test_Constant);

  DOLFIN_TCASE_CREATE("Expression");
  DOLFIN_TCASE_ADD(test_Expression);

  DOLFIN_TCASE_CREATE("Value");
  DOLFIN_TCASE_ADD(test_Value);

  DOLFIN_TCASE_CREATE("Analytic");
  DOLFIN_TCASE_ADD(test_Analytic);

  DOLFIN_TCASE_CREATE("Real");
  DOLFIN_TCASE_ADD(test_Real1);
  DOLFIN_TCASE_ADD(test_Real2);
  DOLFIN_TCASE_ADD(test_Realt);

  DOLFIN_TCASE_CREATE("Operators");
  DOLFIN_TCASE_ADD(test_Operators);

  DOLFIN_TCASE_CREATE("UFCFunction");
  DOLFIN_TCASE_ADD(test_UFCFunction);

  DOLFIN_TCASE_CREATE("Function");
  DOLFIN_TCASE_ADD(test_Function);

  DOLFIN_TCASE_CREATE("Field");
  DOLFIN_TCASE_ADD(test_Field);

  DOLFIN_TCASE_CREATE("SpaceTimeFunction");
  DOLFIN_TCASE_ADD(test_SpaceTimeFunction);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/function", suite)
//-----------------------------------------------------------------------------

#endif
