#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "Constant/Constant.h"
#include "Expression/Expression.h"
#include "Field/Field.h"
#include "Function/Function.h"
#include "Operators/Operators.h"
#include "Real/Real.h"
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

  DOLFIN_TCASE_CREATE("Real");
  DOLFIN_TCASE_ADD(test_Real);

  DOLFIN_TCASE_CREATE("Operators");
  DOLFIN_TCASE_ADD(test_Operators);

  DOLFIN_TCASE_CREATE("UFCFunction");
  DOLFIN_TCASE_ADD(test_UFCFunction);

  DOLFIN_TCASE_CREATE("Function");
  DOLFIN_TCASE_ADD(test_Function);

  DOLFIN_TCASE_CREATE("Field");
  DOLFIN_TCASE_ADD(test_Field);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/function", suite)
//-----------------------------------------------------------------------------

#endif
