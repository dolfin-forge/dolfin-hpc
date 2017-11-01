#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include "DiscreteSpaces/DiscreteSpaces.h"
#include "FiniteElement/FiniteElement.h"
#include "DofNumbering/DofNumbering.h"
#include "DofMap/DofMap.h"
#include "FiniteElementSpace/FiniteElementSpace.h"
#include "Form/Form.h"
#include "BoundaryNormal/BoundaryNormal.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "fem")
{
  DOLFIN_TCASE_CREATE("DiscreteSpaces");
  DOLFIN_TCASE_ADD(test_DiscreteSpaces);

  DOLFIN_TCASE_CREATE("FiniteElement");
  DOLFIN_TCASE_ADD(test_FiniteElement);

  DOLFIN_TCASE_CREATE("DofNumbering");
  DOLFIN_TCASE_ADD(test_DofNumbering);

  DOLFIN_TCASE_CREATE("DofMap");
  DOLFIN_TCASE_ADD(test_DofMap);

  DOLFIN_TCASE_CREATE("FiniteElementSpace");
  DOLFIN_TCASE_ADD(test_FiniteElementSpace);

  DOLFIN_TCASE_CREATE("Form");
  DOLFIN_TCASE_ADD(test_Functional);
  DOLFIN_TCASE_ADD(test_LinearForm);
  DOLFIN_TCASE_ADD(test_BilinearForm);

  DOLFIN_TCASE_CREATE("BoundaryNormal");
  DOLFIN_TCASE_ADD(test_NodeNormal);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/fem", suite)
//-----------------------------------------------------------------------------

#endif
