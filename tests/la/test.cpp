#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "Vector/Vector.h"
#include "Matrix/Matrix.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "la")
{
  DOLFIN_TCASE_CREATE("Vector");
  DOLFIN_TCASE_ADD(test_init_vec);
  DOLFIN_TCASE_ADD(test_get_set_vec);
  DOLFIN_TCASE_ADD(test_add_vec);
  DOLFIN_TCASE_ADD(test_max_vec);
  DOLFIN_TCASE_ADD(test_min_vec);
  DOLFIN_TCASE_ADD(test_op_assign_vec);
  DOLFIN_TCASE_ADD(test_op_mul_vec);
  DOLFIN_TCASE_ADD(test_op_scale_vec);
  DOLFIN_TCASE_ADD(test_zero_vec);
  DOLFIN_TCASE_ADD(test_trilinos);

  DOLFIN_TCASE_CREATE( "Matrix" );
  DOLFIN_TCASE_ADD(test_init_mat);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/la", suite)
//-----------------------------------------------------------------------------

#endif
