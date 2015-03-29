#include <dolfin/config/dolfin_config.h>

#include <iostream>
#include <iomanip>

#ifdef HAVE_CHECK

#include <check.h>

#include "Vector/Vector.h"
#include "Matrix/Matrix.h"


int argc;
char **argv;

int main(void)
{

  Suite* s;
  SRunner *sr;
  int number_failed;

  s = test_suite_vec();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  s = test_suite_mat();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;

}

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/la tests ***\n");
  return 0;
}

#endif
