#include <dolfin/config/dolfin_config.h>

#include <dolfin/mesh/DynamicSubDomain.h>
#include <dolfin/mesh/GeometricSubDomain.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/UnitSquare.h>

#include <iostream>
#include <iomanip>

using namespace dolfin;

#ifdef HAVE_CHECK

#include <check.h>

int argc;
char * argv;

void setup()
{
}

void teardown()
{
}

bool test_create(SubDomain& sd)
{
  uint const dim = sd.mesh().topology().dim();
  for (uint i = 0; i < dim; ++i)
  {
    MeshFunction<uint>& marker = sd.marker(i);
  }
  return false;
}

//-----------------------------------------------------------------------------
START_TEST( test_dynamic_create )
{
  int init_failed = 0;

  uint const N = 16;
  UnitSquare mesh(N, N);
  DynamicSubDomain sd(mesh);
  init_failed = test_create(sd);

  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

Suite * test_suite()
{
  Suite * s = suite_create("mesh/SubDomain");

  // Dynamic
  TCase * tc_dynamic = tcase_create("mesh/SubDomain/DynamicSubDomain");
  tcase_add_test(tc_dynamic, test_dynamic_create);
  suite_add_tcase(s, tc_dynamic);

  // Geometric
//  TCase * tc_geometric = tcase_create("mesh/SubDomain/GeometricSubDomain");
//  tcase_add_test(tc_geometric, test_geometric_create);
//  suite_add_tcase(s, tc_geometric);

  return s;
}

int main(void)
{
  int number_failed;
  Suite* s = test_suite();
  SRunner* sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;
}

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/mesh/SubDomain tests ***\n");
  return 0;
}

#endif
