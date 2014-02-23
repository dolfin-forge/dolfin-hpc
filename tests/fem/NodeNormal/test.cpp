#include <dolfin/config/dolfin_config.h>

#include <dolfin/fem/NodeNormal.h>
#include <dolfin/mesh/Mesh.h>

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

//-----------------------------------------------------------------------------
START_TEST( test_init )
  {
    int init_failed = 0;

    Mesh mesh("heart.bin");
    NodeNormal nn(mesh, VertexNormal::none);

    ufl::VectorElement space(ufl::Family::CG, mesh.type(), 1,
                             mesh.topology().dim());
    nn.init(mesh, space.repr());

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------

Suite * test_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("fem");
  tc = tcase_create("fem_NodeNormal");

  tcase_add_test(tc, test_init);

  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);

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
  fprintf(stderr, "*** Check is required for dolfin/fem tests ***\n");
  return 0;
}

#endif
