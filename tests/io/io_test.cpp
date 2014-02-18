#include <dolfin/config/dolfin_config.h>

#include <iostream>

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
START_TEST( test_init_offset_file )
{
  int init_failed = 0;


  fail_unless( init_failed == 0 );
}END_TEST

Suite *io_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("FEM");
  tc = tcase_create("io");

  tcase_add_test(tc, test_init_offset_file);

  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);

  return s;
}

int main(void)
{
  int number_failed;
  Suite* s = io_suite();
  SRunner* sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;
}

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/io tests ***\n");
  return 0;
}

#endif
