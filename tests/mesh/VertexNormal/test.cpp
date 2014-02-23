#include <dolfin/config/dolfin_config.h>

#include <dolfin/mesh/VertexNormal.h>

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

void test_vertex_normal(std::string filename, VertexNormal::Type type)
{
  Mesh mesh(filename);
  VertexNormal vn(mesh, type);

//  uint const gdim = mesh.geometry().dim();
//  for(uint i = 0 ; i < gdim ; ++i)
//  {
//    std::stringstream ss;
//    ss << "E" << i;
//    size_t pos = filename.find('.');
//    File f(filename.substr(0,pos)+".pvd");
//    f << vn.basis()[d][i];
//  }
}

//-----------------------------------------------------------------------------
START_TEST( test_init_weight_none )
  {
    int init_failed = 0;

    std::string const relpath = "../../../data/meshes/";
    test_vertex_normal(relpath+"square16.bin", VertexNormal::none);
    test_vertex_normal(relpath+"cylinder.xml.gz", VertexNormal::none);
    test_vertex_normal(relpath+"aneurysm.xml.gz", VertexNormal::none);

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------

Suite * test_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("data");
  tc = tcase_create("data");

  tcase_add_test(tc, test_init_weight_none );

  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

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
  fprintf(stderr, "*** Check is required for dolfin/data tests ***\n");
  return 0;
}

#endif
