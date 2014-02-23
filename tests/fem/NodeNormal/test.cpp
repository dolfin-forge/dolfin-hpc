#include <dolfin/config/dolfin_config.h>

#include <dolfin/fem/NodeNormal.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/ufl/UFLVectorElement.h>

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

void test_node_normal(std::string filename, VertexNormal::Type type)
{
  Mesh mesh(filename);
  NodeNormal nn(mesh, type);
  uint const gdim = mesh.geometry().dim();
  ufl::VectorElement space(ufl::Family::CG, mesh.type(), 1, gdim);
  nn.init(mesh, space.repr());
  message("Compute");
  nn.compute();

  std::vector<std::pair<Function *, std::string> > fields;
  for(uint i = 0 ; i < gdim ; ++i)
  {
    std::stringstream ss;
    ss << "E" << i;
    fields.push_back(std::pair<Function *, std::string>(&nn.basis()[i],ss.str()));
  }
  fields.push_back(std::pair<Function *, std::string>(&nn.node_type(),"TYPE"));

  size_t beg = filename.find_last_of('/');
  if(beg != std::string::npos)
  {
    filename.erase(0 ,beg+1);
  }
  size_t pos = filename.find('.');
  message(filename);
  File f(filename.substr(0,pos)+".pvd");
  f << fields;
}

//-----------------------------------------------------------------------------
START_TEST( test_init )
  {
    int init_failed = 0;

    std::string const relpath = "../../../data/meshes/";
    test_node_normal(relpath+"square16.bin", VertexNormal::none);
    test_node_normal(relpath+"cylinder.xml.gz", VertexNormal::none);
    test_node_normal(relpath+"aneurysm.xml.gz", VertexNormal::none);

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
  fprintf(stderr, "*** Check is required for dolfin/fem tests ***\n");
  return 0;
}

#endif
