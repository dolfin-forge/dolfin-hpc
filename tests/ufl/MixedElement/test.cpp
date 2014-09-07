#include <dolfin/config/dolfin_config.h>

#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLVectorElement.h>
#include <dolfin/ufl/UFLMixedElement.h>

#include <iostream>
#include <iomanip>

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

    dolfin::UnitSquare mesh(4,4);

    ufl::VectorElement Uspace(ufl::Family::CG, mesh.type(), 2,mesh.geometry().dim());
    Uspace.display();

    ufl::FiniteElement Pspace(ufl::Family::CG, mesh.type(), 1);
    Pspace.display();

    ufl::FiniteElementBase::List spaces;
    spaces.push_back(&Uspace);
    spaces.push_back(&Pspace);

    ufl::MixedElement UPspace(spaces);
    UPspace.display();

    ufl::MixedElement UPspaceFormRepr(UPspace.repr());
    UPspaceFormRepr.display();

    UPspaceFormRepr.sub_elements()[0]->display();
    UPspaceFormRepr.sub_elements()[1]->display();

    if(UPspace != UPspaceFormRepr)
    {
      std::cout << "Representation string differ" << std::endl;
      init_failed = 1;
    }

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------

Suite *ufl_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("UFL");
  tc = tcase_create("ufl");

  tcase_add_test(tc, test_init );

  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);

  return s;
}

int main(void)
{
  int number_failed;
  Suite* s = ufl_suite();
  SRunner* sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;
}

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/ufl tests ***\n");
  return 0;
}

#endif
