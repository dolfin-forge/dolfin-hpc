#include <dolfin/config/dolfin_config.h>

#include <dolfin/common/Array.h>
#include <dolfin/mesh/UnitInterval.h>
#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/mesh/UnitCube.h>
#include <dolfin/mesh/Box.h>
#include <dolfin/mesh/UnitDisk.h>

#include <iostream>
#include <iomanip>
#include <cstdio>

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

int test_mesh(Mesh& mesh)
{
  message("Test mesh");
  mesh.refine();
  mesh.refine();
  return 0;
}

//-----------------------------------------------------------------------------
START_TEST( test_create )
{
  int init_failed = 0;

  dolfin::uint const Nx = 2;
  dolfin::uint const Ny = 4;
  dolfin::uint const Nz = 8;

  // Interval
  begin("Create UnitInterval");
  {
    UnitInterval mesh(Nx);
    init_failed += test_mesh(mesh);
  }
  end();
  skip();

  // Square
  begin("Create UnitSquare");
  {
    Array<UnitSquare::Type> types;
    types.push_back(UnitSquare::right);
    types.push_back(UnitSquare::left);
    types.push_back(UnitSquare::crisscross);
    for(dolfin::uint tp = 0; tp < types.size(); ++tp)
    {
      UnitSquare mesh(Nx, Ny, types[tp]);
      init_failed += test_mesh(mesh);
    }
  }
  end();
  skip();

  // Cube
  begin("Create UnitCube");
  {
    UnitCube mesh(Nx, Ny, Nz);
    init_failed += test_mesh(mesh);
  }
  end();
  skip();

  // Box
  begin("Create Box");
  {
    real a = 0.0;
    real b = 4.0;
    real c = 0.0;
    real d = 2.0;
    real e = 0.0;
    real f = 1.0;
    Box mesh(a,b,c, d, e, f, Nx, Ny, Nz);
    init_failed += test_mesh(mesh);
  }
  end();
  skip();

  // Disk
  begin("Create UnitDisk");
  {
    Array<UnitDisk::Type> types;
    types.push_back(UnitDisk::right);
    types.push_back(UnitDisk::left);
    types.push_back(UnitDisk::crisscross);
    Array<UnitDisk::Transformation> trans;
    trans.push_back(UnitDisk::maxn);
    trans.push_back(UnitDisk::sumn);
    trans.push_back(UnitDisk::rotsumn);
    for(dolfin::uint tp = 0; tp < types.size(); ++tp)
    {
      for(dolfin::uint tr = 0; tr < trans.size(); ++tr)
      {
        UnitDisk mesh(Nx, types[tp], trans[tr]);
        init_failed += test_mesh(mesh);
      }
    }
  }
  end();
  skip();

  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

Suite * test_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("UnitMeshes");
  tc = tcase_create("UnitMeshes");

  tcase_add_test(tc, test_create);

  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc, 60);

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
  fprintf(stderr, "*** Check is required for dolfin/mesh tests ***\n");
  return 0;
}

#endif
