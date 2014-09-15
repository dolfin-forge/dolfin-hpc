#include <dolfin/config/dolfin_config.h>

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/UnitCube.h>
#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/ufl/UFLFamily.h>
#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLVectorElement.h>

using dolfin::begin;
using dolfin::end;
using dolfin::message;
using dolfin::skip;
using dolfin::Mesh;
using dolfin::ElementLibrary;

using ufl::Cell;
using ufl::Domain;
using ufl::Family;
using ufl::Object;
using ufl::Space;

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

  uint const deg_max = 2;
  std::vector<Family::Type> v;
  v.push_back(Family::DG);
  v.push_back(Family::CG);

  dolfin::UnitSquare sq(4, 4);
  dolfin::UnitCube cb(2, 2, 2);

  for (std::vector<Family::Type>::const_iterator it = v.begin();
      it != v.end(); ++it)
  {
    Family f(*it);
    uint d_min = f.degree_min();
    uint d_max = std::min(f.degree_max(), std::max(d_min, deg_max));
    ufl::Domain::Set domains = f.domains();

    for (ufl::Domain::Set::const_iterator dom_it = domains.begin();
        dom_it != domains.end(); ++dom_it)
    {
      Domain dom(*dom_it);
      Cell cell(dom);
      uint const dim = cell.topological_dimension();
      dolfin::Mesh * m;
      if(dim == 2)
      {
        m = &sq;
      }
      else if (dim == 3)
      {
        m = &cb;
      }
      else
      {
        continue;
      }
      for (uint d = d_min; d <= d_max; ++d)
      {
        begin("DiscreteSpace");
        //
        skip();

        begin("Creating UFLFiniteElement:");
        ufl::FiniteElement uflfem(*it, cell, d);
        message(uflfem.repr());
        end();
        skip();

        begin("Creating UFLFiniteElement from factory function:");
        ufl::FiniteElementBase * factuflfem = ufl::FiniteElementBase::create(uflfem.repr());
        message(factuflfem->repr());
        delete factuflfem;
        end();
        skip();

        begin("Creating FiniteElement from UFL representation:");
        dolfin::FiniteElement fem(uflfem);
        message(fem.signature());
        if(fem.signature() != uflfem.repr())
        {
          init_failed = 1;
          message("WRONG");
        }
        else
        {
          message("MATCHING");
        }
        skip();
        fem.disp();
        end();
        skip();

        begin("Creating UFL representation from FiniteElement:");
        ufl::FiniteElement uflfemd(Object::repr_t(fem.signature()));
        message(uflfemd.repr());
        if(uflfem.repr() != uflfemd.repr())
        {
          init_failed = 1;
          message("WRONG");
        }
        else
        {
          message("MATCHING");
        }
        end();
        skip();

        begin("Creating corresponding DofMap:");
        ufc::dofmap * ufcdm = ElementLibrary::create_dof_map(dolfin::DofMap::dofmap_signature(fem.signature()));
        dolfin::DofMap dm(*m, *ufcdm, true);
        message(dm.signature());
        if(dm.signature() != dolfin::DofMap::dofmap_signature(uflfem.repr()))
        {
          init_failed = 1;
          message("WRONG");
        }
        else
        {
          message("MATCHING");
        }
        skip();
        dm.disp();
        end();
        skip();

        begin("Creating corresponding DiscreteSpace:");
        dolfin::FiniteElementSpace femspace(*m,fem, dm, false); // Cannot set true
        femspace.disp();
        end();
        skip();

        //
        end();
      }
    }
  }

  fail_unless( init_failed == 0 );
}END_TEST

//-----------------------------------------------------------------------------

Suite * test_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("data");
  tc = tcase_create("data");

  tcase_add_test(tc, test_init );

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
  fprintf(stderr, "*** Check is required for dolfin/data tests ***\n");
  return 0;
}

#endif
