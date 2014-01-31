#include <dolfin/config/dolfin_config.h>

#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLDomain.h>
#include <dolfin/ufl/UFLFamily.h>
#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLSpace.h>

using ufl::Cell;
using ufl::CellSurfaceArea;
using ufl::Domain;
using ufl::Family;
using ufl::Object;
using ufl::Space;

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
START_TEST( test_init_cell )
  {
    int init_failed = 0;

    // C++ => Python
    std::vector<Domain::Type> domains;
    domains.push_back(Domain::interval);
    domains.push_back(Domain::triangle);
    domains.push_back(Domain::tetrahedron);

    std::vector<Object::repr_t> domains_repr;
    std::vector<Object::repr_t> spaces_repr;
    std::vector<Object::repr_t> cells_repr;

    std::cout << "======== Create from C++ constructor ========" << std::endl;
    for (std::vector<Domain::Type>::const_iterator it = domains.begin();
        it != domains.end(); ++it)
    {
      Domain d(*it);
      d.display();
      domains_repr.push_back(d.repr());

      Space s(d.dim());
      s.display();
      spaces_repr.push_back(s.repr());

      Cell c(d, s);
      c.display();
      cells_repr.push_back(c.repr());
    }

    // Python => C++
    std::cout << "======== Create from repr constructor ========" << std::endl;
    std::cout << "==== Domains ====" << std::endl;
    for (std::vector<Object::repr_t>::const_iterator it = domains_repr.begin();
        it != domains_repr.end(); ++it)
    {
      Domain d(*it);
      d.display();
      if (*it != d.repr())
      {
        init_failed = 1;
      }
    }
    std::cout << "==== Spaces ====" << std::endl;
    for (std::vector<Object::repr_t>::const_iterator it = spaces_repr.begin();
        it != spaces_repr.end(); ++it)
    {
      Space s(*it);
      s.display();
      if (*it != s.repr())
      {
        init_failed = 1;
      }
    }
    std::cout << "==== Cells ====" << std::endl;
    for (std::vector<Object::repr_t>::const_iterator it = cells_repr.begin();
        it != cells_repr.end(); ++it)
    {
      Cell c(*it);
      c.display();
      if (*it != c.repr())
      {
        init_failed = 1;
      }
    }

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_init_family )
  {
    int init_failed = 0;

    std::vector<Family::Type> v;
    v.push_back(Family::ARG);
    v.push_back(Family::AW);
    v.push_back(Family::BDFM);
    v.push_back(Family::BDM);
    v.push_back(Family::CR);
    v.push_back(Family::DG);
    v.push_back(Family::HER);
    v.push_back(Family::CG);
    v.push_back(Family::MTW);
    v.push_back(Family::MOR);
    v.push_back(Family::N1curl);
    v.push_back(Family::N2curl);
    v.push_back(Family::RT);
    v.push_back(Family::BQ);
    v.push_back(Family::B);
    v.push_back(Family::Q);
    v.push_back(Family::R);
    v.push_back(Family::U);

    std::cout << "======== List families ========" << std::endl;
    for (std::vector<Family::Type>::const_iterator it = v.begin();
        it != v.end(); ++it)
    {
      Family f(*it);
      f.display();
    }
    std::cout << std::endl;

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_init_element_list )
  {
    int init_failed = 0;

    std::cout << "======== Check supported elements ========" << std::endl;
    std::vector<Family::Type> v;
    v.push_back(Family::ARG);
    v.push_back(Family::AW);
    v.push_back(Family::BDFM);
    v.push_back(Family::BDM);
    v.push_back(Family::CR);
    v.push_back(Family::DG);
    v.push_back(Family::HER);
    v.push_back(Family::CG);
    v.push_back(Family::MTW);
    v.push_back(Family::MOR);
    v.push_back(Family::N1curl);
    v.push_back(Family::N2curl);
    v.push_back(Family::RT);
    v.push_back(Family::BQ);
    v.push_back(Family::B);
    v.push_back(Family::Q);
    v.push_back(Family::R);
    v.push_back(Family::U);

    for (std::vector<Family::Type>::const_iterator it = v.begin();
        it != v.end(); ++it)
    {
      bool has_success = Family::has_type(*it);
      init_failed &= !has_success;
      if (has_success)
      {
        Family f(*it);
        std::cout << std::setw(32) << f.name() << " - " << f.short_name()
            << std::endl;
      }
      else
      {
        dolfin::error("Family type does not exist.");
      }
    }
    std::cout << std::endl;

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_init_finite_element )
  {
    int init_failed = 0;

    uint const deg_max = 2;
    std::vector<Family::Type> v;
    v.push_back(Family::DG);
    v.push_back(Family::CG);

    std::cout << "======== Check CG and DG creation (C++) ========"
        << std::endl;
    std::vector<Object::repr_t> fem_repr;
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
        for (uint d = d_min; d <= d_max; ++d)
        {
          ufl::FiniteElement fem(*it, cell, d);
          fem.display();
          fem_repr.push_back(fem.repr());
        }
      }
    }
    std::cout << "======== Check CG and DG creation (repr)========"
        << std::endl;
    for (std::vector<Object::repr_t>::const_iterator it = fem_repr.begin();
        it != fem_repr.end(); ++it)
    {
      ufl::FiniteElement fem(*it);
      fem.display();
      if (*it != fem.repr())
      {
        init_failed = 1;
      }
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

  tcase_add_test(tc, test_init_cell);
  tcase_add_test(tc, test_init_family);
  tcase_add_test(tc, test_init_element_list);
  tcase_add_test(tc, test_init_finite_element);

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
