#include <dolfin/config/dolfin_config.h>

#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLDomain.h>
#include <dolfin/ufl/UFLElementList.h>
#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLSpace.h>

using ufl::Cell;
using ufl::CellSurfaceArea;
using ufl::Domain;
using ufl::ElementList;
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

    std::vector<Domain::Type> domains;
    domains.push_back(Domain::interval);
    domains.push_back(Domain::triangle);
    domains.push_back(Domain::tetrahedron);

    for (std::vector<Domain::Type>::const_iterator it = domains.begin();
        it != domains.end(); ++it)
    {
      Domain d(*it);
      d.display();

      Space s(d.dim());
      s.display();

      Cell c(d, s);
      c.display();
    }

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_init_element_list )
  {
    int init_failed = 0;

    // Display full list of supported elements
    ElementList::Supported().display();

    std::vector<ElementList::FamilyType> v;
    v.push_back(ElementList::ARG);
    v.push_back(ElementList::AW);
    v.push_back(ElementList::BDFM);
    v.push_back(ElementList::BDM);
    v.push_back(ElementList::CR);
    v.push_back(ElementList::DG);
    v.push_back(ElementList::HER);
    v.push_back(ElementList::CG);
    v.push_back(ElementList::MTW);
    v.push_back(ElementList::MOR);
    v.push_back(ElementList::N1curl);
    v.push_back(ElementList::N2curl);
    v.push_back(ElementList::RT);
    v.push_back(ElementList::BQ);
    v.push_back(ElementList::B);
    v.push_back(ElementList::Q);
    v.push_back(ElementList::R);
    v.push_back(ElementList::U);

    for (std::vector<ElementList::FamilyType>::const_iterator it = v.begin();
        it != v.end(); ++it)
    {
      bool has_success = ElementList::Supported().has_family(*it);
      init_failed &= !has_success;
      if (has_success)
      {
        std::cout << std::setw(32) << ElementList::Supported().name(*it)
            << " - " << ElementList::Supported().short_name(*it) << std::endl;
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
    std::vector<ElementList::FamilyType> v;
    v.push_back(ElementList::DG);
    v.push_back(ElementList::CG);

    for (std::vector<ElementList::FamilyType>::const_iterator it = v.begin();
        it != v.end(); ++it)
    {
      uint d_min = ElementList::Supported().degree_min(*it);
      uint d_max = std::min(ElementList::Supported().degree_max(*it),
                            std::max(d_min, deg_max));
      ufl::DomainSet domains = ElementList::Supported().domains(*it);

      for (ufl::DomainSet::const_iterator dom_it = domains.begin();
          dom_it != domains.end(); ++dom_it)
      {
        Domain dom(*dom_it);
        Cell cell(dom);
        for (uint d = d_min; d <= d_max; ++d)
        {
          ufl::FiniteElement fem(*it, cell, d);
          fem.display();
        }
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
