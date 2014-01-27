
#include <dolfin/config/dolfin_config.h>

#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLDomain.h>
#include <dolfin/ufl/UFLSpace.h>

using ufl::Cell;
using ufl::Domain;
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

	for(std::vector<Domain::Type>::const_iterator it = domains.begin();
	    it != domains.end(); ++it)
	{
	  Domain d(*it);
	  d.display();

	  Space s(d.dim());
	  s.display();
	}

	fail_unless( init_failed == 0 );
}
END_TEST
//-----------------------------------------------------------------------------

Suite *ufl_suite()
{
	TCase *tc;
	Suite *s;

	s = suite_create("UFL");
	tc = tcase_create("ufl");

	tcase_add_test(tc, test_init_cell);

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

	return (number_failed == 0 ) ? 0 : 1;
}

#else

int main(void)
{
	fprintf(stderr, "*** Check is required for dolfin/ufl tests ***\n");
	return 0;
}

#endif
