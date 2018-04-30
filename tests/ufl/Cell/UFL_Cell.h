#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/ufl/UFLCell.h>

using namespace dolfin;

using ufl::Cell;
using ufl::Domain;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_UFL_Cell )
  {
    Domain::Set domains;
    domains.insert(domains.begin(), Domain::interval);
    domains.insert(domains.begin(), Domain::triangle);
    domains.insert(domains.begin(), Domain::tetrahedron);

    for (ufl::Domain::Set::const_iterator dom_it = domains.begin();
        dom_it != domains.end(); ++dom_it)
    {
      Domain dom(*dom_it);
      Cell cell(dom);
      cell.display();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
