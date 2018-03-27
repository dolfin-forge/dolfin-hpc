#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/ufl/UFLDomain.h>

using namespace dolfin;
using ufl::Domain;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_UFL_Domain )
  {
    Domain::Set domains;
    domains.insert(domains.begin(), Domain::interval);
    domains.insert(domains.begin(), Domain::triangle);
    domains.insert(domains.begin(), Domain::tetrahedron);

    for (Domain::Set::const_iterator dom_it = domains.begin();
         dom_it != domains.end(); ++dom_it)
    {
      Domain dom(*dom_it);
      dom.display();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
