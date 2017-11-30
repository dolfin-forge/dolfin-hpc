#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLDomain.h>

using namespace dolfin;

using ufl::Cell;
using ufl::Domain;

//-----------------------------------------------------------------------------
START_TEST( test_UFL_Cell )
{
  int init_failed = 0;
  
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
  
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
