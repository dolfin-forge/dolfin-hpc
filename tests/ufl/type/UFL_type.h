#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/common/types.h>

using dolfin::uint;
using dolfin::real;

#include <dolfin/ufl/UFLtype.h>

using ufl::type;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_UFL_type_int )
{
  for (uint d = 0; d < 4; ++d)
  {
    type<uint> t(d);
    t.display();
  }
}DOLFIN_END_TEST
 //-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_UFL_type_real )
{
  for (uint d = 0; d < 4; ++d)
  {
    real val = 1.1 * (real) d;
    type<real> t(val);
    t.display();
  }
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_UFL_type_string )
{
  std::vector<std::string> vals;
  vals.push_back("'interval'");
  vals.push_back("'triangle'");
  vals.push_back("'tetrahedron'");
  
  for (std::vector<std::string>::const_iterator it = vals.begin();
       it != vals.end(); ++it)
  {
    type<std::string> t(*it);
    t.display();
  }
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
