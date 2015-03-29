#ifdef HAVE_CHECK

#include <check.h>

#include <dolfin/config/dolfin_config.h>

#include <dolfin/common/types.h>

using dolfin::uint;
using dolfin::real;

#include <dolfin/ufl/UFLtype.h>

using ufl::type;

//-----------------------------------------------------------------------------
START_TEST( test_UFL_type_int )
{
  int init_failed = 0;
  
  for (uint d = 0; d < 4; ++d)
  {
    type<uint> t(d);
    t.display();
  }
  
  fail_unless( init_failed == 0 );
}END_TEST
 //-----------------------------------------------------------------------------
START_TEST( test_UFL_type_real )
{
  int init_failed = 0;
  
  for (uint d = 0; d < 4; ++d)
  {
    real val = 1.1 * (real) d;
    type<real> t(val);
    t.display();
  }
  
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_UFL_type_string )
{
  int init_failed = 0;
  
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
  
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
