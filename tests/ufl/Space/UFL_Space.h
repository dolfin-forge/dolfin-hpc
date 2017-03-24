#ifdef HAVE_CHECK

#include <check.h>

#include <dolfin/config/dolfin_config.h>

#include <dolfin/ufl/UFLSpace.h>

//-----------------------------------------------------------------------------
START_TEST( test_UFL_Space )
{
  int init_failed = 0;
  
  for (dolfin::uint d = 0; d < 4; ++d)
  {
    ufl::Space s(d);
    s.display();
  }
  
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
