#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/ufl/UFLSpace.h>

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_UFL_Space )
  {
    for (dolfin::uint d = 0; d < 4; ++d)
    {
      ufl::Space s(d);
      s.display();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
