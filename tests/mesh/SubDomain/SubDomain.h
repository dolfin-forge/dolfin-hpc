#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/UnitSquare.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_SubDomain )
{
  int init_failed = 0;
  //---
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
