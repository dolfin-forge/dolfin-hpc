#ifdef HAVE_CHECK

#include <check.h>

#include <dolfin/config/dolfin_config.h>

#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/UnitSquare.h>

using namespace dolfin;



//bool test_create(SubDomain& sd)
//{
 // uint const dim = sd.mesh().topology().dim();
//  for (uint i = 0; i < dim; ++i)
//  {
//    MeshFunction<uint>& marker = sd.marker(i);
//  }
//  return false;
//}

//-----------------------------------------------------------------------------
START_TEST( test_geometric_create )
{
  int init_failed = 0;

  dolfin::uint const N = 16;
  UnitSquare mesh(N, N);
  //init_failed = test_create(sd);

  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
