#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Array.h>
#include <dolfin/mesh/UnitInterval.h>
#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/mesh/UnitCube.h>
#include <dolfin/mesh/Box.h>
#include <dolfin/mesh/UnitDisk.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
int test_mesh(Mesh& mesh)
{
  message("Test mesh");
  mesh.refine();
  mesh.refine();
  return 0;
}
//-----------------------------------------------------------------------------
START_TEST( test_UnitInterval )
{
  dolfin::uint const Nx = 2;

 // Interval
  UnitInterval mesh(Nx);  
  ck_assert(test_mesh(mesh) == 0);

}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_UnitSquare )
{ 
  dolfin::uint const Nx = 2;
  dolfin::uint const Ny = 4;

  // Square
  Array<UnitSquare::Type> types;
  types.push_back(UnitSquare::right);
  types.push_back(UnitSquare::left);
  types.push_back(UnitSquare::crisscross);
  for(dolfin::uint tp = 0; tp < types.size(); ++tp)
  {
    UnitSquare mesh(Nx, Ny, types[tp]);
    ck_assert(test_mesh(mesh) == 0);
  }
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_UnitCube )
{
  dolfin::uint const Nx = 2;
  dolfin::uint const Ny = 4;
  dolfin::uint const Nz = 8;

  // Cube
  UnitCube mesh(Nx, Ny, Nz);
  ck_assert(test_mesh(mesh) == 0);

}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_Box)
{
  // Box
  dolfin::uint const Nx = 2;
  dolfin::uint const Ny = 4;
  dolfin::uint const Nz = 8;

  real a = 0.0;
  real b = 4.0;
  real c = 0.0;
  real d = 2.0;
  real e = 0.0;
  real f = 1.0;
  Box mesh(a,b,c, d, e, f, Nx, Ny, Nz);
  ck_assert(test_mesh(mesh) == 0);


}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_UnitDisk )
{
  // Disk
  dolfin::uint const Nx = 2;
  Array<UnitDisk::Type> types;
  types.push_back(UnitDisk::right);
  types.push_back(UnitDisk::left);
  types.push_back(UnitDisk::crisscross);
  Array<UnitDisk::Transformation> trans;
  trans.push_back(UnitDisk::maxn);
  trans.push_back(UnitDisk::sumn);
  trans.push_back(UnitDisk::rotsumn);
  for(dolfin::uint tp = 0; tp < types.size(); ++tp)
  {
    for(dolfin::uint tr = 0; tr < trans.size(); ++tr)
    {
      UnitDisk mesh(Nx, types[tp], trans[tr]);
      ck_assert(test_mesh(mesh) == 0);
    }
  }
}END_TEST
//-----------------------------------------------------------------------------

#endif
