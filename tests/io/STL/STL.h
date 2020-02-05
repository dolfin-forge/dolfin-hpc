#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/io/STLFile.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/Mesh.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_STLMesh )
{
  Mesh M( path(testdata("io", "STL"), "mug.stl") );
  dolfin_assert(M.is_distributed() == false);
  pushd(testresu("io", "STL"));
  pushd(strcounter(PE::rank()));
  {
    File G("mug.pvd");
    G << M;
  }
  popd();
  popd();
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
