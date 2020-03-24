#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/VertexNormal.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void test(std::string filename, VertexNormal::Type type)
{
  Mesh mesh(filename);
  VertexNormal vn(mesh, type);
}
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_VertexNormal )
  {
    test(mesh_file("cylinder.bin"), VertexNormal::none);
    test(mesh_file("aneurysm.bin"), VertexNormal::none);
    test(mesh_file("sphere.bin"), VertexNormal::none);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
