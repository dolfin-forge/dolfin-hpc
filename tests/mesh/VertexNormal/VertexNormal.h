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
#ifdef HAVE_XML
    test(mesh_file("cylinder.xml.gz"), VertexNormal::none);
    test(mesh_file("aneurysm.xml.gz"), VertexNormal::none);
    test(mesh_file("sphere.xml.gz"), VertexNormal::none);
#endif
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
  
#endif
