#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/VertexNormal.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void test_vertex_normal(std::string filename, VertexNormal::Type type)
{
  Mesh mesh(filename);
  VertexNormal vn(mesh, type);

//  uint const gdim = mesh.geometry().dim();
//  for(uint i = 0 ; i < gdim ; ++i)
//  {
//    std::stringstream ss;
//    ss << "E" << i;
//    size_t pos = filename.find('.');
//    File f(filename.substr(0,pos)+".pvd");
//    f << vn.basis()[d][i];
//  }
}
//-----------------------------------------------------------------------------
START_TEST( test_VertexNormal )
{
  int init_failed = 0;
  
  // NOTE the test is initiated from @top_srcdir@/tests
  std::string const relpath = "../data/meshes/";
#ifdef HAVE_XML
  test_vertex_normal(relpath+"cylinder.xml.gz", VertexNormal::none);
  test_vertex_normal(relpath+"aneurysm.xml.gz", VertexNormal::none);
  test_vertex_normal(relpath+"sphere.xml.gz", VertexNormal::none);
#endif
  
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
