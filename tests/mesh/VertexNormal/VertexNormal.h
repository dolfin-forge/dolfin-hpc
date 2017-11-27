#include "../../tests.h"

#ifdef HAVE_CHECK

#include <dolfin/mesh/VertexNormal.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void test(std::string filename, VertexNormal::Type type)
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
  
#ifdef HAVE_XML
  test(mesh_file("cylinder.xml.gz"), VertexNormal::none);
  test(mesh_file("aneurysm.xml.gz"), VertexNormal::none);
  test(mesh_file("sphere.xml.gz")  , VertexNormal::none);
#endif
  
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
