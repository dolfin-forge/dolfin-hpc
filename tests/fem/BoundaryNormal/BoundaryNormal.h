#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/Elements.h>
#include <dolfin/fem/NodeNormal.h>
#include <dolfin/mesh/Mesh.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void test(std::string file)
{
  Mesh mesh(file);
  uint const gdim = mesh.geometry().dim();
  ufl::VectorElement space(ufl::Family::CG, mesh.type(), 1, gdim);
  FiniteElementSpace Vh(mesh, space);
  NodeNormal nn(mesh);
  nn.init(Vh);
  nn.compute();
}
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_NodeNormal )
  {
#ifdef HAVE_XML
    test(mesh_file("cylinder.xml.gz"));
    test(mesh_file("aneurysm.xml.gz"));
    test(mesh_file("sphere.xml.gz"  ));
#endif
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
