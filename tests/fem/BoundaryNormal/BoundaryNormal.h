#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/Elements.h>
#include <dolfin/fem/NodeNormal.h>
#include <dolfin/mesh/Mesh.h>

using namespace dolfin;

#include <check.h>

void check_NodeNormal_create(std::string file)
{
  Mesh mesh(file);
  dolfin::uint const gdim = mesh.geometry().dim();
  ufl::VectorElement space(ufl::Family::CG, mesh.type(), 1, gdim);
  FiniteElementSpace Vh(mesh, space);
  NodeNormal nn(mesh);
  nn.init(Vh);
  nn.compute();

  size_t beg = file.find_last_of('/');
  if(beg != std::string::npos)
  {
    file.erase(0 ,beg+1);
  }
  size_t pos = file.find('.');
  nn.write(file.substr(0,pos)+".pvd");
}

//-----------------------------------------------------------------------------
START_TEST( test_NodeNormal )
{
  int init_failed = 0;
  
  std::string const relpath = "../data/meshes/";
#ifdef HAVE_XML
  check_NodeNormal_create(relpath+"cylinder.xml.gz");
  check_NodeNormal_create(relpath+"aneurysm.xml.gz");
  check_NodeNormal_create(relpath+"sphere.xml.gz");
#endif
  
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
