#ifdef HAVE_CHECK

#include <check.h>

#include <dolfin/config/dolfin_config.h>

#include <dolfin/fem/NodeNormal.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/ufl/UFLVectorElement.h>

using namespace dolfin;

void test_node_normal(std::string filename)
{
//  Mesh mesh(filename);
//  NodeNormal nn(mesh);
//  uint const gdim = mesh.geometry().dim();
//  ufl::VectorElement space(ufl::Family::CG, mesh.type(), 1, gdim);
//  nn.init(mesh, space.repr());
//  message("Compute");
//  nn.compute();
//
//  std::vector<std::pair<Function *, std::string> > fields;
//  for(uint i = 0 ; i < gdim ; ++i)
//  {
//    std::stringstream ss;
//    ss << "E" << i;
//    fields.push_back(std::pair<Function *, std::string>(&nn.basis()[i],ss.str()));
//  }
//  fields.push_back(std::pair<Function *, std::string>(&nn.node_type(),"TYPE"));
//
//  size_t beg = filename.find_last_of('/');
//  if(beg != std::string::npos)
//  {
//    filename.erase(0 ,beg+1);
//  }
//  size_t pos = filename.find('.');
//  message(filename);
//  File f(filename.substr(0,pos)+".pvd");
//  f << fields;
}

//-----------------------------------------------------------------------------
START_TEST( test_nodenormal_create )
{
  int init_failed = 0;
  
  std::string const relpath = "../data/meshes/";
  test_node_normal(relpath+"square16.bin");
#ifdef HAVE_XML
  test_node_normal(relpath+"cylinder.xml.gz");
  test_node_normal(relpath+"aneurysm.xml.gz");
  test_node_normal(relpath+"sphere.xml.gz");
#endif
  
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
