#include <dolfin.h>

using namespace dolfin;

//
uint const DEGMIN = 1;
uint const DEGMAX = 2;
void test(std::string name, NodeNormal& vn)
{
  for (uint i = DEGMIN; i <= DEGMAX; ++i)
  {
    message("Computing for CG%u", i);
    uint const gdim = vn.mesh().geometry().dim();
    ufl::VectorElement cgd(ufl::Family::CG, vn.mesh().type(), i, gdim);
    FiniteElementSpace Vh(vn.mesh(), cgd);
    vn.init(Vh);
    vn.compute();
    std::stringstream ss;
    ss << "CG" << i << name << ".pvd";
    vn.write(ss.str());
  }
}

//
int main(int argc, char** argv)
{
  dolfin_init(argc, argv);
  //---------------------------------------------------------------------------
  {
    {
      Mesh mesh("../../data/meshes/squareN100R.xml.gz");
      NodeNormal vn(mesh, NodeNormal::none);
      test("square_none", vn);
    }

    {
      Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
      NodeNormal vn(mesh, NodeNormal::none);
      test("cube_none", vn);
    }

    {
      Mesh mesh("../../data/meshes/squareN100R.xml.gz");
      NodeNormal vn(mesh, NodeNormal::facet);
      test("square_facet", vn);
    }

    {
      Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
      NodeNormal vn(mesh, NodeNormal::facet);
      test("cube_facet", vn);
    }

  }
  //---------------------------------------------------------------------------
  dolfin_finalize();
  return 0;
}

