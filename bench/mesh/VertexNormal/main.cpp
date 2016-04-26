#include <dolfin.h>

#include <dolfin/mesh/VertexNormal.h>

using namespace dolfin;

//
void test(std::string name, VertexNormal& vn)
{
  File vt(name + "_vtype.pvd");
  vt << vn.vertex_type();

  Mesh& mesh = vn.mesh();
  uint gdim = mesh.geometry().dim();
  for (uint e = 0; e < gdim; ++e)
 {
    for (uint d = 0; d < gdim; ++d)
    {
      std::stringstream ssb;
      ssb << "B" << e << d;
      //File ved(name + "_"+ssb.str()+".pvd");
      //ved << vn.basis()[e][d];
    }
    ufl::VectorElement cgd(ufl::Family::CG, mesh.type(), 1, gdim);
    FiniteElementSpace Vh(mesh, cgd);
    Function F(Vh);
    uint ii = 0;
    real * block = F.create_block();
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      for (uint d = 0; d < gdim; ++d)
      {
        for (VertexIterator v(*c); !v.end(); ++v)
        {
          block[ii] = vn.basis()[e][d].get(*v);
          ++ii;
        }
      }
    }
    F.set_block(block);
    delete[] block;
    std::stringstream ssf;
    ssf << "F" << e;
    File ved(name + "_" + ssf.str() + ".pvd");
    ved << F;
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
      VertexNormal vn(mesh, VertexNormal::none);
      test("square_none", vn);
    }

    {
      Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
      VertexNormal vn(mesh, VertexNormal::none);
      test("cube_none", vn);
    }

    {
      Mesh mesh("../../data/meshes/squareN100R.xml.gz");
      VertexNormal vn(mesh, VertexNormal::facet);
      test("square_facet", vn);
    }

    {
      Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
      VertexNormal vn(mesh, VertexNormal::facet);
      test("cube_facet", vn);
    }

  }
  //---------------------------------------------------------------------------
  dolfin_finalize();
  return 0;
}

