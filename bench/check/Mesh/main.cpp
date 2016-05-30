#include <dolfin.h>

#include <dolfin/common/AdjacentMapping.h>

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);
  {
    Mesh mesh(t.args.mesh_file);

    if(mesh.is_distributed())
    {

      uint const tdim = mesh.topology().dim();
      mesh.init(tdim - 1);
      SharedMapping sm(mesh.distdata()[tdim - 1]);
    }
  }
  //---------------------------------------------------------------------------
  return 0;
}
