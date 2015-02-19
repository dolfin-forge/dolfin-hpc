#include <dolfin.h>

#include "MeshChecks.h"

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);
  Mesh mesh(t.args.mesh_file);
  bool throw_error = true;

  t.begin_test("Check ghosted and shared entities");
  {
    uint const tdim = mesh.topology().dim();
    for (uint i = 0; i <= tdim; ++i)
    {
      bool ok = true;
      ok &= ghosted_entities_check(mesh, i, throw_error);
      ok &= shared_entities_check(mesh, i, throw_error);
      if (!ok)
      {
        ++ret;
      }
    }
  }
  t.end_test();

  t.begin_test("Check interior boundary entities");
  {
    uint const tdim = mesh.topology().dim();
    for (uint i = 0; i < tdim; ++i)
    {
      bool ok = true;
      ok &= interior_boundary_entities_check(mesh, i, throw_error);
      if (!ok)
      {
        ++ret;
      }
    }
  }
  t.end_test();

  return ret;
}

