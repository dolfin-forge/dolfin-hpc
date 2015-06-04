#include <dolfin.h>

#include <dolfin/mesh/GlobalFacetMap.h>

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);
  Mesh mesh(t.args.mesh_file);
  bool throw_error = true;

  t.begin_test("Global facet map");
  {
    uint tdim = mesh.topology().dim();
    GlobalFacetMap GFM(mesh);
    mesh.init(tdim - 1, 0);
    for (FacetIterator f(mesh); !f.end(); ++f)
    {
      bool global = GFM.globalFacet(*f);
      if (global && f->is_shared())
      {
        error("A global facet is marked as shared.");
      }

      //
      uint shared_vertex = 0;
      for (VertexIterator v(*f); !v.end(); ++v)
      {
        if(v->is_shared())
        {
          ++shared_vertex;
          _set<uint> const& adjs = mesh.distdata().get_shared_adj(*v);
          if (adjs.empty())
          {
            error("A shared entity should have adjacents.");
          }
        }
      }
    }
  }
  t.end_test();

  return ret;
}

